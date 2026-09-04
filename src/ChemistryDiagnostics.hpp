#ifndef CHEMISTRYDIAGNOSTICS_HPP
#define CHEMISTRYDIAGNOSTICS_HPP

#include "Configuration.hpp"
#include "IonizationVariables.hpp"
#include "SafeGslOde.hpp"
#ifdef HAVE_HDF5
#include "HDF5Tools.hpp"
#endif
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <vector>

/** Optional diagnostics: complete cell masks, thread-local event totals and
 * bounded examples. No event-list truncation can bias the spatial masks.
 * Cell order is exactly the original task-subgrid output order. */
class ChemistryDiagnostics {
  bool _enabled;
  std::vector<uint32_t> _restored, _corrected;
  std::vector<std::array<uint64_t, 2*6*5*32>> _counts;
  std::vector<std::vector<std::string>> _samples;
  std::vector<std::array<double, 2>> _worker_seconds;
  std::ofstream _count_file, _sample_file;

public:
  ChemistryDiagnostics(bool enabled, size_t cells, int threads)
      : _enabled(enabled) {
    if (!enabled) return;
    _restored.resize(cells, 0);
    _corrected.resize(cells, 0);
    _counts.resize(threads);
    for (auto &counts : _counts) counts.fill(0);
    _samples.resize(threads);
    _worker_seconds.resize(threads, {{0., 0.}});
    _count_file.open("chemistry_gsl_counts.tsv", std::ios::app);
    _sample_file.open("chemistry_gsl_samples.tsv", std::ios::app);
    if (!_count_file || !_sample_file) cmac_error("Cannot open chemistry diagnostics.");
    _count_file << "# time_s phase(0=hydro,1=trial) element(0=HHe,1=C,2=N,3=O,4=Ne,5=S) kind(0=restore,1=input,2=output_floor,3=simplex,4=attempt) gsl_status events\n";
    _sample_file << "# interval_end_s cell phase element kind status reason nH_m-3 T_K dt_s reached_s accepted_steps initial_ODE_state output_ODE_state initial_all_ions current_all_ions photo_rates_s-1\n";
  }

  class Scope {
    ChemistryDiagnostics &_owner;
    size_t _cell;
    int _thread, _phase;
    const IonizationVariables &_vars;
    double _dt, _jfac;
    SafeGslOde::DiagnosticContext _previous;
    std::chrono::steady_clock::time_point _start;

    static void record(void *data, SafeGslOde::DiagnosticKind kind,
                       const char *reason, int status,
                       const SafeGslOde::DriverState &state) {
      Scope &s = *static_cast<Scope *>(data);
      const int element = std::max(0, std::min(5, state.element));
      const int status_bin = status >= 0 && status < 31 ? status : 31;
      ++s._owner._counts[s._thread][((s._phase*6+element)*5+kind)*32+status_bin];
      if (kind == SafeGslOde::ATTEMPTED) return;
      const uint32_t bit = uint32_t(1) << (element + 8*s._phase);
      if (kind == SafeGslOde::RESTORED) s._owner._restored[s._cell] |= bit;
      else s._owner._corrected[s._cell] |= bit;
      // Routine trace-ion floors should not displace examples of real failures.
      auto &samples = s._owner._samples[s._thread];
      if (kind == SafeGslOde::OUTPUT_CLAMPED || samples.size() >= 32) return;
      std::ostringstream row;
      row.precision(17);
      row << s._cell << '\t' << s._phase << '\t' << element << '\t'
          << kind << '\t' << status << '\t' << reason << '\t'
          << s._vars.get_number_density() << '\t' << s._vars.get_temperature()
          << '\t' << s._dt << '\t' << state.reached_time << '\t' << state.steps;
      for (int which = 0; which < 5; ++which) {
        row << '\t';
        const size_t count = which < 2 ? std::min(state.dimension,
            SafeGslOde::MAXIMUM_TRACKED_DIMENSION) : NUMBER_OF_IONNAMES;
        for (size_t i = 0; i < count; ++i) {
          if (i) row << ',';
          row << (which == 0 ? state.initial[i] : which == 1 ? state.output[i]
              : which == 2 ? s._vars.get_prev_ionic_fraction(i)
              : which == 3 ? s._vars.get_ionic_fraction(i)
                           : s._jfac*s._vars.get_mean_intensity(i));
        }
      }
      samples.push_back(row.str());
    }
  public:
    Scope(ChemistryDiagnostics &owner, size_t cell, int thread,
          const IonizationVariables &vars, double dt, bool trial = false,
          double jfac = 1.)
        : _owner(owner), _cell(cell), _thread(thread), _phase(trial),
          _vars(vars), _dt(dt), _jfac(jfac),
          _previous(SafeGslOde::diagnostic_context()) {
      if (owner._enabled) {
        _start = std::chrono::steady_clock::now();
        SafeGslOde::diagnostic_context().callback = record;
        SafeGslOde::diagnostic_context().data = this;
      }
    }
    ~Scope() {
      if (_owner._enabled) _owner._worker_seconds[_thread][_phase] +=
          std::chrono::duration<double>(std::chrono::steady_clock::now()-_start).count();
      SafeGslOde::diagnostic_context() = _previous;
    }
  };

  // Called serially after the worker region; no locks in the cell/ODE path.
  void flush(double time) {
    if (!_enabled) return;
    for (int phase=0; phase<2; ++phase) {
      double seconds = 0.;
      for (auto &worker : _worker_seconds) { seconds += worker[phase]; worker[phase] = 0.; }
      if (seconds) _count_file << "# worker_seconds " << time << ' ' << phase << ' ' << seconds << '\n';
    }
    for (int phase=0; phase<2; ++phase) for (int element=0; element<6; ++element)
      for (int kind=0; kind<5; ++kind) for (int status=0; status<32; ++status) {
        const int index = ((phase*6+element)*5+kind)*32+status;
        uint64_t total = 0;
        for (auto &counts : _counts) { total += counts[index]; counts[index] = 0; }
        if (total) _count_file << time << '\t' << phase << '\t' << element
            << '\t' << kind << '\t' << status << '\t' << total << '\n';
      }
    for (auto &samples : _samples) {
      for (const auto &row : samples) _sample_file << time << '\t' << row << '\n';
      samples.clear();
    }
    _count_file.flush();
    _sample_file.flush();
  }

  void write_snapshot(const std::string &filename) {
    if (!_enabled || filename.empty()) return;
#ifdef HAVE_HDF5
    auto file = HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_APPEND);
    auto group = HDF5Tools::create_group(file, "ChemistryDiagnostics");
    std::string description = "Since previous output or restart; bits 0..5: hydro HHe,C,N,O,Ne,S; bits 8..13: radiation trials. Original subgrid/cell order.";
    HDF5Tools::write_attribute<std::string>(group, "MaskMeaning", description);
    HDF5Tools::write_dataset<uint32_t>(group, "Restored", _restored);
    HDF5Tools::write_dataset<uint32_t>(group, "Corrected", _corrected);
    HDF5Tools::close_group(group);
    HDF5Tools::close_file(file);
    std::fill(_restored.begin(), _restored.end(), 0);
    std::fill(_corrected.begin(), _corrected.end(), 0);
#endif
  }
};
#endif
