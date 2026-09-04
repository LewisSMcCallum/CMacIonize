/** Focused regression checks; no hydrodynamic simulation is required. */
#include "Assert.hpp"
#include "Abundances.hpp"
#include "ChargeTransferRates.hpp"
#include "ChemistryDiagnostics.hpp"
#include "ChiantiRecombinationRates.hpp"
#include "CollisionalRates.hpp"
#include "FrequencyImportanceSampling.hpp"
#include "HydroStepChemistry.hpp"
#include "Pegase3PhotonSourceSpectrum.hpp"
#include "WMBasicPhotonSourceSpectrum.hpp"
#include <chrono>
#include <limits>

static void sampling() {
  // Exactly 1% of the photon luminosity is in the upper half of this spectrum.
  const std::vector<double> nu{1., 2., 3.}, cdf{0., .99, 1.};
  RandomGenerator rng(12345);
  double total = 0., hard = 0.;
  size_t hard_packets = 0;
  const size_t samples = 300000;
  for (size_t i = 0; i < samples; ++i) {
    double weight;
    const double frequency = sample_frequency_importance(nu, cdf, rng, .5, weight);
    total += weight;
    if (frequency >= 2.) { hard += weight; ++hard_packets; }
  }
  assert_condition(std::abs(total/samples - 1.) < .005);
  assert_condition(std::abs(hard/samples - .01) < .00015);
  assert_condition(hard_packets > samples/5);

  // Check the real spectrum, including exact old RNG behaviour when disabled.
  WMBasicPhotonSourceSpectrum spectrum(40000., 25.);
  RandomGenerator old_rng(42), disabled_rng(42), weighted_rng(234);
  double natural_hard = 0., weighted_hard = 0., weighted_total = 0.;
  size_t weighted_hard_packets = 0;
  const double threshold = 35.12 * 1.602176634e-19 /
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK);
  for (size_t i = 0; i < samples; ++i) {
    double weight;
    const double old_nu = spectrum.get_random_frequency(old_rng);
    assert_condition(old_nu == spectrum.get_random_frequency_weighted(disabled_rng, 0., weight));
    assert_condition(weight == 1.);
    natural_hard += old_nu > threshold;
    const double new_nu = spectrum.get_random_frequency_weighted(weighted_rng, .5, weight);
    weighted_total += weight;
    if (new_nu > threshold) { weighted_hard += weight; ++weighted_hard_packets; }
  }
  assert_condition(std::abs(weighted_total/samples - 1.) < .01);
  assert_condition(std::abs(weighted_hard-natural_hard)/samples < .0015);
  assert_condition(weighted_hard_packets > 2*natural_hard);
  std::cout << "WMBasic >35.12 eV: natural packets=" << natural_hard
            << ", importance packets=" << weighted_hard_packets
            << ", weighted equivalent=" << weighted_hard << '\n';
  Pegase3PhotonSourceSpectrum holmes(1.e10, .02);
  double holmes_weight = 0., natural_mean = 0., weighted_mean = 0.;
  for (size_t i=0; i<samples; ++i) {
    double weight;
    natural_mean += holmes.get_random_frequency(old_rng);
    const double frequency = holmes.get_random_frequency_weighted(weighted_rng, .5, weight);
    weighted_mean += frequency*weight;
    holmes_weight += weight;
  }
  assert_condition(std::abs(holmes_weight/samples - 1.) < .01);
  assert_condition(std::abs(weighted_mean/natural_mean - 1.) < .01);
}

static IonizationVariables initial_cell() {
  IonizationVariables cell;
  cell.set_number_density(1.e6);
  cell.set_temperature(8000.);
  for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
    cell.set_ionic_fraction(i, .05);
    cell.set_mean_intensity(i, 1.e-12);
  }
  cell.set_ionic_fraction(ION_H_n, .8);
#ifdef HAS_HELIUM
  cell.set_ionic_fraction(ION_He_n, .8);
  cell.set_ionic_fraction(ION_He_p1, .1);
#endif
  cell.copy_previous_fractions();
  return cell;
}

static void chemistry() {
  Abundances abundances(.1, 1.4e-4, 7.5e-5, 3.19e-4, 1.17e-4, 1.86e-5);
  ChiantiRecombinationRates rr;
  ChargeTransferRates ctr;
  CollisionalRates cr;
  IonizationStateCalculator calculator(1., abundances, rr, ctr, cr);
  const double dt = 1.e10;
  const IonizationVariables initial = initial_cell();
  IonizationVariables cell = initial;
  // Repeated radiation trials must neither evolve the physical initial state
  // repeatedly nor scale already normalized rates again.
  for (int iteration = 0; iteration < 10; ++iteration) {
    for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) cell.set_mean_intensity(i, 2.e-12);
    cell.set_heating(HEATINGTERM_H, 3.);
    HydroStepChemistry::radiation_trial(cell, calculator, .5, dt, iteration == 9);
  }
  for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
    assert_condition(cell.get_ionic_fraction(i) == initial.get_ionic_fraction(i));
    assert_condition(cell.get_mean_intensity(i) == 1.e-12);
  }
  const double heating = 1.5 * PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK);
  assert_values_equal_rel(cell.get_heating(HEATINGTERM_H), heating, 1.e-14);
  IonizationVariables expected = cell;
  calculator.calculate_ionization_state(1., 1., expected, dt, true, true);
  HydroStepChemistry::advance(cell, calculator, dt);
  for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
    assert_condition(cell.get_ionic_fraction(i) == expected.get_ionic_fraction(i));
    assert_condition(cell.get_mean_intensity(i) == 1.e-12);
  }
  assert_values_equal_rel(cell.get_heating(HEATINGTERM_H), heating, 1.e-14);
  HydroStepChemistry::advance(cell, calculator, dt);
  assert_condition(cell.get_ionic_fraction(ION_H_n) < expected.get_ionic_fraction(ION_H_n));

  // Hydrogen recombination has an exact solution: x_H+ = 1/(1+alpha*n*t).
  const double alpha = 2.6e-19, n = 1.e6, t = 1.e12;
  const double neutral = IonizationStateCalculator::compute_time_dependent_hydrogen(alpha, 0., n, 0., 0., t);
  assert_values_equal_tol(neutral, 1.-1./(1.+alpha*n*t), 1.e-7);
  double partial = 0.;
  for (int i=0; i<20; ++i)
    partial = IonizationStateCalculator::compute_time_dependent_hydrogen(alpha, 0., n, 0., partial, t/20);
  assert_values_equal_tol(partial, neutral, 1.e-7);

#ifdef HAS_HELIUM
  double h0 = 0., he0 = 0., hep = 0.;
  IonizationStateCalculator::compute_time_dependent_hydrogen_helium(
      alpha, alpha, alpha, 0., 0., n, .1, 10000., h0, he0, hep, 0., 0., 0., dt);
  assert_condition(std::isfinite(h0) && h0 > 0. && std::isfinite(he0));
  assert_condition(he0 + hep <= 1. + 1.e-15);
  h0 = 1.; he0 = 1.; hep = 1.e-14;
  IonizationStateCalculator::compute_time_dependent_hydrogen_helium(
      alpha, alpha, alpha, 0., 0., n, .1, 100., h0, he0, hep, 0., 0., 0., dt);
  assert_condition(he0 + hep <= 1. + 1.e-15);
#endif

  // A small cost comparison at fixed total elapsed physical time (not a promise
  // about full-run runtime). Short calls repeat allocations/rate evaluations.
  for (double temperature : {8000., 1.e6}) for (int steps : {1, 20}) {
    size_t restores = 0;
    SafeGslOde::diagnostic_context().data = &restores;
    SafeGslOde::diagnostic_context().callback = [](void *data,
        SafeGslOde::DiagnosticKind kind, const char *, int,
        const SafeGslOde::DriverState &) {
      if (kind == SafeGslOde::RESTORED) ++*static_cast<size_t *>(data);
    };
    const auto begin = std::chrono::steady_clock::now();
    for (int repeat=0; repeat<10; ++repeat) {
      IonizationVariables sample = initial;
      sample.set_temperature(temperature);
      for (int step=0; step<steps; ++step)
        HydroStepChemistry::advance(sample, calculator, 1.57788e12/steps);
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    SafeGslOde::diagnostic_context() = SafeGslOde::DiagnosticContext();
    std::cout << "Chemistry T=" << temperature << " steps=" << steps
              << " seconds/10 cells=" << seconds << " restored solves=" << restores << '\n';
  }
}

static int bad_rhs(double, const double[], double f[], void *) {
  f[0] = std::numeric_limits<double>::quiet_NaN();
  return GSL_SUCCESS;
}

static void diagnostics() {
  // Exceed the example cap: every cell must still appear in the full mask.
  ChemistryDiagnostics diagnostics(true, 100, 1);
  const auto cell = initial_cell();
  for (size_t i=0; i<100; ++i) {
    ChemistryDiagnostics::Scope scope(diagnostics, i, 0, cell, 1.);
    SafeGslOde::set_element(0);
    gsl_odeiv2_system system = {bad_rhs, nullptr, 1, nullptr};
    auto driver = gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45, .01, 0., 0.);
    double time = 0., y[] = {.25};
    gsl_odeiv2_driver_apply(driver, &time, 1., y);
    assert_condition(y[0] == .25 && time == 0.);
    gsl_odeiv2_driver_free(driver);
  }
  diagnostics.flush(1.);
#ifdef HAVE_HDF5
  auto file = HDF5Tools::open_file("test_chemistry_masks.hdf5", HDF5Tools::HDF5FILEMODE_WRITE);
  HDF5Tools::close_file(file);
  diagnostics.write_snapshot("test_chemistry_masks.hdf5");
  file = HDF5Tools::open_file("test_chemistry_masks.hdf5", HDF5Tools::HDF5FILEMODE_READ);
  auto group = HDF5Tools::open_group(file, "ChemistryDiagnostics");
  const auto mask = HDF5Tools::read_dataset<uint32_t>(group, "Restored");
  assert_condition(mask.size() == 100);
  for (auto value : mask) assert_condition(value == 1);
  HDF5Tools::close_group(group);
  HDF5Tools::close_file(file);
#endif
}

int main() {
  sampling();
  chemistry();
  diagnostics();
  std::cout << "Sampling, hydro-step chemistry and complete diagnostic masks: PASS\n";
}
