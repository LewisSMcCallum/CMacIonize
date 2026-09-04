/** Small chemistry regression/benchmark; never launches a simulation. */
#include "Assert.hpp"
#include "Abundances.hpp"
#include "ChargeTransferRates.hpp"
#include "ChiantiRecombinationRates.hpp"
#include "CollisionalRates.hpp"
#include "HydroStepChemistry.hpp"
#include <chrono>
#include <iomanip>
#include <limits>

// Compare the wrapper against fresh, unwrapped GSL drivers.
#undef gsl_odeiv2_driver_alloc_y_new
#undef gsl_odeiv2_driver_apply
#undef gsl_odeiv2_driver_free

struct Reaction { double rate, equilibrium; size_t dimension; };
static int rhs(double, const double y[], double f[], void *data) {
  const Reaction &r = *static_cast<Reaction *>(data);
  for (size_t i = 0; i < r.dimension; ++i)
    f[i] = r.rate * (r.equilibrium - y[i]);
  return GSL_SUCCESS;
}

static void reuse() {
  // Changing dimension/parameters and revisiting earlier systems must not
  // retain another cell's derivative or integration history, even on threads.
#ifdef HAVE_OPENMP
#pragma omp parallel for num_threads(4)
#endif
  for (int sample = 0; sample < 200; ++sample) {
    Reaction r{.1 + sample % 7, .01 * (sample % 90), 1 + size_t(sample % 4)};
    gsl_odeiv2_system system{rhs, nullptr, r.dimension, &r};
    double reference[4] = {.1, .2, .3, .4}, actual[4] = {.1, .2, .3, .4};
    double t = 0.;
    auto fresh = gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45,
        .01, SafeGslOde::ABSOLUTE_TOLERANCE, SafeGslOde::RELATIVE_TOLERANCE);
    assert_condition(gsl_odeiv2_driver_apply(fresh, &t, 1., reference) == GSL_SUCCESS);
    gsl_odeiv2_driver_free(fresh);
    auto cached = SafeGslOde::driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45, .01, 0., 0.);
    t = 0.;
    SafeGslOde::driver_apply(cached, &t, 1., actual);
    for (size_t i = 0; i < r.dimension; ++i) {
      assert_condition(reference[i] == actual[i]);
      const double exact = r.equilibrium + (.1 * (i+1) - r.equilibrium) * std::exp(-r.rate);
      assert_values_equal_tol(actual[i], exact, 1.e-7);
    }
    SafeGslOde::driver_free(cached);
    // A bad cell must restore, and must not poison the next use of its slot.
    r.rate = std::numeric_limits<double>::quiet_NaN();
    cached = SafeGslOde::driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45, .01, 0., 0.);
    t = 0.; actual[0] = .3;
    SafeGslOde::driver_apply(cached, &t, 1., actual);
    assert_condition(t == 0. && actual[0] == .3);
    SafeGslOde::driver_free(cached);
  }
  Reaction r{1.e6, .5, 1};
  gsl_odeiv2_system system{rhs, nullptr, 1, &r};
  auto cached = SafeGslOde::driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45, 1., 0., 0.);
  // Force exhaustion cheaply, then reuse the same workspace successfully.
  gsl_odeiv2_driver_set_nmax(cached, 2);
  double t = 0., y[] = {.1};
  SafeGslOde::driver_apply(cached, &t, 1., y);
  assert_condition(t == 0. && y[0] == .1);
  SafeGslOde::driver_free(cached);
  r.rate = .1;
  cached = SafeGslOde::driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45, 1., 0., 0.);
  assert_condition(cached->nmax == SafeGslOde::MAXIMUM_DRIVER_STEPS);
  SafeGslOde::driver_apply(cached, &t, 1., y);
  assert_condition(t == 1.);
  assert_values_equal_tol(y[0], .5 - .4 * std::exp(-.1), 1.e-7);
  SafeGslOde::driver_free(cached);
}

static IonizationVariables cell(double n, double T, double j, double fraction) {
  IonizationVariables value;
  value.set_number_density(n);
  value.set_temperature(T);
  for (int ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
    value.set_ionic_fraction(ion, .01);
    value.set_mean_intensity(ion, j);
  }
  value.set_ionic_fraction(ION_H_n, fraction);
#ifdef HAS_HELIUM
  value.set_ionic_fraction(ION_He_n, fraction);
  value.set_ionic_fraction(ION_He_p1, .5 * (1. - fraction));
#endif
  return value;
}

int main() {
  reuse();
  Abundances abundances(.1, 1.4e-4, 7.5e-5, 3.19e-4, 1.17e-4, 1.86e-5);
  ChiantiRecombinationRates rr;
  ChargeTransferRates ctr;
  CollisionalRates cr;
  IonizationStateCalculator calculator(1., abundances, rr, ctr, cr);
  size_t restores = 0;
  SafeGslOde::diagnostic_context().data = &restores;
  SafeGslOde::diagnostic_context().callback = [](void *data,
      SafeGslOde::DiagnosticKind kind, const char *, int,
      const SafeGslOde::DriverState &) {
    if (kind == SafeGslOde::RESTORED) ++*static_cast<size_t *>(data);
  };
  // Actual atomic rates: cold, warm, ionized and SN-hot cells, at both hydro
  // and radiation timesteps. Printed rows also support before/after comparison.
  std::cout << std::setprecision(17);
  size_t row = 0;
  for (double n : {1.e4, 1.e6, 1.e8})
    for (double T : {100., 8000., 50000., 1.e6})
      for (double j : {0., 1.e-12})
        for (double fraction : {1.e-6, .999})
          for (double dt : {3.e9, 1.5e12}) {
            auto sample = cell(n, T, j, fraction);
            restores = 0;
            HydroStepChemistry::advance(sample, calculator, dt);
            std::cout << "state " << row++ << ' ' << restores;
            for (int ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
              const double x = sample.get_ionic_fraction(ion);
              assert_condition(std::isfinite(x) && x >= 0. && x <= 1.);
              std::cout << ' ' << x;
            }
            std::cout << '\n';
          }
  // A many-cell, short-hydro-step benchmark, not a full-run speed prediction.
  for (double T : {100., 8000., 50000., 1.e6}) {
    const auto initial = cell(1.e6, T, 1.e-12, .1);
    restores = 0;
    const auto start = std::chrono::steady_clock::now();
    double checksum = 0.;
    for (int repeat = 0; repeat < 5000; ++repeat) {
      auto sample = initial;
      HydroStepChemistry::advance(sample, calculator, 3.e9);
      checksum += sample.get_ionic_fraction(ION_H_n);
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "benchmark " << T << ' ' << seconds << ' ' << restores << ' ' << checksum << '\n';
  }
  SafeGslOde::diagnostic_context() = SafeGslOde::DiagnosticContext();
  std::cout << "Chemistry solver reuse: PASS\n";
}
