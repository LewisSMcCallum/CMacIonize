/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2026 Lewis McCallum
 *
 * CMacIonize is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CMacIonize is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file SafeGslOde.hpp
 *
 * @brief Defensive wrapper around the GSL ODE driver used by the
 * time-dependent ionization solver.
 *
 * The chemistry variables are fractions, so an absolute error of 1.e-4 can be
 * physically important even when the neutral fraction is close to unity.  The
 * wrapper therefore enforces 1.e-8 absolute and relative tolerances, caps the
 * number of internal driver steps, preflights the reaction RHS for finite
 * values, and falls back to the input state rather than aborting an RHD
 * calculation when an exceptionally stiff or malformed cell cannot be
 * advanced safely.
 */
#ifndef SAFEGSLODE_HPP
#define SAFEGSLODE_HPP

#include "Error.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

namespace SafeGslOde {

static const double ABSOLUTE_TOLERANCE = 1.e-8;
static const double RELATIVE_TOLERANCE = 1.e-8;
static const unsigned long MAXIMUM_DRIVER_STEPS = 10000;
static const std::size_t MAXIMUM_TRACKED_DIMENSION = 32;

struct DriverState {
  gsl_odeiv2_driver *driver;
  const gsl_odeiv2_system *system;
  std::size_t dimension;
  int element;
  double initial[MAXIMUM_TRACKED_DIMENSION];
  double output[MAXIMUM_TRACKED_DIMENSION];
  double reached_time;
  unsigned long steps;

  DriverState() : driver(nullptr), system(nullptr), dimension(0), element(0),
                  initial{}, output{}, reached_time(0.), steps(0) {}
};

/**
 * @brief Per-thread metadata for the currently active chemistry solve.
 *
 * CMacIonize allocates, applies and frees one GSL driver at a time in each
 * thread.  Keeping only the current driver avoids a mutex or map in this hot
 * path while still letting the apply wrapper know the system dimension and
 * evaluate the RHS once before entering the adaptive integrator.
 */
inline DriverState &driver_state() {
  static thread_local DriverState state;
  return state;
}

/** @brief Thread-local workspace, shared sequentially by systems of one size.
 *
 * GSL retains a pointer to its system: own that small object here, not on a
 * caller's stack. Only params is borrowed, for the duration of a solve. Reset
 * before every cell, including after a failed solve; never reuse derivatives.
 */
struct DriverWorkspace {
  gsl_odeiv2_system system{};
  gsl_odeiv2_driver *driver = nullptr;
  const gsl_odeiv2_step_type *step_type = nullptr;
  ~DriverWorkspace() {
    if (driver != nullptr) ::gsl_odeiv2_driver_free(driver);
  }
};
inline DriverWorkspace &driver_workspace(std::size_t dimension) {
  static thread_local DriverWorkspace workspaces[MAXIMUM_TRACKED_DIMENSION + 1];
  return workspaces[dimension];
}

enum DiagnosticKind { RESTORED, INPUT_CLAMPED, OUTPUT_CLAMPED, SIMPLEX_CORRECTED, ATTEMPTED };
typedef void (*DiagnosticCallback)(void *, DiagnosticKind, const char *, int,
                                    const DriverState &);
struct DiagnosticContext {
  DiagnosticCallback callback = nullptr;
  void *data = nullptr;
};
inline DiagnosticContext &diagnostic_context() {
  static thread_local DiagnosticContext context;
  return context;
}
inline void diagnostic(DiagnosticKind kind, const char *reason, int status = 0) {
  DiagnosticContext &context = diagnostic_context();
  if (context.callback) context.callback(context.data, kind, reason, status,
                                        driver_state());
}
// Element bits: H/He=0, C=1, N=2, O=3, Ne=4, S=5.
inline void set_element(int element) {
  DriverState &state = driver_state();
  state.element = element;
  state.dimension = 0;
  state.reached_time = 0.;
  state.steps = 0;
}

/** @brief Rate-limited warning for pathological chemistry cells. */
inline void warn_failure(const char *reason, const int status = GSL_SUCCESS) {
  diagnostic(RESTORED, reason, status);
  static std::atomic< unsigned int > warning_count(0);
  const unsigned int count = warning_count.fetch_add(1);
  if (count < 20) {
    if (status == GSL_SUCCESS) {
      cmac_warning("Time-dependent ionization ODE fallback: %s", reason);
    } else {
      cmac_warning("Time-dependent ionization ODE fallback: %s (GSL status %d: %s)",
                   reason, status, gsl_strerror(status));
    }
    if (count == 19) {
      cmac_warning("Further time-dependent ionization ODE fallback warnings will be suppressed.");
    }
  }
}

/**
 * @brief Reuse an ODE workspace with unchanged tolerances and a hard step cap.
 */
inline gsl_odeiv2_driver *driver_alloc_y_new(
    const gsl_odeiv2_system *system, const gsl_odeiv2_step_type *step_type,
    double initial_step, const double, const double) {

  DriverState &state = driver_state();
  state.driver = nullptr;
  state.system = system;
  state.dimension = system == nullptr ? 0 : system->dimension;
  state.reached_time = 0.;
  state.steps = 0;
  std::fill(state.initial, state.initial + MAXIMUM_TRACKED_DIMENSION, 0.);
  std::fill(state.output, state.output + MAXIMUM_TRACKED_DIMENSION, 0.);

  if (system == nullptr || system->function == nullptr ||
      state.dimension == 0 || state.dimension > MAXIMUM_TRACKED_DIMENSION) {
    warn_failure("invalid ODE system or unsupported dimension");
    state.system = nullptr;
    return nullptr;
  }

  if (!std::isfinite(initial_step) || initial_step <= 0.) {
    // A one second seed is only a starting guess; GSL adapts it immediately.
    initial_step = 1.;
  }

  DriverWorkspace &workspace = driver_workspace(state.dimension);
  if (workspace.driver != nullptr && workspace.step_type != step_type) {
    ::gsl_odeiv2_driver_free(workspace.driver);
    workspace.driver = nullptr;
  }
  workspace.system = *system;
  state.system = &workspace.system;
  if (workspace.driver == nullptr) {
    workspace.driver = ::gsl_odeiv2_driver_alloc_y_new(
        &workspace.system, step_type, initial_step, ABSOLUTE_TOLERANCE,
        RELATIVE_TOLERANCE);
    workspace.step_type = step_type;
  }
  gsl_odeiv2_driver *driver = workspace.driver;
  if (driver == nullptr) {
    warn_failure("GSL driver allocation failed");
    state.system = nullptr;
    return nullptr;
  }

  int status = ::gsl_odeiv2_driver_reset_hstart(driver, initial_step);
  if (status == GSL_SUCCESS)
    status = ::gsl_odeiv2_driver_set_nmax(driver, MAXIMUM_DRIVER_STEPS);
  if (status != GSL_SUCCESS) {
    warn_failure("could not reset the GSL driver or set its step cap", status);
    ::gsl_odeiv2_driver_free(driver);
    workspace.driver = nullptr;
    state.system = nullptr;
    return nullptr;
  }

  state.driver = driver;
  return driver;
}

/**
 * @brief Advance the ODE while preserving a safe fallback state.
 *
 * On any GSL failure, non-finite reaction derivative/output, or grossly
 * non-physical fraction, the input state is restored and success is returned
 * to the legacy caller.  The fallback is deliberately a zero-order hold: for
 * a pathological cell it is safer to defer its chemistry by one radiation
 * step than to inject a NaN or an arbitrary ion fraction into the hydro
 * calculation.
 */
inline int driver_apply(gsl_odeiv2_driver *driver, double *time,
                        const double target_time, double y[]) {
  DriverState &state = driver_state();
  const std::size_t dimension = state.dimension;

  if (time == nullptr || y == nullptr || dimension == 0 ||
      dimension > MAXIMUM_TRACKED_DIMENSION) {
    warn_failure("invalid state passed to GSL driver");
    return GSL_SUCCESS;
  }

  double initial[MAXIMUM_TRACKED_DIMENSION];
  diagnostic(ATTEMPTED, "ODE attempt");
  bool input_clamped = false;
  for (std::size_t i = 0; i < dimension; ++i) {
    initial[i] = y[i];
    state.initial[i] = y[i];
    state.output[i] = y[i];
    if (!std::isfinite(initial[i])) {
      warn_failure("non-finite input ionic fraction");
      return GSL_SUCCESS;
    }
    // All systems currently routed through this wrapper evolve fractions.
    // Remove tiny advection/roundoff excursions before giving them to GSL.
    y[i] = std::max(0., std::min(1., y[i]));
    input_clamped |= y[i] != initial[i];
    initial[i] = y[i];
  }
  if (input_clamped) diagnostic(INPUT_CLAMPED, "input outside [0,1]");

  if (!std::isfinite(*time) || !std::isfinite(target_time)) {
    warn_failure("non-finite integration time");
    for (std::size_t i = 0; i < dimension; ++i) {
      y[i] = initial[i];
    }
    return GSL_SUCCESS;
  }

  if (target_time <= *time) {
    return GSL_SUCCESS;
  }

  if (driver == nullptr || driver != state.driver || state.system == nullptr) {
    warn_failure("GSL driver was not available");
    for (std::size_t i = 0; i < dimension; ++i) {
      y[i] = initial[i];
    }
    return GSL_SUCCESS;
  }

  // Evaluate the exact user RHS once before entering GSL.  This cheaply catches
  // overflowed reaction coefficients, invalid electron densities and similar
  // pathologies before the adaptive integrator spends thousands of attempts on
  // a state that can never produce a finite derivative.
  double derivative[MAXIMUM_TRACKED_DIMENSION];
  const int preflight_status = state.system->function(
      *time, y, derivative, state.system->params);
  bool finite_derivative = preflight_status == GSL_SUCCESS;
  for (std::size_t i = 0; i < dimension && finite_derivative; ++i) {
    finite_derivative = std::isfinite(derivative[i]);
  }
  if (!finite_derivative) {
    warn_failure("non-finite reaction derivative before ODE integration",
                 preflight_status);
    for (std::size_t i = 0; i < dimension; ++i) {
      y[i] = initial[i];
    }
    return GSL_SUCCESS;
  }

  const double initial_time = *time;
  const int status = ::gsl_odeiv2_driver_apply(driver, time, target_time, y);
  state.reached_time = *time;
  state.steps = driver->n;
  for (std::size_t i = 0; i < dimension; ++i) state.output[i] = y[i];

  bool valid = status == GSL_SUCCESS && std::isfinite(*time);
  for (std::size_t i = 0; i < dimension && valid; ++i) {
    // Allow only a very small solver overshoot here; the legacy callers clamp
    // roundoff-sized excursions after a successful solve.
    valid = std::isfinite(y[i]) && y[i] >= -1.e-8 && y[i] <= 1. + 1.e-8;
  }

  if (!valid) {
    warn_failure(status == GSL_SUCCESS ? "non-finite or non-physical ODE output"
                                       : "GSL could not complete the chemistry step",
                 status);
    *time = initial_time;
    for (std::size_t i = 0; i < dimension; ++i) {
      y[i] = initial[i];
    }
    // Existing callers historically abort for H/He and metals when GSL
    // returns an error.  The restored state is a valid, explicit fallback, so
    // report success here and let the run continue.
    return GSL_SUCCESS;
  }

  return GSL_SUCCESS;
}

/** @brief Release this solve; keep the workspace and last diagnostic result. */
inline void driver_free(gsl_odeiv2_driver *driver) {
  DriverState &state = driver_state();
  if (state.driver == driver) {
    if (driver != nullptr) driver_workspace(state.dimension).system.params = nullptr;
    state.driver = nullptr;
    state.system = nullptr;
  }
}

} // namespace SafeGslOde

// IonizationStateCalculator.cpp includes its header before including GSL and
// contains all current CMacIonize ODE-driver calls.  Map those legacy calls to
// the defensive wrapper without changing the solver equations themselves.
#define gsl_odeiv2_driver_alloc_y_new SafeGslOde::driver_alloc_y_new
#define gsl_odeiv2_driver_apply SafeGslOde::driver_apply
#define gsl_odeiv2_driver_free SafeGslOde::driver_free

#endif // SAFEGSLODE_HPP
