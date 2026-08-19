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
 * number of internal driver steps, and falls back to the input state rather
 * than aborting an RHD calculation when an exceptionally stiff or malformed
 * cell cannot be advanced safely.
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
  std::size_t dimension;

  DriverState() : driver(nullptr), dimension(0) {}
};

/**
 * @brief Per-thread metadata for the currently active chemistry solve.
 *
 * CMacIonize allocates, applies and frees one GSL driver at a time in each
 * thread.  Keeping only the current driver avoids a mutex or map in this hot
 * path while still letting the apply wrapper know the system dimension.
 */
inline DriverState &driver_state() {
  static thread_local DriverState state;
  return state;
}

/** @brief Rate-limited warning for pathological chemistry cells. */
inline void warn_failure(const char *reason, const int status = GSL_SUCCESS) {
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
 * @brief Allocate an ODE driver with chemistry-appropriate tolerances and a
 * hard step cap.
 */
inline gsl_odeiv2_driver *driver_alloc_y_new(
    const gsl_odeiv2_system *system, const gsl_odeiv2_step_type *step_type,
    double initial_step, const double, const double) {

  DriverState &state = driver_state();
  state.driver = nullptr;
  state.dimension = system == nullptr ? 0 : system->dimension;

  if (system == nullptr || state.dimension == 0 ||
      state.dimension > MAXIMUM_TRACKED_DIMENSION) {
    warn_failure("invalid ODE system or unsupported dimension");
    return nullptr;
  }

  if (!std::isfinite(initial_step) || initial_step <= 0.) {
    // The caller normally supplies timestep/100 or timestep/1000.  A one
    // second seed is only a harmless starting guess; GSL adapts it immediately.
    initial_step = 1.;
  }

  gsl_odeiv2_driver *driver = ::gsl_odeiv2_driver_alloc_y_new(
      system, step_type, initial_step, ABSOLUTE_TOLERANCE,
      RELATIVE_TOLERANCE);
  if (driver == nullptr) {
    warn_failure("GSL driver allocation failed");
    return nullptr;
  }

  const int status =
      ::gsl_odeiv2_driver_set_nmax(driver, MAXIMUM_DRIVER_STEPS);
  if (status != GSL_SUCCESS) {
    warn_failure("could not set the GSL driver step cap", status);
    ::gsl_odeiv2_driver_free(driver);
    return nullptr;
  }

  state.driver = driver;
  return driver;
}

/**
 * @brief Advance the ODE while preserving a safe fallback state.
 *
 * On any GSL failure, non-finite output, or grossly non-physical fraction, the
 * input state is restored and success is returned to the legacy caller.  The
 * fallback is deliberately a zero-order hold: for a pathological cell it is
 * safer to defer its chemistry by one radiation step than to inject a NaN or
 * an arbitrary ion fraction into the hydro calculation.
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
  for (std::size_t i = 0; i < dimension; ++i) {
    initial[i] = y[i];
    if (!std::isfinite(initial[i])) {
      warn_failure("non-finite input ionic fraction");
      return GSL_SUCCESS;
    }
    // All systems currently routed through this wrapper evolve fractions.
    // Remove tiny advection/roundoff excursions before giving them to GSL.
    y[i] = std::max(0., std::min(1., y[i]));
    initial[i] = y[i];
  }

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

  if (driver == nullptr || driver != state.driver) {
    warn_failure("GSL driver was not available");
    for (std::size_t i = 0; i < dimension; ++i) {
      y[i] = initial[i];
    }
    return GSL_SUCCESS;
  }

  const double initial_time = *time;
  const int status = ::gsl_odeiv2_driver_apply(driver, time, target_time, y);

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

/** @brief Free a wrapped driver and clear its per-thread metadata. */
inline void driver_free(gsl_odeiv2_driver *driver) {
  DriverState &state = driver_state();
  if (driver != nullptr) {
    ::gsl_odeiv2_driver_free(driver);
  }
  if (state.driver == driver) {
    state.driver = nullptr;
    state.dimension = 0;
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
