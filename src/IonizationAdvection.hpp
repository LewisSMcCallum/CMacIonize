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
 * @file IonizationAdvection.hpp
 *
 * @brief Bounded multidimensional second-order reconstruction and conservative
 * positivity-preserving flux helpers for advected ionic fractions.
 */
#ifndef IONIZATIONADVECTION_HPP
#define IONIZATIONADVECTION_HPP

#include "CoordinateVector.hpp"
#include "IonizationVariables.hpp"

#include <algorithm>
#include <cmath>

namespace IonizationAdvection {

/** @brief Two-argument minmod limiter. */
inline double minmod(const double a, const double b) {
  if (a * b <= 0.) {
    return 0.;
  }
  return std::copysign(std::min(std::abs(a), std::abs(b)), a);
}

/** @brief Monotonized-central slope on a uniform three-cell stencil. */
inline double mc_slope(const double xm, const double x0, const double xp) {
  const double dl = x0 - xm;
  const double dr = xp - x0;
  return minmod(0.5 * (xp - xm), minmod(2. * dl, 2. * dr));
}

/** @brief Limited one-sided slope on three consecutive cell centres. */
inline double one_sided_slope(const double x0, const double x1,
                              const double x2) {
  return minmod(x1 - x0, x2 - x1);
}

/**
 * @brief Limited slope for one cell and one coordinate direction.
 *
 * Prefer the symmetric MC stencil. At a local subgrid/domain edge, retain a
 * bounded one-sided second-order slope when two cells are available on the
 * interior side. Falling back to zero is deliberately first order but safe.
 */
inline double limited_cell_slope(
    const IonizationVariables &cell, const IonizationVariables *minus,
    const IonizationVariables *plus, const IonizationVariables *minus_two,
    const IonizationVariables *plus_two, const int_fast32_t ion) {
  const double x0 = cell.get_ionic_fraction(ion);
  if (minus != nullptr && plus != nullptr) {
    return mc_slope(minus->get_ionic_fraction(ion), x0,
                    plus->get_ionic_fraction(ion));
  }
  if (minus == nullptr && plus != nullptr && plus_two != nullptr) {
    return one_sided_slope(x0, plus->get_ionic_fraction(ion),
                           plus_two->get_ionic_fraction(ion));
  }
  if (plus == nullptr && minus != nullptr && minus_two != nullptr) {
    return one_sided_slope(minus_two->get_ionic_fraction(ion),
                           minus->get_ionic_fraction(ion), x0);
  }
  return 0.;
}

/** @brief Clamp a fraction to supplied monotonic/physical bounds. */
inline double bound_face_fraction(const double value, const double lower,
                                  const double upper) {
  return std::max(std::max(0., lower),
                  std::min(std::min(1., upper), value));
}

/** @brief Renormalize one element's explicitly stored ion stages if needed. */
inline void normalize_group(double fractions[NUMBER_OF_IONNAMES],
                            const int_fast32_t first,
                            const int_fast32_t last) {
  double sum = 0.;
  for (int_fast32_t ion = first; ion <= last; ++ion) {
    fractions[ion] = std::max(0., std::min(1., fractions[ion]));
    sum += fractions[ion];
  }
  if (sum > 1.) {
    const double inverse_sum = 1. / sum;
    for (int_fast32_t ion = first; ion <= last; ++ion) {
      fractions[ion] *= inverse_sum;
    }
  }
}

/**
 * @brief Enforce the physical simplex for every explicitly stored element.
 *
 * The highest ionization stage is implicit in CMacIonize, so the sum of the
 * stored stages for an element must not exceed unity.
 */
inline void enforce_ionic_simplex(double fractions[NUMBER_OF_IONNAMES]) {
  fractions[ION_H_n] = std::max(0., std::min(1., fractions[ION_H_n]));
#ifdef HAS_HELIUM
  normalize_group(fractions, ION_He_n, ION_He_p1);
#endif
#ifdef HAS_CARBON
  normalize_group(fractions, ION_C_p1, ION_C_p2);
#endif
#ifdef HAS_NITROGEN
  normalize_group(fractions, ION_N_n, ION_N_p2);
#endif
#ifdef HAS_OXYGEN
  normalize_group(fractions, ION_O_n, ION_O_p3);
#endif
#ifdef HAS_NEON
  normalize_group(fractions, ION_Ne_n, ION_Ne_p3);
#endif
#ifdef HAS_SULPHUR
  normalize_group(fractions, ION_S_p1, ION_S_p3);
#endif
#ifdef HAS_ARGON
  normalize_group(fractions, ION_Ar_n, ION_Ar_p3);
#endif
#ifdef HAS_MAGNESIUM
  fractions[ION_Mg_p1] =
      std::max(0., std::min(1., fractions[ION_Mg_p1]));
#endif
}

/** @brief Enforce the ionic simplex on a cell state. */
inline void enforce_ionic_simplex(IonizationVariables &variables) {
  double fractions[NUMBER_OF_IONNAMES];
  for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
    fractions[ion] = variables.get_ionic_fraction(ion);
  }
  enforce_ionic_simplex(fractions);
  for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
    variables.set_ionic_fraction(ion, fractions[ion]);
  }
}

/**
 * @brief Reconstruct one donor-cell ionic state at a face using a genuinely
 * multidimensional MUSCL-Hancock predictor.
 *
 * Slopes are dimensionless cell-to-cell changes. The cell state is first
 * predicted to half time with -0.5 dt v.grad(x), including all three
 * directions, then extrapolated by half a slope to the requested face. The
 * result is bounded by the local multidimensional stencil and physical simplex.
 */
inline void reconstruct_face_3d(
    const IonizationVariables &cell,
    const double slopes[3][NUMBER_OF_IONNAMES],
    const double lower[NUMBER_OF_IONNAMES],
    const double upper[NUMBER_OF_IONNAMES], const CoordinateVector<> &velocity,
    const double cell_size[3], const uint_fast8_t direction,
    const bool high_side, const double dt,
    double face[NUMBER_OF_IONNAMES]) {

  double courant[3];
  for (uint_fast8_t d = 0; d < 3; ++d) {
    courant[d] = cell_size[d] > 0. ? velocity[d] * dt / cell_size[d] : 0.;
    if (!std::isfinite(courant[d])) {
      courant[d] = 0.;
    }
    // A healthy hydro step is already inside its CFL limit. This is only a
    // guard against a damaged velocity feeding an unbounded scalar predictor.
    courant[d] = std::max(-1., std::min(1., courant[d]));
  }

  for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
    double advective_prediction = 0.;
    for (uint_fast8_t d = 0; d < 3; ++d) {
      advective_prediction += courant[d] * slopes[d][ion];
    }
    const double half_time =
        cell.get_ionic_fraction(ion) - 0.5 * advective_prediction;
    const double reconstructed =
        half_time + (high_side ? 0.5 : -0.5) * slopes[direction][ion];
    face[ion] = bound_face_fraction(reconstructed, lower[ion], upper[ion]);
  }
  enforce_ionic_simplex(face);
}

/**
 * @brief Legacy directional interface reconstruction.
 *
 * Retained for the specialised shearing-periodic remap. Ordinary Cartesian
 * hydro faces use reconstruct_face_3d().
 */
inline void reconstruct_interface(
    const IonizationVariables *left_minus, const IonizationVariables &left,
    const IonizationVariables &right, const IonizationVariables *right_plus,
    double courant_left, double courant_right,
    double left_face[NUMBER_OF_IONNAMES],
    double right_face[NUMBER_OF_IONNAMES]) {

  if (!std::isfinite(courant_left)) {
    courant_left = 0.;
  }
  if (!std::isfinite(courant_right)) {
    courant_right = 0.;
  }
  courant_left = std::max(-1., std::min(1., courant_left));
  courant_right = std::max(-1., std::min(1., courant_right));

  for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
    const double xl = left.get_ionic_fraction(ion);
    const double xr = right.get_ionic_fraction(ion);

    double slope_left = 0.;
    if (left_minus != nullptr) {
      slope_left = mc_slope(left_minus->get_ionic_fraction(ion), xl, xr);
    } else if (right_plus != nullptr) {
      slope_left = one_sided_slope(
          xl, xr, right_plus->get_ionic_fraction(ion));
    }

    double slope_right = 0.;
    if (right_plus != nullptr) {
      slope_right = mc_slope(xl, xr, right_plus->get_ionic_fraction(ion));
    } else if (left_minus != nullptr) {
      slope_right = one_sided_slope(
          left_minus->get_ionic_fraction(ion), xl, xr);
    }

    const double reconstructed_left =
        xl + 0.5 * (1. - courant_left) * slope_left;
    const double reconstructed_right =
        xr - 0.5 * (1. + courant_right) * slope_right;

    left_face[ion] = bound_face_fraction(
        reconstructed_left, std::min(xl, xr), std::max(xl, xr));
    right_face[ion] = bound_face_fraction(
        reconstructed_right, std::min(xl, xr), std::max(xl, xr));
  }

  enforce_ionic_simplex(left_face);
  enforce_ionic_simplex(right_face);
}

/**
 * @brief Legacy one-sided physical-boundary reconstruction.
 */
inline void reconstruct_boundary(
    const IonizationVariables &boundary,
    const IonizationVariables *inside_one,
    const IonizationVariables *inside_two, const bool high_side,
    double courant, double face[NUMBER_OF_IONNAMES]) {

  if (!std::isfinite(courant)) {
    courant = 0.;
  }
  courant = std::max(-1., std::min(1., courant));

  for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
    const double x0 = boundary.get_ionic_fraction(ion);
    double slope = 0.;
    double neighbour = x0;
    if (inside_one != nullptr) {
      neighbour = inside_one->get_ionic_fraction(ion);
    }
    if (inside_one != nullptr && inside_two != nullptr) {
      if (high_side) {
        slope = one_sided_slope(inside_two->get_ionic_fraction(ion),
                                inside_one->get_ionic_fraction(ion), x0);
      } else {
        slope = one_sided_slope(x0, inside_one->get_ionic_fraction(ion),
                                inside_two->get_ionic_fraction(ion));
      }
    }

    const double reconstructed =
        high_side ? x0 + 0.5 * (1. - courant) * slope
                  : x0 - 0.5 * (1. + courant) * slope;
    face[ion] = bound_face_fraction(reconstructed, std::min(x0, neighbour),
                                    std::max(x0, neighbour));
  }
  enforce_ionic_simplex(face);
}

/**
 * @brief Bound one element's outgoing face composition by the ion mass that is
 * actually available in the donor cell.
 *
 * The implicit highest stage is included in the limiter. Candidate face
 * fractions are projected onto the simplex subject to per-stage upper bounds
 * q_stage/(|F_rho| dt). Therefore no explicit or implicit stage can export more
 * ion mass than the donor contains after already accumulated conservative
 * fluxes. Any correction is made to the face composition before the equal and
 * opposite flux is deposited, so conservation is retained.
 */
inline void limit_group_outflow(
    const IonizationVariables &donor, const double donor_mass,
    const double donor_mass_delta_before, const double outflow_mass_rate,
    const double dt, const int_fast32_t first, const int_fast32_t last,
    double face[NUMBER_OF_IONNAMES]) {

  if (!(outflow_mass_rate > 0.) || !(dt > 0.) || !(donor_mass > 0.)) {
    return;
  }
  const double transferred_mass = outflow_mass_rate * dt;
  if (!(transferred_mass > 0.) || !std::isfinite(transferred_mass)) {
    return;
  }

  const int_fast32_t nexplicit = last - first + 1;
  double candidate[NUMBER_OF_IONNAMES + 1];
  double upper_bound[NUMBER_OF_IONNAMES + 1];
  double limited[NUMBER_OF_IONNAMES + 1];

  double cell_explicit_sum = 0.;
  double face_explicit_sum = 0.;
  double delta_explicit_sum = 0.;
  for (int_fast32_t k = 0; k < nexplicit; ++k) {
    const int_fast32_t ion = first + k;
    const double cell_fraction =
        std::max(0., std::min(1., donor.get_ionic_fraction(ion)));
    candidate[k] = std::max(0., std::min(1., face[ion]));
    cell_explicit_sum += cell_fraction;
    face_explicit_sum += candidate[k];
    const double delta = donor.get_delta_ionic_fraction(ion);
    delta_explicit_sum += delta;
    const double available =
        std::max(0., donor_mass * cell_fraction + dt * delta);
    upper_bound[k] = std::min(1., available / transferred_mass);
  }

  const double implicit_cell =
      std::max(0., 1. - std::min(1., cell_explicit_sum));
  const double implicit_face =
      std::max(0., 1. - std::min(1., face_explicit_sum));
  const double implicit_delta =
      donor_mass_delta_before - delta_explicit_sum;
  const double implicit_available =
      std::max(0., donor_mass * implicit_cell + dt * implicit_delta);
  candidate[nexplicit] = implicit_face;
  upper_bound[nexplicit] =
      std::min(1., implicit_available / transferred_mass);

  double limited_sum = 0.;
  for (int_fast32_t k = 0; k <= nexplicit; ++k) {
    limited[k] = std::min(candidate[k], upper_bound[k]);
    limited_sum += limited[k];
  }

  // Clipping a stage to its available budget creates a deficit in the face
  // simplex. Redistribute that deficit only to stages with remaining physical
  // capacity. This is a bounded projection and changes nothing when the
  // high-order flux was already positivity preserving.
  double deficit = std::max(0., 1. - limited_sum);
  if (deficit > 0.) {
    double spare_sum = 0.;
    for (int_fast32_t k = 0; k <= nexplicit; ++k) {
      spare_sum += std::max(0., upper_bound[k] - limited[k]);
    }
    if (spare_sum > 0.) {
      const double amount = std::min(deficit, spare_sum);
      for (int_fast32_t k = 0; k <= nexplicit; ++k) {
        const double spare = std::max(0., upper_bound[k] - limited[k]);
        limited[k] += amount * spare / spare_sum;
      }
      deficit -= amount;
    }
  }

  // A remaining deficit would mean that this one hydro face exports more total
  // mass than the donor provisionally owns. There is then no composition-only
  // limiter that can make every stage positive. The hydro CFL normally makes
  // this impossible; retain a physical simplex as a fail-safe and let the
  // existing cell-level guard catch any pathological hydro state.
  if (deficit > 1.e-12) {
    double sum_available = 0.;
    for (int_fast32_t k = 0; k <= nexplicit; ++k) {
      sum_available += upper_bound[k];
    }
    if (sum_available > 0.) {
      for (int_fast32_t k = 0; k <= nexplicit; ++k) {
        limited[k] = upper_bound[k] / sum_available;
      }
    }
  }

  for (int_fast32_t k = 0; k < nexplicit; ++k) {
    face[first + k] = std::max(0., std::min(1., limited[k]));
  }
}

/** @brief Apply the conservative positivity limiter to every element. */
inline void limit_outflow_face(
    const IonizationVariables &donor, const double donor_mass,
    const double donor_mass_delta_before, const double outflow_mass_rate,
    const double dt, double face[NUMBER_OF_IONNAMES]) {
  limit_group_outflow(donor, donor_mass, donor_mass_delta_before,
                      outflow_mass_rate, dt, ION_H_n, ION_H_n, face);
#ifdef HAS_HELIUM
  limit_group_outflow(donor, donor_mass, donor_mass_delta_before,
                      outflow_mass_rate, dt, ION_He_n, ION_He_p1, face);
#endif
#ifdef HAS_CARBON
  limit_group_outflow(donor, donor_mass, donor_mass_delta_before,
                      outflow_mass_rate, dt, ION_C_p1, ION_C_p2, face);
#endif
#ifdef HAS_NITROGEN
  limit_group_outflow(donor, donor_mass, donor_mass_delta_before,
                      outflow_mass_rate, dt, ION_N_n, ION_N_p2, face);
#endif
#ifdef HAS_OXYGEN
  limit_group_outflow(donor, donor_mass, donor_mass_delta_before,
                      outflow_mass_rate, dt, ION_O_n, ION_O_p3, face);
#endif
#ifdef HAS_NEON
  limit_group_outflow(donor, donor_mass, donor_mass_delta_before,
                      outflow_mass_rate, dt, ION_Ne_n, ION_Ne_p3, face);
#endif
#ifdef HAS_SULPHUR
  limit_group_outflow(donor, donor_mass, donor_mass_delta_before,
                      outflow_mass_rate, dt, ION_S_p1, ION_S_p3, face);
#endif
#ifdef HAS_ARGON
  limit_group_outflow(donor, donor_mass, donor_mass_delta_before,
                      outflow_mass_rate, dt, ION_Ar_n, ION_Ar_p3, face);
#endif
#ifdef HAS_MAGNESIUM
  limit_group_outflow(donor, donor_mass, donor_mass_delta_before,
                      outflow_mass_rate, dt, ION_Mg_p1, ION_Mg_p1, face);
#endif
  enforce_ionic_simplex(face);
}

/** @brief Add a conservative ion-mass flux to a pair of cells. */
inline void add_interface_flux(
    IonizationVariables &left, IonizationVariables &right,
    const double mflux, const double left_face[NUMBER_OF_IONNAMES],
    const double right_face[NUMBER_OF_IONNAMES]) {
  if (mflux == 0.) {
    return;
  }
  const double *upwind = mflux > 0. ? left_face : right_face;
  for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
    const double ion_flux = mflux * upwind[ion];
    left.increase_delta_ionic_fraction(ion, -ion_flux);
    right.increase_delta_ionic_fraction(ion, ion_flux);
  }
}

/** @brief Add one already-limited donor face flux to a pair of cells. */
inline void add_donor_interface_flux(
    IonizationVariables &left, IonizationVariables &right, const double mflux,
    const double donor_face[NUMBER_OF_IONNAMES]) {
  if (mflux == 0.) {
    return;
  }
  for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
    const double ion_flux = mflux * donor_face[ion];
    left.increase_delta_ionic_fraction(ion, -ion_flux);
    right.increase_delta_ionic_fraction(ion, ion_flux);
  }
}

/** @brief Add the ionic flux through a physical boundary to its interior cell. */
inline void add_boundary_flux(
    IonizationVariables &cell, const double mflux,
    const double face[NUMBER_OF_IONNAMES]) {
  for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
    cell.increase_delta_ionic_fraction(ion, -mflux * face[ion]);
  }
}

} // namespace IonizationAdvection

#endif // IONIZATIONADVECTION_HPP
