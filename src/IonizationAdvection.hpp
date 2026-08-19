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
 * @brief Bounded second-order reconstruction helpers for advected ionic
 * fractions.
 */
#ifndef IONIZATIONADVECTION_HPP
#define IONIZATIONADVECTION_HPP

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

/**
 * @brief Limited one-sided slope, used only where a four-cell stencil is not
 * available inside a subgrid.
 */
inline double one_sided_slope(const double x0, const double x1,
                              const double x2) {
  return minmod(x1 - x0, x2 - x1);
}

/** @brief Clamp a reconstructed fraction to the local interface bounds. */
inline double bound_face_fraction(const double value, const double left,
                                  const double right) {
  const double lower = std::max(0., std::min(left, right));
  const double upper = std::min(1., std::max(left, right));
  return std::max(lower, std::min(upper, value));
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
 * @brief Reconstruct the two ionic states at an interface.
 *
 * This is a directionally time-centred MUSCL-Hancock reconstruction for the
 * passive mass fractions. A monotonized-central slope is used when the full
 * four-cell stencil is available. At a subgrid edge a limited one-sided slope
 * retains second-order accuracy in smooth flow without storing three gradients
 * for every ion in every cell. Reconstructed states are kept inside the local
 * interface extrema and the physical ionic simplex.
 *
 * @param left_minus Cell immediately before left, or nullptr at a local edge.
 * @param left Left cell.
 * @param right Right cell.
 * @param right_plus Cell immediately after right, or nullptr at a local edge.
 * @param courant_left v_n dt/dx for the left cell.
 * @param courant_right v_n dt/dx for the right cell.
 * @param left_face Output state reconstructed from the left.
 * @param right_face Output state reconstructed from the right.
 */
inline void reconstruct_interface(
    const IonizationVariables *left_minus, const IonizationVariables &left,
    const IonizationVariables &right, const IonizationVariables *right_plus,
    double courant_left, double courant_right,
    double left_face[NUMBER_OF_IONNAMES],
    double right_face[NUMBER_OF_IONNAMES]) {

  // A healthy hydro step already satisfies this through its CFL condition.
  // The clamp prevents a damaged velocity from amplifying a scalar slope.
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

    left_face[ion] = bound_face_fraction(reconstructed_left, xl, xr);
    right_face[ion] = bound_face_fraction(reconstructed_right, xl, xr);
  }

  enforce_ionic_simplex(left_face);
  enforce_ionic_simplex(right_face);
}

/**
 * @brief Reconstruct an ionic state at a physical domain boundary from the
 * three interior cells nearest that boundary.
 *
 * @param boundary Boundary cell.
 * @param inside_one First cell moving inward from the boundary, if available.
 * @param inside_two Second cell moving inward, if available.
 * @param high_side True for the positive-coordinate boundary.
 * @param courant v_n dt/dx in the positive coordinate direction.
 * @param face Output reconstructed boundary state.
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
    face[ion] = bound_face_fraction(reconstructed, x0, neighbour);
  }
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
