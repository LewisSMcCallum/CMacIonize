/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2026 Lewis McCallum
 *
 * CMacIonize is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 ******************************************************************************/

/**
 * @file GalacticShearingBox.hpp
 *
 * @brief Local Galactic rotation source terms for a co-rotating Cartesian box.
 */
#ifndef GALACTICSHEARINGBOX_HPP
#define GALACTICSHEARINGBOX_HPP

#include "DensitySubGridCreator.hpp"
#include "Hydro.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

#include <cmath>

/**
 * @brief Coriolis and radial tidal terms in the local Galactic frame.
 *
 * CMacIonize uses x toward the Galactic centre and y along rotation. This is
 * the opposite radial sign to the common outward-x shearing-sheet convention:
 *
 *   du_x/dt = -2 Omega u_y + 2 q Omega^2 x
 *   du_y/dt =  2 Omega u_x
 *
 * The source update is an exact rotation about the local background shear
 * u_y = q Omega x while x is held fixed during this operator-split step.
 */
class GalacticShearingBox {
private:
  const bool _enabled;
  const bool _initialize_background_shear;
  const double _omega;
  const double _shear_parameter;
  const double _radial_centre;

public:
  GalacticShearingBox(ParameterFile &params, Log *log = nullptr)
      : _enabled(params.get_value< bool >("GalacticShearingBox:enabled", false)),
        _initialize_background_shear(params.get_value< bool >(
            "GalacticShearingBox:initialize background shear", true)),
        _omega(params.get_physical_value< QUANTITY_VELOCITY >(
                   "GalacticShearingBox:circular velocity", "232. km s^-1") /
               params.get_physical_value< QUANTITY_LENGTH >(
                   "GalacticShearingBox:solar radius", "8.2 kpc")),
        _shear_parameter(params.get_value< double >(
            "GalacticShearingBox:shear parameter", 1.)),
        _radial_centre(params.get_physical_value< QUANTITY_LENGTH >(
            "GalacticShearingBox:radial centre", "0. pc")) {
    if (_enabled && log) {
      log->write_status("Enabled local Galactic rotation with Omega = ", _omega,
                        " s^-1 and q = ", _shear_parameter, ".");
      log->write_warning(
          "GalacticShearingBox currently supplies the local source terms and "
          "initial linear shear, but x boundaries are ordinary periodic "
          "boundaries without the time-dependent y remap. Treat this as a "
          "development setup, not a complete shearing-periodic box.");
    }
  }

  inline bool enabled() const { return _enabled; }

  /**
   * @brief Apply the Galactic-frame velocity update to a collisionless source.
   *
   * This is the same exact Coriolis/tidal rotation used for the gas.  The
   * position is held fixed during this operator-split source update.
   */
  inline void apply_to_source(const CoordinateVector<> &position,
                              CoordinateVector<> &velocity,
                              const double timestep) const {
    if (!_enabled || timestep <= 0.) {
      return;
    }
    const double angle = 2. * _omega * timestep;
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    const double equilibrium_y =
        _shear_parameter * _omega * (position.x() - _radial_centre);
    const double residual_y = velocity.y() - equilibrium_y;

    const double old_x = velocity.x();
    velocity[0] = old_x * cosine - residual_y * sine;
    velocity[1] = equilibrium_y + old_x * sine + residual_y * cosine;
  }

  /** @brief Add the equilibrium linear shear to the initial velocity field. */
  inline void initialize(
      DensitySubGridCreator< HydroDensitySubGrid > &grid_creator) const {
    if (!_enabled || !_initialize_background_shear) {
      return;
    }
    for (auto grid = grid_creator.begin(); grid != grid_creator.original_end();
         ++grid) {
      for (auto cell = (*grid).hydro_begin(); cell != (*grid).hydro_end();
           ++cell) {
        CoordinateVector<> velocity =
            cell.get_hydro_variables().get_primitives_velocity();
        velocity[1] += _shear_parameter * _omega *
                       (cell.get_cell_midpoint().x() - _radial_centre);
        cell.get_hydro_variables().set_primitives_velocity(velocity);
      }
    }
  }

  /**
   * @brief Apply one operator-split Coriolis/tidal source update.
   *
   * Momentum and kinetic energy are changed together, leaving internal energy
   * unchanged by the frame transformation.
   */
  inline void apply(HydroDensitySubGrid &subgrid, const Hydro &hydro,
                    const double timestep) const {
    if (!_enabled || timestep <= 0.) {
      return;
    }
    const double angle = 2. * _omega * timestep;
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    for (auto cell = subgrid.hydro_begin(); cell != subgrid.hydro_end(); ++cell) {
      HydroVariables &variables = cell.get_hydro_variables();
      const double mass = variables.get_conserved_mass();
      if (mass <= 0.) {
        continue;
      }
      const CoordinateVector<> old_velocity =
          variables.get_conserved_momentum() / mass;
      const double equilibrium_y =
          _shear_parameter * _omega *
          (cell.get_cell_midpoint().x() - _radial_centre);
      const double residual_y = old_velocity.y() - equilibrium_y;

      CoordinateVector<> new_velocity = old_velocity;
      new_velocity[0] = old_velocity.x() * cosine - residual_y * sine;
      new_velocity[1] = equilibrium_y + old_velocity.x() * sine +
                        residual_y * cosine;

      const double kinetic_change =
          0.5 * mass * (new_velocity.norm2() - old_velocity.norm2());
      variables.set_conserved_momentum(mass * new_velocity);
      variables.set_conserved_total_energy(
          variables.get_conserved_total_energy() + kinetic_change);
      hydro.set_primitive_variables(variables,
                                    cell.get_ionization_variables(),
                                    1. / cell.get_volume());
    }
  }
};

#endif // GALACTICSHEARINGBOX_HPP
