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
#include "IonizationAdvection.hpp"
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
  const bool _shearing_periodic_boundaries;
  const double _omega;
  const double _shear_parameter;
  const double _radial_centre;

  inline static int_fast32_t wrap_index(const int_fast32_t index,
                                        const int_fast32_t size) {
    const int_fast32_t remainder = index % size;
    return remainder < 0 ? remainder + size : remainder;
  }

  /**
   * @brief Locate a Cartesian cell using whole-grid integer coordinates.
   */
  inline static HydroDensitySubGrid::hydroiterator get_cell(
      DensitySubGridCreator< HydroDensitySubGrid > &grid_creator,
      const int_fast32_t ix, const int_fast32_t iy, const int_fast32_t iz) {
    const CoordinateVector< int_fast32_t > subgrid_cells =
        grid_creator.get_subgrid_cell_layout();
    const CoordinateVector< int_fast32_t > subgrid(
        ix / subgrid_cells.x(), iy / subgrid_cells.y(),
        iz / subgrid_cells.z());
    const CoordinateVector< int_fast32_t > local(
        ix % subgrid_cells.x(), iy % subgrid_cells.y(),
        iz % subgrid_cells.z());
    const CoordinateVector< int_fast32_t > layout =
        grid_creator.get_subgrid_layout();
    const size_t subgrid_index =
        subgrid.x() * layout.y() * layout.z() +
        subgrid.y() * layout.z() + subgrid.z();
    HydroDensitySubGrid &grid = *grid_creator.get_subgrid(subgrid_index);
    return grid.hydro_begin() +
           local.x() * subgrid_cells.y() * subgrid_cells.z() +
           local.y() * subgrid_cells.z() + local.z();
  }

  /**
   * @brief Get the two periodic donor rows and their upper-row weight.
   */
  inline static void get_y_weights(const double donor_row,
                                   const int_fast32_t number_of_rows,
                                   int_fast32_t &row0, int_fast32_t &row1,
                                   double &weight1) {
    const double lower = std::floor(donor_row);
    weight1 = donor_row - lower;
    row0 = wrap_index(static_cast< int_fast32_t >(lower), number_of_rows);
    row1 = wrap_index(row0 + 1, number_of_rows);
  }

  inline static HydroVariables
  interpolate_hydro(const HydroVariables &state0,
                    const HydroVariables &state1, const double weight1) {
    const double weight0 = 1. - weight1;
    HydroVariables state;
    for (uint_fast8_t i = 0; i < 5; ++i) {
      state.primitives(i) =
          weight0 * state0.primitives(i) + weight1 * state1.primitives(i);
      state.conserved(i) =
          weight0 * state0.conserved(i) + weight1 * state1.conserved(i);
      state.primitive_gradients(i) =
          weight0 * state0.primitive_gradients(i) +
          weight1 * state1.primitive_gradients(i);
    }
    return state;
  }

  inline static IonizationVariables interpolate_ionization(
      const IonizationVariables &state0, const IonizationVariables &state1,
      const double weight1) {
    const double weight0 = 1. - weight1;
    IonizationVariables state;
    for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
      state.set_ionic_fraction(
          ion, weight0 * state0.get_ionic_fraction(ion) +
                   weight1 * state1.get_ionic_fraction(ion));
    }
    return state;
  }

  /**
   * @brief Express a state in a frame whose y velocity is higher by delta_v.
   */
  inline static void shift_velocity_frame(HydroVariables &state,
                                          const double delta_v) {
    state.primitives(2) += delta_v;
    const double mass = state.conserved(0);
    const double old_y_momentum = state.conserved(2);
    state.conserved(2) = old_y_momentum + delta_v * mass;
    state.conserved(4) +=
        delta_v * old_y_momentum + 0.5 * delta_v * delta_v * mass;
  }

  inline double wrapped_shift(const Box<> &box, const double time) const {
    const double y_size = box.get_sides().y();
    double shift =
        std::fmod(_shear_parameter * _omega * box.get_sides().x() * time,
                  y_size);
    if (shift < 0.) {
      shift += y_size;
    }
    return shift;
  }

public:
  GalacticShearingBox(ParameterFile &params, Log *log = nullptr)
      : _enabled(params.get_value< bool >("GalacticShearingBox:enabled", false)),
        _initialize_background_shear(params.get_value< bool >(
            "GalacticShearingBox:initialize background shear", true)),
        _shearing_periodic_boundaries(params.get_value< bool >(
            "GalacticShearingBox:shearing periodic boundaries", false)),
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
      if (_shearing_periodic_boundaries) {
        log->write_status(
            "Enabled conservative hydrodynamic shearing-periodic x "
            "boundaries.");
      } else {
        log->write_warning(
            "GalacticShearingBox x boundaries are ordinary periodic "
            "boundaries. Set \"shearing periodic boundaries: true\" for the "
            "time-dependent hydrodynamic y remap.");
      }
    }
  }

  inline bool enabled() const { return _enabled; }

  inline bool shearing_periodic_boundaries() const {
    return _enabled && _shearing_periodic_boundaries;
  }

  /**
   * @brief Check the geometrical requirements of the shearing remap.
   */
  inline void validate_boundaries(
      const CoordinateVector< bool > &periodicity,
      DensitySubGridCreator< HydroDensitySubGrid > &grid_creator) const {
    if (!shearing_periodic_boundaries()) {
      return;
    }
    if (grid_creator.is_spherical() || !periodicity.x() || !periodicity.y()) {
      cmac_error("Galactic shearing-periodic boundaries require a Cartesian "
                 "grid with periodic x and y.");
    }
    const CoordinateVector< int_fast32_t > layout =
        grid_creator.get_subgrid_layout();
    const size_t high_x_grid =
        (layout.x() - 1) * layout.y() * layout.z();
    const size_t high_y_grid = (layout.y() - 1) * layout.z();
    if ((*grid_creator.get_subgrid(high_x_grid))
                .get_neighbour(TRAVELDIRECTION_FACE_X_P) ==
            NEIGHBOUR_OUTSIDE ||
        (*grid_creator.get_subgrid(high_y_grid))
                .get_neighbour(TRAVELDIRECTION_FACE_Y_P) ==
            NEIGHBOUR_OUTSIDE) {
      cmac_error("Galactic shearing-periodic boundaries also require "
                 "DensitySubGridCreator periodicity in x and y.");
    }
  }

  /**
   * @brief Is this the ordinary periodic task replaced by the shear remap?
   */
  inline bool replaces_task(
      const size_t subgrid, const int_fast32_t direction,
      const DensitySubGridCreator< HydroDensitySubGrid > &grid_creator) const {
    return shearing_periodic_boundaries() &&
           direction == TRAVELDIRECTION_FACE_X_P &&
           grid_creator.get_grid_position(subgrid).x() ==
               grid_creator.get_subgrid_layout().x() - 1;
  }

  /**
   * @brief Add the remapped x-boundary contribution to cell gradients.
   *
   * The low-x state sampled by a high-x cell at y is taken from
   * y-q*Omega*Lx*t and shifted into the high-x velocity frame. The inverse
   * mapping is applied independently at x low.
   */
  inline void add_boundary_gradients(
      DensitySubGridCreator< HydroDensitySubGrid > &grid_creator,
      const Hydro &hydro, const double time) const {
    if (!shearing_periodic_boundaries()) {
      return;
    }
    const CoordinateVector< int_fast32_t > layout =
        grid_creator.get_subgrid_layout();
    const CoordinateVector< int_fast32_t > subgrid_cells =
        grid_creator.get_subgrid_cell_layout();
    const int_fast32_t nx = layout.x() * subgrid_cells.x();
    const int_fast32_t ny = layout.y() * subgrid_cells.y();
    const int_fast32_t nz = layout.z() * subgrid_cells.z();
    const Box<> box = grid_creator.get_box();
    const double shift_rows =
        wrapped_shift(box, time) * ny / box.get_sides().y();
    const double velocity_jump =
        _shear_parameter * _omega * box.get_sides().x();

    for (int_fast32_t iy = 0; iy < ny; ++iy) {
      int_fast32_t low0, low1, high0, high1;
      double low_weight1, high_weight1;
      get_y_weights(iy - shift_rows, ny, low0, low1, low_weight1);
      get_y_weights(iy + shift_rows, ny, high0, high1, high_weight1);
      for (int_fast32_t iz = 0; iz < nz; ++iz) {
        auto high_cell = get_cell(grid_creator, nx - 1, iy, iz);
        const auto low_cell0 = get_cell(grid_creator, 0, low0, iz);
        const auto low_cell1 = get_cell(grid_creator, 0, low1, iz);
        HydroVariables low_ghost =
            interpolate_hydro(low_cell0.get_hydro_variables(),
                              low_cell1.get_hydro_variables(), low_weight1);
        shift_velocity_frame(low_ghost, velocity_jump);
        HydroDensitySubGrid &high_grid =
            *grid_creator.get_subgrid(high_cell.get_cell_midpoint());
        high_grid.add_x_remapped_ghost_gradient(
            high_cell.get_index(), hydro, low_ghost, true);

        auto low_cell = get_cell(grid_creator, 0, iy, iz);
        const auto high_cell0 = get_cell(grid_creator, nx - 1, high0, iz);
        const auto high_cell1 = get_cell(grid_creator, nx - 1, high1, iz);
        HydroVariables high_ghost =
            interpolate_hydro(high_cell0.get_hydro_variables(),
                              high_cell1.get_hydro_variables(), high_weight1);
        shift_velocity_frame(high_ghost, -velocity_jump);
        HydroDensitySubGrid &low_grid =
            *grid_creator.get_subgrid(low_cell.get_cell_midpoint());
        low_grid.add_x_remapped_ghost_gradient(
            low_cell.get_index(), hydro, high_ghost, false);
      }
    }
  }

  /**
   * @brief Add conservative, linearly remapped fluxes at the radial boundary.
   *
   * One flux is computed per high-x face and deposited back onto the two
   * overlapping low-x cells. Mass and the non-shearing momentum components
   * are therefore conserved to round-off. The y momentum and total energy are
   * transformed between the two local background-shear frames.
   */
  inline void add_boundary_fluxes(
      DensitySubGridCreator< HydroDensitySubGrid > &grid_creator,
      const Hydro &hydro, const double time, const double timestep,
      const bool advect_ionization) const {
    if (!shearing_periodic_boundaries()) {
      return;
    }
    const CoordinateVector< int_fast32_t > layout =
        grid_creator.get_subgrid_layout();
    const CoordinateVector< int_fast32_t > subgrid_cells =
        grid_creator.get_subgrid_cell_layout();
    const int_fast32_t nx = layout.x() * subgrid_cells.x();
    const int_fast32_t ny = layout.y() * subgrid_cells.y();
    const int_fast32_t nz = layout.z() * subgrid_cells.z();
    const Box<> box = grid_creator.get_box();
    const double dx = box.get_sides().x() / nx;
    const double area =
        (box.get_sides().y() / ny) * (box.get_sides().z() / nz);
    const double shift_rows =
        wrapped_shift(box, time) * ny / box.get_sides().y();
    const double velocity_jump =
        _shear_parameter * _omega * box.get_sides().x();

    for (int_fast32_t iy = 0; iy < ny; ++iy) {
      int_fast32_t low0, low1;
      double weight1;
      get_y_weights(iy - shift_rows, ny, low0, low1, weight1);
      const double weight0 = 1. - weight1;
      for (int_fast32_t iz = 0; iz < nz; ++iz) {
        auto high_cell = get_cell(grid_creator, nx - 1, iy, iz);
        auto low_cell0 = get_cell(grid_creator, 0, low0, iz);
        auto low_cell1 = get_cell(grid_creator, 0, low1, iz);

        HydroVariables high_state;
        high_state.copy_all(high_cell.get_hydro_variables());
        HydroVariables low_state =
            interpolate_hydro(low_cell0.get_hydro_variables(),
                              low_cell1.get_hydro_variables(), weight1);
        for (uint_fast8_t i = 0; i < 5; ++i) {
          high_state.delta_conserved(i) = 0.;
          low_state.delta_conserved(i) = 0.;
        }
        shift_velocity_frame(low_state, velocity_jump);

        IonizationVariables high_ionization;
        IonizationVariables low_ionization = interpolate_ionization(
            low_cell0.get_ionization_variables(),
            low_cell1.get_ionization_variables(), weight1);
        for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
          high_ionization.set_ionic_fraction(
              ion, high_cell.get_ionization_variables().get_ionic_fraction(
                       ion));
        }

        double high_face[NUMBER_OF_IONNAMES];
        double low_face[NUMBER_OF_IONNAMES];
        if (advect_ionization) {
          const IonizationVariables *high_minus = nullptr;
          const IonizationVariables *low_plus = nullptr;
          IonizationVariables low_plus_state;
          if (nx > 1) {
            auto high_inside = get_cell(grid_creator, nx - 2, iy, iz);
            high_minus = &high_inside.get_ionization_variables();
            auto low_plus0 = get_cell(grid_creator, 1, low0, iz);
            auto low_plus1 = get_cell(grid_creator, 1, low1, iz);
            low_plus_state = interpolate_ionization(
                low_plus0.get_ionization_variables(),
                low_plus1.get_ionization_variables(), weight1);
            low_plus = &low_plus_state;
          }
          const double high_courant =
              high_state.get_primitives_velocity().x() * timestep / dx;
          const double low_courant =
              low_state.get_primitives_velocity().x() * timestep / dx;
          IonizationAdvection::reconstruct_interface(
              high_minus, high_ionization, low_ionization, low_plus,
              high_courant, low_courant, high_face, low_face);
        }

        // Do not let Hydro perform its legacy donor-cell ion advection here.
        // We use the final, limited hydro mass flux below with the same
        // second-order reconstructed ionic face state used everywhere else.
        hydro.do_flux_calculation(
            0, high_state, high_ionization, low_state, low_ionization, dx,
            area, timestep, false);

        HydroVariables &high_variables = high_cell.get_hydro_variables();
        HydroVariables &low_variables0 = low_cell0.get_hydro_variables();
        HydroVariables &low_variables1 = low_cell1.get_hydro_variables();
        for (uint_fast8_t i = 0; i < 5; ++i) {
          high_variables.delta_conserved(i) +=
              high_state.delta_conserved(i);
        }

        // Transform the receiving low-x flux out of the high-x frame.
        double low_delta[5] = {
            low_state.delta_conserved(0), low_state.delta_conserved(1),
            low_state.delta_conserved(2) -
                velocity_jump * low_state.delta_conserved(0),
            low_state.delta_conserved(3),
            low_state.delta_conserved(4) -
                velocity_jump * low_state.delta_conserved(2) +
                0.5 * velocity_jump * velocity_jump *
                    low_state.delta_conserved(0)};
        for (uint_fast8_t i = 0; i < 5; ++i) {
          low_variables0.delta_conserved(i) += weight0 * low_delta[i];
          low_variables1.delta_conserved(i) += weight1 * low_delta[i];
        }

        if (advect_ionization) {
          // high_state started with zero deltas, so its mass change is exactly
          // -mflux after the hydro solver and all hydro flux limiters.
          const double mflux = -high_state.delta_conserved(0);
          const double *upwind_face = mflux > 0. ? high_face : low_face;
          for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
            const double ion_flux = mflux * upwind_face[ion];
            high_cell.get_ionization_variables()
                .increase_delta_ionic_fraction(ion, -ion_flux);
            low_cell0.get_ionization_variables()
                .increase_delta_ionic_fraction(ion, weight0 * ion_flux);
            low_cell1.get_ionization_variables()
                .increase_delta_ionic_fraction(ion, weight1 * ion_flux);
          }
        }
      }
    }
  }

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
