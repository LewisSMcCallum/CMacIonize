/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2026 Lewis McCallum
 *
 * CMacIonize is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *******************************************************************************/

/**
 * @file InitialTurbulence.hpp
 *
 * @brief One-time initialization of an Alvelius turbulent velocity field.
 */
#ifndef INITIALTURBULENCE_HPP
#define INITIALTURBULENCE_HPP

#include "AlveliusTurbulenceForcing.hpp"
#include "DensitySubGridCreator.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

/** @brief Diagnostics for a one-time initial turbulent velocity field. */
struct InitialTurbulenceStatistics {
  CoordinateVector<> raw_mean_velocity;
  double raw_rms_velocity;
  double normalization_factor;
  CoordinateVector<> final_mean_velocity;
  double final_rms_velocity;
  CoordinateVector<> final_component_dispersions;

  InitialTurbulenceStatistics()
      : raw_mean_velocity(0.), raw_rms_velocity(0.),
        normalization_factor(0.), final_mean_velocity(0.),
        final_rms_velocity(0.), final_component_dispersions(0.) {}
};

/** @brief Add a mass-centred, RMS-normalized Alvelius velocity field once. */
class InitialTurbulence {
public:
  static InitialTurbulenceStatistics initialize(
      const bool enabled,
      DensitySubGridCreator< HydroDensitySubGrid > &grid_creator,
      const Hydro &hydro, const double target_rms_velocity,
      const int_fast32_t seed, Log *log = nullptr) {
    InitialTurbulenceStatistics statistics;
    if (!enabled) {
      return statistics;
    }
    if (target_rms_velocity < 0.) {
      cmac_error("Initial turbulence target 3-D RMS velocity must be "
                 "non-negative.");
    }

    // The raw normalization is arbitrary; the requested RMS sets its scale.
    AlveliusTurbulenceForcing forcing(
        grid_creator.get_subgrid_layout(),
        grid_creator.get_subgrid_cell_layout(), grid_creator.get_box(), 1., 3.,
        2.5, 0.2, 1., seed, 1., 0., nullptr);
    forcing.update_turbulence(1.);

    double total_mass = 0.;
    double raw_second_moment = 0.;
    CoordinateVector<> raw_weighted_sum;
    std::vector< CoordinateVector<> > field;
    for (size_t igrid = 0;
         igrid < grid_creator.number_of_original_subgrids(); ++igrid) {
      forcing.get_turbulent_field(igrid, field);
      HydroDensitySubGrid &subgrid = *grid_creator.get_subgrid(igrid);
      size_t icell = 0;
      for (auto cell = subgrid.hydro_begin(); cell != subgrid.hydro_end();
           ++cell, ++icell) {
        const double mass =
            cell.get_hydro_variables().get_conserved_mass();
        total_mass += mass;
        raw_weighted_sum += mass * field[icell];
        raw_second_moment += mass * field[icell].norm2();
      }
    }
    if (!(total_mass > 0.)) {
      cmac_error("Cannot initialize turbulence in a grid with zero total "
                 "mass.");
    }
    statistics.raw_mean_velocity = raw_weighted_sum / total_mass;
    const double raw_variance =
        std::max(0., raw_second_moment / total_mass -
                         statistics.raw_mean_velocity.norm2());
    statistics.raw_rms_velocity = std::sqrt(raw_variance);
    if (!(statistics.raw_rms_velocity > 0.)) {
      cmac_error("Generated initial turbulent field has zero RMS velocity.");
    }
    statistics.normalization_factor =
        target_rms_velocity / statistics.raw_rms_velocity;

    CoordinateVector<> final_weighted_sum;
    CoordinateVector<> final_second_moments;
    for (size_t igrid = 0;
         igrid < grid_creator.number_of_original_subgrids(); ++igrid) {
      forcing.get_turbulent_field(igrid, field);
      HydroDensitySubGrid &subgrid = *grid_creator.get_subgrid(igrid);
      size_t icell = 0;
      for (auto cell = subgrid.hydro_begin(); cell != subgrid.hydro_end();
           ++cell, ++icell) {
        HydroVariables &variables = cell.get_hydro_variables();
        const double mass = variables.get_conserved_mass();
        const CoordinateVector<> turbulent_velocity =
            statistics.normalization_factor *
            (field[icell] - statistics.raw_mean_velocity);
        variables.set_primitives_velocity(
            variables.get_primitives_velocity() + turbulent_velocity);
        hydro.set_conserved_variables(variables, cell.get_volume());

        final_weighted_sum += mass * turbulent_velocity;
        final_second_moments +=
            mass * CoordinateVector<>(
                       turbulent_velocity.x() * turbulent_velocity.x(),
                       turbulent_velocity.y() * turbulent_velocity.y(),
                       turbulent_velocity.z() * turbulent_velocity.z());
      }
    }
    statistics.final_mean_velocity = final_weighted_sum / total_mass;
    statistics.final_component_dispersions = CoordinateVector<>(
        std::sqrt(std::max(
            0., final_second_moments.x() / total_mass -
                    statistics.final_mean_velocity.x() *
                        statistics.final_mean_velocity.x())),
        std::sqrt(std::max(
            0., final_second_moments.y() / total_mass -
                    statistics.final_mean_velocity.y() *
                        statistics.final_mean_velocity.y())),
        std::sqrt(std::max(
            0., final_second_moments.z() / total_mass -
                    statistics.final_mean_velocity.z() *
                        statistics.final_mean_velocity.z())));
    statistics.final_rms_velocity =
        statistics.final_component_dispersions.norm();

    if (log) {
      log->write_status(
          "Initial turbulence raw mass-weighted mean velocity: ",
          statistics.raw_mean_velocity.x(), " ",
          statistics.raw_mean_velocity.y(), " ",
          statistics.raw_mean_velocity.z(), " m s^-1.");
      log->write_status("Initial turbulence raw 3-D RMS velocity: ",
                        statistics.raw_rms_velocity, " m s^-1.");
      log->write_status("Initial turbulence normalization factor: ",
                        statistics.normalization_factor, ".");
      log->write_status(
          "Initial turbulence final mass-weighted mean turbulent velocity: ",
          statistics.final_mean_velocity.x(), " ",
          statistics.final_mean_velocity.y(), " ",
          statistics.final_mean_velocity.z(), " m s^-1.");
      log->write_status("Initial turbulence final 3-D RMS velocity: ",
                        statistics.final_rms_velocity, " m s^-1.");
      log->write_status(
          "Initial turbulence final component dispersions: ",
          statistics.final_component_dispersions.x(), " ",
          statistics.final_component_dispersions.y(), " ",
          statistics.final_component_dispersions.z(), " m s^-1.");
    }
    return statistics;
  }
};

#endif // INITIALTURBULENCE_HPP
