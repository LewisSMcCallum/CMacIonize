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
 * @brief One-time initialization of a turbulent velocity field.
 */
#ifndef INITIALTURBULENCE_HPP
#define INITIALTURBULENCE_HPP

#include "AlveliusTurbulenceForcing.hpp"
#include "DensitySubGridCreator.hpp"
#include "ParameterFile.hpp"
#include "SampledInitialTurbulence.hpp"

#include <algorithm>
#include <cmath>
#include <string>
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

/** @brief Add a mass-centred, RMS-normalized turbulent velocity field once. */
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

    const Box<> box = grid_creator.get_box();
    const double mode_length =
        std::max(box.get_sides().x(),
                 std::max(box.get_sides().y(), box.get_sides().z()));

    double minimum_wave_number = 1.;
    double maximum_wave_number = 3.;
    double peak_wave_number = 2.5;
    double concentration_factor = 0.2;
    uint_fast32_t number_of_modes = 0;
    std::string components = "xyz";
    bool wavelengths_set = false;

    ParameterFile *params = ParameterFile::get_active_parameter_file();
    if (params != nullptr) {
      components = params->get_value< std::string >(
          "InitialTurbulence:components", components);
      number_of_modes = params->get_value< uint_fast32_t >(
          "InitialTurbulence:number of modes", number_of_modes);

      const bool has_largest =
          params->has_value("InitialTurbulence:largest wavelength");
      const bool has_smallest =
          params->has_value("InitialTurbulence:smallest wavelength");
      const bool has_peak =
          params->has_value("InitialTurbulence:peak wavelength");
      wavelengths_set = has_largest || has_smallest || has_peak;
      if (wavelengths_set) {
        if (!(has_largest && has_smallest && has_peak)) {
          cmac_error("Initial turbulence physical scales require largest, "
                     "smallest and peak wavelength to all be specified.");
        }
        const double largest_wavelength =
            params->get_physical_value< QUANTITY_LENGTH >(
                "InitialTurbulence:largest wavelength");
        const double smallest_wavelength =
            params->get_physical_value< QUANTITY_LENGTH >(
                "InitialTurbulence:smallest wavelength");
        const double peak_wavelength =
            params->get_physical_value< QUANTITY_LENGTH >(
                "InitialTurbulence:peak wavelength");
        if (!(smallest_wavelength > 0.) ||
            peak_wavelength < smallest_wavelength ||
            largest_wavelength < peak_wavelength) {
          cmac_error("Invalid initial turbulence wavelength interval: require "
                     "0 < smallest <= peak <= largest.");
        }
        minimum_wave_number = mode_length / largest_wavelength;
        maximum_wave_number = mode_length / smallest_wavelength;
        peak_wave_number = mode_length / peak_wavelength;
        if (params->has_value("InitialTurbulence:concentration factor")) {
          concentration_factor = params->get_value< double >(
              "InitialTurbulence:concentration factor", concentration_factor);
        } else {
          // Keep the spectral width fixed in physical rather than box units.
          concentration_factor =
              0.5 * (maximum_wave_number - minimum_wave_number);
        }
      } else {
        minimum_wave_number = params->get_value< double >(
            "InitialTurbulence:minimum wave number", minimum_wave_number);
        maximum_wave_number = params->get_value< double >(
            "InitialTurbulence:maximum wave number", maximum_wave_number);
        peak_wave_number = params->get_value< double >(
            "InitialTurbulence:peak wave number", peak_wave_number);
        concentration_factor = params->get_value< double >(
            "InitialTurbulence:concentration factor", concentration_factor);
      }
    }

    const bool vertical_only = components == "z" || components == "vertical";
    if (!vertical_only && components != "xyz") {
      cmac_error("Unknown InitialTurbulence:components value '%s'. Expected "
                 "'xyz' or 'z'.",
                 components.c_str());
    }
    if (!(minimum_wave_number > 0.) ||
        maximum_wave_number < minimum_wave_number ||
        peak_wave_number < minimum_wave_number ||
        peak_wave_number > maximum_wave_number) {
      cmac_error("Invalid initial turbulence spectrum: kmin=%g, kpeak=%g, "
                 "kmax=%g.",
                 minimum_wave_number, peak_wave_number, maximum_wave_number);
    }
    if (!(concentration_factor > 0.)) {
      cmac_error("Initial turbulence concentration factor must be positive.");
    }

    // Vertical-only turbulence needs the sampled-mode implementation, because
    // it deliberately uses k_z=0 modes with amplitudes polarized along z.
    // This gives v=(0,0,v_z(x,y)): no imposed x-y velocity and zero initial
    // divergence, while retaining a finite vertical velocity dispersion.
    if (vertical_only && number_of_modes == 0) {
      number_of_modes = 32;
    }

    if (log) {
      log->write_status("Initial turbulence components: ",
                        vertical_only ? "z only (k_z=0)" : "xyz", ".");
      log->write_status("Initial turbulence spectrum: |k| = ",
                        minimum_wave_number, "--", maximum_wave_number,
                        ", peak = ", peak_wave_number,
                        ", concentration = ", concentration_factor,
                        ", sampled modes = ", number_of_modes,
                        number_of_modes > 0 ? "." : " (legacy full shell).",
                        wavelengths_set ? " Physical wavelengths requested."
                                        : "");
    }

    AlveliusTurbulenceForcing *forcing = nullptr;
    SampledInitialTurbulence *sampled_forcing = nullptr;
    if (number_of_modes > 0) {
      sampled_forcing = new SampledInitialTurbulence(
          grid_creator.get_subgrid_layout(),
          grid_creator.get_subgrid_cell_layout(), box, minimum_wave_number,
          maximum_wave_number, peak_wave_number, concentration_factor,
          number_of_modes, seed, log, vertical_only);
    } else {
      // A zero mode count preserves the previous full-shell implementation.
      forcing = new AlveliusTurbulenceForcing(
          grid_creator.get_subgrid_layout(),
          grid_creator.get_subgrid_cell_layout(), box, minimum_wave_number,
          maximum_wave_number, peak_wave_number, concentration_factor, 1., seed,
          1., 0., nullptr);
      forcing->update_turbulence(1.);
    }

    double total_mass = 0.;
    double raw_second_moment = 0.;
    CoordinateVector<> raw_weighted_sum;
    std::vector< CoordinateVector<> > field;
    for (size_t igrid = 0;
         igrid < grid_creator.number_of_original_subgrids(); ++igrid) {
      if (sampled_forcing != nullptr) {
        sampled_forcing->get_turbulent_field(igrid, field);
      } else {
        forcing->get_turbulent_field(igrid, field);
      }
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
      if (sampled_forcing != nullptr) {
        sampled_forcing->get_turbulent_field(igrid, field);
      } else {
        forcing->get_turbulent_field(igrid, field);
      }
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

    delete forcing;
    delete sampled_forcing;

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
