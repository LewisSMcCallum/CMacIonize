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
 * @file testSwiggumFilePhotonSourceDistribution.cpp
 *
 * @brief Restart test for SwiggumFilePhotonSourceDistribution.
 */

#include "Assert.hpp"
#include "DensitySubGridCreator.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "ParameterFile.hpp"
#include "RestartReader.hpp"
#include "RestartWriter.hpp"
#include "SwiggumFilePhotonSourceDistribution.hpp"
#include "UnitConverter.hpp"

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

/** @brief Narrow test access to placement diagnostics and shifted tracks. */
class SwiggumFilePhotonSourceDistributionTestAccess {
private:
  static const SwiggumFilePhotonSourceDistribution::Star &
  star(const SwiggumFilePhotonSourceDistribution &distribution,
       const uint_fast64_t id) {
    for (const auto &entry : distribution._stars) {
      if (entry.id == id) {
        return entry;
      }
    }
    cmac_error("Missing test star.");
    return distribution._stars[0];
  }

  static const SwiggumFilePhotonSourceDistribution::SnapCluster &
  cluster(const SwiggumFilePhotonSourceDistribution &distribution,
          const uint_fast32_t id) {
    for (const auto &entry : distribution._snap_clusters) {
      if (entry.id == id) {
        return entry;
      }
    }
    cmac_error("Missing test cluster.");
    return distribution._snap_clusters[0];
  }

public:
  static CoordinateVector<> offset(
      const SwiggumFilePhotonSourceDistribution &distribution,
      const uint_fast32_t id) {
    return cluster(distribution, id).offset;
  }

  static std::string status(
      const SwiggumFilePhotonSourceDistribution &distribution,
      const uint_fast32_t id) {
    return cluster(distribution, id).status;
  }

  static CoordinateVector<> raw_position(
      const SwiggumFilePhotonSourceDistribution &distribution,
      const uint_fast64_t id, const double time) {
    return distribution.interpolate_position(star(distribution, id), time);
  }

  static CoordinateVector<> shifted_position(
      const SwiggumFilePhotonSourceDistribution &distribution,
      const uint_fast64_t id, const double time) {
    return distribution.get_shifted_position(star(distribution, id), time);
  }

  static CoordinateVector<> feedback_position(
      const SwiggumFilePhotonSourceDistribution &distribution,
      const size_t index) {
    return distribution._feedback_positions[index];
  }

  static CoordinateVector<> source_position(
      const SwiggumFilePhotonSourceDistribution &distribution,
      const uint_fast64_t id) {
    for (size_t i = 0; i < distribution._source_ids.size(); ++i) {
      if (distribution._source_ids[i] == id) {
        return distribution._source_positions[i];
      }
    }
    cmac_error("Missing active test source.");
    return CoordinateVector<>(0.);
  }

  static size_t feedback_count(
      const SwiggumFilePhotonSourceDistribution &distribution) {
    return distribution._feedback_positions.size();
  }

  static size_t inspected_subgrids(
      const SwiggumFilePhotonSourceDistribution &distribution) {
    return distribution._last_snap_subgrids_inspected;
  }
};

typedef DensitySubGridCreator< HydroDensitySubGrid > HydroGridCreator;

static void assert_vector(const CoordinateVector<> value,
                          const CoordinateVector<> expected,
                          const double tolerance = 1.e-8) {
  assert_values_equal_tol(value.x(), expected.x(), tolerance);
  assert_values_equal_tol(value.y(), expected.y(), tolerance);
  assert_values_equal_tol(value.z(), expected.z(), tolerance);
}

static void set_original_gas(HydroGridCreator &grid, const double density,
                             const double temperature,
                             const double neutral_fraction) {
  for (auto subgrid = grid.begin(); subgrid != grid.original_end(); ++subgrid) {
    for (auto cell = (*subgrid).begin(); cell != (*subgrid).end(); ++cell) {
      IonizationVariables &variables = cell.get_ionization_variables();
      variables.set_number_density(density);
      variables.set_temperature(temperature);
      variables.set_ionic_fraction(ION_H_n, neutral_fraction);
    }
  }
}

static void set_cell_gas(HydroGridCreator &grid,
                         const CoordinateVector<> position,
                         const double density, const double temperature,
                         const double neutral_fraction) {
  HydroDensitySubGrid &subgrid = *grid.get_subgrid(position);
  auto cell = subgrid.get_cell(position);
  IonizationVariables &variables = cell.get_ionization_variables();
  variables.set_number_density(density);
  variables.set_temperature(temperature);
  variables.set_ionic_fraction(ION_H_n, neutral_fraction);
}

static size_t count_lines(const std::string &filename) {
  std::ifstream input(filename.c_str());
  size_t lines = 0;
  std::string line;
  while (std::getline(input, line)) {
    ++lines;
  }
  return lines;
}

int main() {
  const double pc = UnitConverter::to_SI< QUANTITY_LENGTH >(1., "pc");
  const double Myr = UnitConverter::to_SI< QUANTITY_TIME >(1., "Myr");
  const double per_cm3 = UnitConverter::to_SI< QUANTITY_NUMBER_DENSITY >(
      1., "cm^-3");

  // Legacy/disabled behaviour remains exactly unchanged, including restart.
  ParameterFile params("test_swiggum_file_photon_source_distribution.yml");
  SwiggumFilePhotonSourceDistribution distribution(params);

  assert_condition(distribution.get_number_of_sources() == 1);
  const CoordinateVector<> position = distribution.get_position(0);
  assert_vector(position,
                CoordinateVector<>(20. * pc, 30. * pc, 40. * pc),
                1.e-12 * pc);
  const double luminosity = distribution.get_total_luminosity();

  {
    RestartWriter restart_writer(
        "test_swiggum_file_photon_source_distribution.dump");
    distribution.write_restart_file(restart_writer);
  }

  RestartReader restart_reader(
      "test_swiggum_file_photon_source_distribution.dump");
  SwiggumFilePhotonSourceDistribution restarted_distribution(restart_reader);

  assert_condition(restarted_distribution.get_number_of_sources() == 1);
  assert_values_equal(restarted_distribution.get_position(0).x(), position.x());
  assert_values_equal(restarted_distribution.get_position(0).y(), position.y());
  assert_values_equal(restarted_distribution.get_position(0).z(), position.z());
  assert_values_equal(restarted_distribution.get_total_luminosity(),
                      luminosity);

  // The filename prefix provides a version marker without invalidating old
  // Swiggum restart records, which ended immediately after these six fields.
  {
    RestartWriter legacy_writer("test_swiggum_legacy.dump");
    legacy_writer.write(
        std::string("test_swiggum_file_photon_source_distribution.csv"));
    CoordinateVector<>(-100. * pc, -200. * pc, -300. * pc)
        .write_restart_file(legacy_writer);
    CoordinateVector<>(500. * pc, 600. * pc, 700. * pc)
        .write_restart_file(legacy_writer);
    legacy_writer.write(0.);
    legacy_writer.write(0.);
    legacy_writer.write(0.1);
    legacy_writer.write(
        UnitConverter::to_SI< QUANTITY_ENERGY >(1.e51, "erg"));
  }
  {
    RestartReader legacy_reader("test_swiggum_legacy.dump");
    SwiggumFilePhotonSourceDistribution legacy(legacy_reader);
    assert_vector(legacy.get_position(0), position, 1.e-12 * pc);
  }

  // Snap seven synthetic clusters through distinct placement outcomes.
  std::remove("test_cluster_snap_events.csv");
  ParameterFile snap_params("test_swiggum_snap.yml");
  SwiggumFilePhotonSourceDistribution snapped(snap_params);
  const Box<> box(CoordinateVector<>(0.),
                  CoordinateVector<>(80. * pc, 80. * pc, 40. * pc));
  HydroGridCreator grid(box, CoordinateVector< int_fast32_t >(16, 16, 8),
                        CoordinateVector< int_fast32_t >(8, 8, 4),
                        CoordinateVector< bool >(false));
  HomogeneousDensityFunction density_function;
  grid.initialize(density_function);
  set_original_gas(grid, 0.5 * per_cm3, 8000., 0.5);

  // Add radiation copies and make them attractive. Placement must still use
  // only the original subgrids.
  const size_t original_subgrids = grid.number_of_original_subgrids();
  std::vector< uint_fast8_t > copy_levels(original_subgrids, 1);
  grid.create_copies(copy_levels);
  for (size_t igrid = original_subgrids;
       igrid < grid.number_of_actual_subgrids(); ++igrid) {
    HydroDensitySubGrid &copy = *grid.get_subgrid(igrid);
    for (auto cell = copy.begin(); cell != copy.end(); ++cell) {
      cell.get_ionization_variables().set_number_density(100. * per_cm3);
      cell.get_ionization_variables().set_temperature(100.);
      cell.get_ionization_variables().set_ionic_fraction(ION_H_n, 1.);
    }
  }

  // Low-density birth site: snap to the nearby cold, neutral dense cell.
  set_cell_gas(grid, CoordinateVector<>(22.5 * pc, 12.5 * pc, 12.5 * pc),
               10. * per_cm3, 100., 1.);
  snapped.update(&grid, 1.01 * Myr);
  assert_condition(SwiggumFilePhotonSourceDistributionTestAccess::status(
                       snapped, 0) == "snapped");
  assert_vector(SwiggumFilePhotonSourceDistributionTestAccess::offset(
                    snapped, 0),
                CoordinateVector<>(10. * pc, 0., 0.), 1.e-7 * pc);
  assert_vector(
      SwiggumFilePhotonSourceDistributionTestAccess::source_position(snapped,
                                                                     1),
      SwiggumFilePhotonSourceDistributionTestAccess::shifted_position(
          snapped, 1, 1.01 * Myr),
      1.e-8 * pc);
  assert_condition(
      SwiggumFilePhotonSourceDistributionTestAccess::inspected_subgrids(
          snapped) < original_subgrids);

  // All members receive one constant translation: separations and velocities
  // are therefore exactly preserved.
  const CoordinateVector<> raw1 =
      SwiggumFilePhotonSourceDistributionTestAccess::raw_position(snapped, 1,
                                                                  2. * Myr);
  const CoordinateVector<> raw2 =
      SwiggumFilePhotonSourceDistributionTestAccess::raw_position(snapped, 2,
                                                                  2. * Myr);
  const CoordinateVector<> shifted1 =
      SwiggumFilePhotonSourceDistributionTestAccess::shifted_position(
          snapped, 1, 2. * Myr);
  const CoordinateVector<> shifted2 =
      SwiggumFilePhotonSourceDistributionTestAccess::shifted_position(
          snapped, 2, 2. * Myr);
  assert_vector(shifted1 - raw1, shifted2 - raw2, 1.e-8 * pc);
  assert_vector(shifted2 - shifted1, raw2 - raw1, 1.e-8 * pc);
  const CoordinateVector<> raw_step =
      SwiggumFilePhotonSourceDistributionTestAccess::raw_position(
          snapped, 1, 2.5 * Myr) -
      SwiggumFilePhotonSourceDistributionTestAccess::raw_position(
          snapped, 1, 2. * Myr);
  const CoordinateVector<> shifted_step =
      SwiggumFilePhotonSourceDistributionTestAccess::shifted_position(
          snapped, 1, 2.5 * Myr) -
      SwiggumFilePhotonSourceDistributionTestAccess::shifted_position(
          snapped, 1, 2. * Myr);
  assert_vector(shifted_step, raw_step, 1.e-8 * pc);

  // Restart preserves the offset and does not process/log the cluster again.
  {
    RestartWriter writer("test_swiggum_snap.dump");
    snapped.write_restart_file(writer);
  }
  const size_t log_lines_before_restart =
      count_lines("test_cluster_snap_events.csv");
  {
    RestartReader reader("test_swiggum_snap.dump");
    SwiggumFilePhotonSourceDistribution restarted(reader);
    assert_vector(SwiggumFilePhotonSourceDistributionTestAccess::offset(
                      restarted, 0),
                  CoordinateVector<>(10. * pc, 0., 0.), 1.e-7 * pc);
    restarted.update(&grid, 0.1 * Myr);
    assert_condition(count_lines("test_cluster_snap_events.csv") ==
                     log_lines_before_restart);
  }

  // Suitable nominal gas: no translation.
  set_original_gas(grid, 0.5 * per_cm3, 8000., 0.5);
  set_cell_gas(grid, CoordinateVector<>(32.5 * pc, 12.5 * pc, 12.5 * pc),
               5. * per_cm3, 100., 1.);
  snapped.update(&grid, 1. * Myr);
  assert_condition(SwiggumFilePhotonSourceDistributionTestAccess::status(
                       snapped, 1) == "already_in_gas");
  assert_vector(SwiggumFilePhotonSourceDistributionTestAccess::offset(
                    snapped, 1),
                CoordinateVector<>(0.));

  // Outside-box births and searches with no eligible original cell are safe.
  set_original_gas(grid, 0.5 * per_cm3, 8000., 0.5);
  snapped.update(&grid, 1. * Myr);
  assert_condition(SwiggumFilePhotonSourceDistributionTestAccess::status(
                       snapped, 2) == "outside_box_at_birth");
  assert_vector(SwiggumFilePhotonSourceDistributionTestAccess::offset(
                    snapped, 2),
                CoordinateVector<>(0.));
  snapped.update(&grid, 1. * Myr);
  assert_condition(SwiggumFilePhotonSourceDistributionTestAccess::status(
                       snapped, 3) == "no_eligible_target");
  assert_vector(SwiggumFilePhotonSourceDistributionTestAccess::offset(
                    snapped, 3),
                CoordinateVector<>(0.));
  assert_condition(
      SwiggumFilePhotonSourceDistributionTestAccess::inspected_subgrids(
          snapped) < original_subgrids);

  // A valid target beyond the allowed displacement is rejected.
  set_original_gas(grid, 0.5 * per_cm3, 8000., 0.5);
  set_cell_gas(grid, CoordinateVector<>(67.5 * pc, 12.5 * pc, 12.5 * pc),
               10. * per_cm3, 100., 1.);
  snapped.update(&grid, 1. * Myr);
  assert_condition(SwiggumFilePhotonSourceDistributionTestAccess::status(
                       snapped, 4) == "maximum_displacement_exceeded");
  assert_vector(SwiggumFilePhotonSourceDistributionTestAccess::offset(
                    snapped, 4),
                CoordinateVector<>(0.));

  // Births at different times inspect the live grid state at their own birth.
  set_original_gas(grid, 0.5 * per_cm3, 8000., 0.5);
  set_cell_gas(grid, CoordinateVector<>(22.5 * pc, 52.5 * pc, 12.5 * pc),
               10. * per_cm3, 100., 1.);
  snapped.update(&grid, 1. * Myr);
  set_original_gas(grid, 0.5 * per_cm3, 8000., 0.5);
  set_cell_gas(grid, CoordinateVector<>(2.5 * pc, 52.5 * pc, 12.5 * pc),
               10. * per_cm3, 100., 1.);
  snapped.update(&grid, 1. * Myr);
  assert_vector(SwiggumFilePhotonSourceDistributionTestAccess::offset(
                    snapped, 5),
                CoordinateVector<>(10. * pc, 0., 0.), 1.e-7 * pc);
  assert_vector(SwiggumFilePhotonSourceDistributionTestAccess::offset(
                    snapped, 6),
                CoordinateVector<>(-10. * pc, 0., 0.), 1.e-7 * pc);

  // The later supernova follows the same permanently shifted stellar track.
  snapped.update(&grid, 1. * Myr);
  assert_condition(
      SwiggumFilePhotonSourceDistributionTestAccess::feedback_count(snapped) ==
      1);
  assert_vector(
      SwiggumFilePhotonSourceDistributionTestAccess::feedback_position(snapped,
                                                                       0),
      SwiggumFilePhotonSourceDistributionTestAccess::shifted_position(
          snapped, 1, 8. * Myr),
      1.e-8 * pc);

  return 0;
}
