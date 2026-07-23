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
#include "ParameterFile.hpp"
#include "RestartReader.hpp"
#include "RestartWriter.hpp"
#include "SwiggumFilePhotonSourceDistribution.hpp"

int main() {
  ParameterFile params("test_swiggum_file_photon_source_distribution.yml");
  SwiggumFilePhotonSourceDistribution distribution(params);

  assert_condition(distribution.get_number_of_sources() == 1);
  const CoordinateVector<> position = distribution.get_position(0);
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

  return 0;
}
