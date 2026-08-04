/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * CMacIonize is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CMacIonize is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file testSingleSupernovaPhotonSourceDistribution.cpp
 *
 * @brief Unit test for the SingleSupernovaPhotonSourceDistribution class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "HydroDensitySubGrid.hpp"
#include "SingleSupernovaPhotonSourceDistribution.hpp"
#include "SupernovaHandler.hpp"

/**
 * @brief Unit test for the SingleSupernovaPhotonSourceDistribution class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const double Myr_in_s = 1.e6 * 365.25 * 24. * 3600.;
  const double pc_in_m = 3.086e16;

  CoordinateVector<> position;
  SingleSupernovaPhotonSourceDistribution distribution(position, 10. * Myr_in_s,
                                                       1.e49, 1.e44);

  // restart test
  {
    RestartWriter restart_writer(
        "singlesupernovaphotonsourcedistribution.dump");
    distribution.write_restart_file(restart_writer);
  }

  const double box[6] = {-5. * pc_in_m, -5. * pc_in_m, -5. * pc_in_m,
                         10. * pc_in_m, 10. * pc_in_m, 10. * pc_in_m};
  const CoordinateVector< int_fast32_t > cell_layout(10, 10, 10);
  const Abundances abundances;
  Hydro hydro(5. / 3., 100., 1.e4, 1.e99, false, abundances);
  SupernovaHandler handler(1.e44);
  assert_condition(handler.uses_tigress_like_injection());

  HydroDensitySubGrid thermal_grid(box, cell_layout);
  for (auto cell = thermal_grid.hydro_begin();
       cell != thermal_grid.hydro_end(); ++cell) {
    cell.get_hydro_variables().set_primitives_density(1.67262192e-21);
    cell.get_hydro_variables().set_primitives_pressure(1.e-13);
  }
  thermal_grid.initialize_hydrodynamic_variables(hydro, false);

  const double injection_radius = 3. * pc_in_m;
  const double injection_volume =
      4. * M_PI * std::pow(injection_radius, 3.) / 3.;
  double thermal_feedback_mass = 0.;
  int thermal_cell_count = 0;
  for (auto cell = thermal_grid.hydro_begin();
       cell != thermal_grid.hydro_end(); ++cell) {
    if ((cell.get_cell_midpoint() - position).norm() <= injection_radius) {
      thermal_feedback_mass +=
          cell.get_hydro_variables().get_conserved_mass();
      ++thermal_cell_count;
    }
  }
  const double thermal_nbar =
      1.e-6 * thermal_feedback_mass /
      (injection_volume * 1.67262192369e-27);

  handler.inject_sne(thermal_grid, hydro, position, injection_radius,
                     20. * pc_in_m, thermal_nbar, thermal_cell_count);
  handler.inject_sne(thermal_grid, hydro, position, injection_radius,
                     20. * pc_in_m, thermal_nbar, thermal_cell_count);

  double injected_thermal_energy = 0.;
  double injected_kinetic_energy = 0.;
  for (auto cell = thermal_grid.hydro_begin();
       cell != thermal_grid.hydro_end(); ++cell) {
    injected_thermal_energy +=
        cell.get_hydro_variables().get_energy_term();
    injected_kinetic_energy +=
        0.5 * CoordinateVector<>::dot_product(
                  cell.get_hydro_variables().get_primitives_velocity(),
                  cell.get_hydro_variables().get_conserved_momentum());
  }
  assert_values_equal_rel(injected_thermal_energy, 2. * 0.72e44, 1.e-12);
  assert_values_equal_rel(injected_kinetic_energy, 2. * 0.28e44, 1.e-12);

  HydroDensitySubGrid momentum_grid(box, cell_layout);
  for (auto cell = momentum_grid.hydro_begin();
       cell != momentum_grid.hydro_end(); ++cell) {
    cell.get_hydro_variables().set_primitives_density(1.67262192e-19);
    cell.get_hydro_variables().set_primitives_pressure(1.e-13);
  }
  momentum_grid.initialize_hydrodynamic_variables(hydro, false);

  double momentum_feedback_mass = 0.;
  int momentum_cell_count = 0;
  for (auto cell = momentum_grid.hydro_begin();
       cell != momentum_grid.hydro_end(); ++cell) {
    if ((cell.get_cell_midpoint() - position).norm() <= injection_radius) {
      momentum_feedback_mass +=
          cell.get_hydro_variables().get_conserved_mass();
      ++momentum_cell_count;
    }
  }
  const double momentum_nbar =
      1.e-6 * momentum_feedback_mass /
      (injection_volume * 1.67262192369e-27);
  handler.inject_sne(momentum_grid, hydro, position, injection_radius,
                     5. * pc_in_m, momentum_nbar, momentum_cell_count);

  double injected_radial_momentum = 0.;
  for (auto cell = momentum_grid.hydro_begin();
       cell != momentum_grid.hydro_end(); ++cell) {
    const CoordinateVector<> displacement =
        cell.get_cell_midpoint() - position;
    if (displacement.norm() <= injection_radius) {
      injected_radial_momentum += CoordinateVector<>::dot_product(
          cell.get_hydro_variables().get_conserved_momentum(),
          displacement / displacement.norm());
    }
    assert_condition(cell.get_hydro_variables().get_energy_term() == 0.);
  }
  const double expected_terminal_momentum =
      2.8e5 * 1.98847e30 * 1.e3 * std::pow(momentum_nbar, -0.17);
  assert_values_equal_rel(injected_radial_momentum,
                          expected_terminal_momentum, 1.e-12);

  // The optional SILCC-like branch retains the previous pure-thermal versus
  // terminal-momentum behaviour and old r_sf < 4*dx criterion.
  SupernovaHandler silcc_handler(1.e44, false);
  assert_condition(!silcc_handler.uses_tigress_like_injection());

  HydroDensitySubGrid silcc_thermal_grid(box, cell_layout);
  for (auto cell = silcc_thermal_grid.hydro_begin();
       cell != silcc_thermal_grid.hydro_end(); ++cell) {
    cell.get_hydro_variables().set_primitives_density(1.67262192e-21);
    cell.get_hydro_variables().set_primitives_pressure(1.e-13);
  }
  silcc_thermal_grid.initialize_hydrodynamic_variables(hydro, false);
  silcc_handler.inject_sne(silcc_thermal_grid, hydro, position,
                           injection_radius, 5. * pc_in_m, thermal_nbar,
                           thermal_cell_count);
  double silcc_thermal_energy = 0.;
  for (auto cell = silcc_thermal_grid.hydro_begin();
       cell != silcc_thermal_grid.hydro_end(); ++cell) {
    silcc_thermal_energy += cell.get_hydro_variables().get_energy_term();
  }
  assert_values_equal_rel(silcc_thermal_energy, 1.e44, 1.e-12);

  HydroDensitySubGrid silcc_momentum_grid(box, cell_layout);
  for (auto cell = silcc_momentum_grid.hydro_begin();
       cell != silcc_momentum_grid.hydro_end(); ++cell) {
    cell.get_hydro_variables().set_primitives_density(1.67262192e-21);
    cell.get_hydro_variables().set_primitives_pressure(1.e-13);
  }
  silcc_momentum_grid.initialize_hydrodynamic_variables(hydro, false);
  silcc_handler.inject_sne(silcc_momentum_grid, hydro, position,
                           injection_radius, 2. * pc_in_m, thermal_nbar,
                           thermal_cell_count);
  double silcc_radial_momentum = 0.;
  for (auto cell = silcc_momentum_grid.hydro_begin();
       cell != silcc_momentum_grid.hydro_end(); ++cell) {
    const CoordinateVector<> displacement =
        cell.get_cell_midpoint() - position;
    if (displacement.norm() <= injection_radius) {
      silcc_radial_momentum += CoordinateVector<>::dot_product(
          cell.get_hydro_variables().get_conserved_momentum(),
          displacement / displacement.norm());
    }
    assert_condition(cell.get_hydro_variables().get_energy_term() == 0.);
  }
  const double expected_silcc_momentum =
      2.6e5 * 1.98847e30 * 1.e3 * std::pow(thermal_nbar, -2. / 17.);
  assert_values_equal_rel(silcc_radial_momentum,
                          expected_silcc_momentum, 1.e-12);

  // restart test
  {
    RestartReader restart_reader(
        "singlesupernovaphotonsourcedistribution.dump");
    SingleSupernovaPhotonSourceDistribution restart_distribution(
        restart_reader);
  }

  return 0;
}
