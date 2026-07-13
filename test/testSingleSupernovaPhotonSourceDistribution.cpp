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

  const double box[6] = {-2. * pc_in_m, -2. * pc_in_m, -2. * pc_in_m,
                         4. * pc_in_m,  4. * pc_in_m,  4. * pc_in_m};
  const CoordinateVector< int_fast32_t > cell_layout(4, 4, 4);
  const Abundances abundances;
  Hydro hydro(5. / 3., 100., 1.e4, 1.e99, false, abundances);
  SupernovaHandler handler(1.e44);

  HydroDensitySubGrid thermal_grid(box, cell_layout);
  for (auto cell = thermal_grid.hydro_begin();
       cell != thermal_grid.hydro_end(); ++cell) {
    cell.get_hydro_variables().set_primitives_density(1.67262192e-21);
    cell.get_hydro_variables().set_primitives_pressure(1.e-13);
  }
  thermal_grid.initialize_hydrodynamic_variables(hydro, false);
  handler.inject_sne(thermal_grid, hydro, position, 10. * pc_in_m,
                     40. * pc_in_m, 1., 64);
  handler.inject_sne(thermal_grid, hydro, position, 10. * pc_in_m,
                     40. * pc_in_m, 1., 64);
  for (auto cell = thermal_grid.hydro_begin();
       cell != thermal_grid.hydro_end(); ++cell) {
    assert_values_equal_rel(cell.get_hydro_variables().get_energy_term(),
                            2.e44 / 64., 1.e-12);
  }

  HydroDensitySubGrid momentum_grid(box, cell_layout);
  for (auto cell = momentum_grid.hydro_begin();
       cell != momentum_grid.hydro_end(); ++cell) {
    cell.get_hydro_variables().set_primitives_density(1.67262192e-21);
    cell.get_hydro_variables().set_primitives_pressure(1.e-13);
  }
  momentum_grid.initialize_hydrodynamic_variables(hydro, false);
  // r_inj is not below r_st/3, so this remnant is unresolved.
  handler.inject_sne(momentum_grid, hydro, position, 2.5 * pc_in_m,
                     6. * pc_in_m, 1., 56);
  bool injected_momentum = false;
  for (auto cell = momentum_grid.hydro_begin();
       cell != momentum_grid.hydro_end(); ++cell) {
    injected_momentum =
        injected_momentum ||
        cell.get_hydro_variables().get_primitives_velocity().norm() > 0.;
    assert_condition(cell.get_hydro_variables().get_energy_term() == 0.);
  }
  assert_condition(injected_momentum);

  // restart test
  {
    RestartReader restart_reader(
        "singlesupernovaphotonsourcedistribution.dump");
    SingleSupernovaPhotonSourceDistribution restart_distribution(
        restart_reader);
  }

  return 0;
}
