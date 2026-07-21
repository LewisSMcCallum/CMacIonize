/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testHydroDensitySubGrid.cpp
 *
 * @brief Unit test for the HydroDensitySubGrid class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Assert.hpp"
#include "HydroDensitySubGrid.hpp"

#include <fstream>
#include <vector>

/**
 * @brief Check conservative first-order advection around a periodic 1-D grid.
 *
 * Uniform hydro states give an exactly known upwind update:
 * x_i(new) = (1-C) x_i(old) + C x_(i-1)(old), where C is the Courant number.
 */
void test_ionization_advection(const Hydro &hydro) {
  const int_fast32_t number_of_cells = 16;
  const int_fast32_t cells_per_subgrid = number_of_cells / 2;
  const double left_box[6] = {0., 0., 0., 0.5, 1., 1.};
  const double right_box[6] = {0.5, 0., 0., 0.5, 1., 1.};
  HydroDensitySubGrid left_grid(
      left_box, CoordinateVector< int_fast32_t >(cells_per_subgrid, 1, 1));
  HydroDensitySubGrid right_grid(
      right_box, CoordinateVector< int_fast32_t >(cells_per_subgrid, 1, 1));

  std::vector< double > initial_fraction(number_of_cells);
  int_fast32_t index = 0;
  HydroDensitySubGrid *grids[2] = {&left_grid, &right_grid};
  for (int_fast32_t igrid = 0; igrid < 2; ++igrid) {
    for (auto cell = grids[igrid]->hydro_begin();
         cell != grids[igrid]->hydro_end(); ++cell) {
      cell.get_hydro_variables().set_primitives_density(1.);
      cell.get_hydro_variables().set_primitives_velocity(
          CoordinateVector<>(1., 0., 0.));
      cell.get_hydro_variables().set_primitives_pressure(1.);
      initial_fraction[index] = (index < number_of_cells / 2) ? 0.2 : 0.8;
      cell.get_ionization_variables().set_ionic_fraction(
          ION_H_n, initial_fraction[index]);
      ++index;
    }
    grids[igrid]->initialize_hydrodynamic_variables(hydro, false);
  }

  const double courant_number = 0.25;
  const double timestep = courant_number / number_of_cells;
  left_grid.inner_flux_sweep(hydro, timestep, true);
  right_grid.inner_flux_sweep(hydro, timestep, true);
  left_grid.outer_flux_sweep(TRAVELDIRECTION_FACE_X_P, hydro, right_grid,
                             timestep, true);
  // Connect the second subgrid back to the first across the periodic boundary.
  right_grid.outer_flux_sweep(TRAVELDIRECTION_FACE_X_P, hydro, left_grid,
                              timestep, true);
  left_grid.update_conserved_variables(timestep, true);
  right_grid.update_conserved_variables(timestep, true);
  left_grid.update_primitive_variables(hydro);
  right_grid.update_primitive_variables(hydro);

  double initial_neutral_mass = 0.;
  double final_neutral_mass = 0.;
  index = 0;
  for (int_fast32_t igrid = 0; igrid < 2; ++igrid) {
    for (auto cell = grids[igrid]->hydro_begin();
         cell != grids[igrid]->hydro_end(); ++cell) {
      const int_fast32_t upwind =
          (index + number_of_cells - 1) % number_of_cells;
      const double expected =
          (1. - courant_number) * initial_fraction[index] +
          courant_number * initial_fraction[upwind];
      const double actual =
          cell.get_ionization_variables().get_ionic_fraction(ION_H_n);
      assert_values_equal_tol(actual, expected, 1.e-12);
      assert_condition(actual >= 0. && actual <= 1.);
      initial_neutral_mass += initial_fraction[index] / number_of_cells;
      final_neutral_mass +=
          cell.get_hydro_variables().get_conserved_mass() * actual;
      ++index;
    }
  }
  assert_values_equal_tol(final_neutral_mass, initial_neutral_mass, 1.e-12);
}

/**
 * @brief Unit test for the HydroDensitySubGrid class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const double box1[6] = {-0.5, -0.25, -0.25, 0.5, 0.5, 0.5};
  const double box2[6] = {0., -0.25, -0.25, 0.5, 0.5, 0.5};
  const CoordinateVector< int_fast32_t > ncell(50, 3, 3);
  HydroDensitySubGrid test_grid1(box1, ncell);
  HydroDensitySubGrid test_grid2(box2, ncell);

  for (auto cellit = test_grid1.hydro_begin(); cellit != test_grid1.hydro_end();
       ++cellit) {
    cellit.get_hydro_variables().set_primitives_density(1.);
    cellit.get_hydro_variables().set_primitives_pressure(1.);
    cellit.get_ionization_variables().set_ionic_fraction(ION_H_n, 1.);
  }

  for (auto cellit = test_grid2.hydro_begin(); cellit != test_grid2.hydro_end();
       ++cellit) {
    cellit.get_hydro_variables().set_primitives_density(0.125);
    cellit.get_hydro_variables().set_primitives_pressure(0.1);
    cellit.get_ionization_variables().set_ionic_fraction(ION_H_n, 0.1);
  }

  const double dt = 0.001;
  const Abundances abundances;
  const Hydro hydro(5. / 3., 100., 1.e4, 1.e99, false, abundances);
  const InflowHydroBoundary inflow_boundary;
  const ReflectiveHydroBoundary reflective_boundary;

  test_grid1.initialize_hydrodynamic_variables(hydro, false);
  test_grid2.initialize_hydrodynamic_variables(hydro, false);

  test_ionization_advection(hydro);

  // A constant gravity kick should change kinetic, but not internal, energy.
  const double gravity_box[6] = {0., 0., 0., 1., 1., 1.};
  HydroDensitySubGrid gravity_grid(
      gravity_box, CoordinateVector< int_fast32_t >(1, 1, 1));
  auto gravity_cell = gravity_grid.hydro_begin();
  gravity_cell.get_hydro_variables().set_primitives_density(2.);
  gravity_cell.get_hydro_variables().set_primitives_velocity(
      CoordinateVector<>(1., 2., 3.));
  gravity_cell.get_hydro_variables().set_primitives_pressure(4.);
  gravity_cell.get_hydro_variables().set_gravitational_acceleration(
      CoordinateVector<>(4., 5., 6.));
  gravity_grid.initialize_hydrodynamic_variables(hydro, false);
  const double old_internal_energy =
      gravity_cell.get_hydro_variables().get_conserved_total_energy() -
      0.5 * gravity_cell.get_hydro_variables().get_conserved_momentum().norm2() /
          gravity_cell.get_hydro_variables().get_conserved_mass();
  gravity_grid.update_conserved_variables(0.01);
  const double new_internal_energy =
      gravity_cell.get_hydro_variables().get_conserved_total_energy() -
      0.5 * gravity_cell.get_hydro_variables().get_conserved_momentum().norm2() /
          gravity_cell.get_hydro_variables().get_conserved_mass();
  assert_values_equal_tol(old_internal_energy, new_internal_energy, 1.e-12);

  for (uint_fast32_t istep = 0; istep < 300; ++istep) {

    // gradient calculations
    test_grid1.inner_gradient_sweep(hydro);
    test_grid2.inner_gradient_sweep(hydro);
    test_grid1.outer_gradient_sweep(TRAVELDIRECTION_FACE_X_P, hydro,
                                    test_grid2);
    test_grid1.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_X_N, hydro,
                                          inflow_boundary);
    test_grid2.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_X_P, hydro,
                                          reflective_boundary);

    test_grid1.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_Y_N, hydro,
                                          inflow_boundary);
    test_grid1.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_Y_P, hydro,
                                          inflow_boundary);
    test_grid1.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_Z_N, hydro,
                                          inflow_boundary);
    test_grid1.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_Z_P, hydro,
                                          inflow_boundary);
    test_grid2.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_Y_N, hydro,
                                          inflow_boundary);
    test_grid2.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_Y_P, hydro,
                                          inflow_boundary);
    test_grid2.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_Z_N, hydro,
                                          inflow_boundary);
    test_grid2.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_Z_P, hydro,
                                          inflow_boundary);

    // slope limiting
    test_grid1.apply_slope_limiter(hydro);
    test_grid2.apply_slope_limiter(hydro);

    // second order time prediction
    test_grid1.predict_primitive_variables(hydro, 0.5 * dt);
    test_grid2.predict_primitive_variables(hydro, 0.5 * dt);

    // flux exchanges
    test_grid1.inner_flux_sweep(hydro, dt, true);
    test_grid2.inner_flux_sweep(hydro, dt, true);
    test_grid1.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_X_N, hydro,
                                      inflow_boundary, dt, true);
    test_grid1.outer_flux_sweep(TRAVELDIRECTION_FACE_X_P, hydro, test_grid2,
                                dt, true);
    test_grid2.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_X_P, hydro,
                                      reflective_boundary, dt, true);

    test_grid1.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_Y_N, hydro,
                                      inflow_boundary, dt, true);
    test_grid1.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_Y_P, hydro,
                                      inflow_boundary, dt, true);
    test_grid1.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_Z_N, hydro,
                                      inflow_boundary, dt, true);
    test_grid1.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_Z_P, hydro,
                                      inflow_boundary, dt, true);
    test_grid2.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_Y_N, hydro,
                                      inflow_boundary, dt, true);
    test_grid2.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_Y_P, hydro,
                                      inflow_boundary, dt, true);
    test_grid2.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_Z_N, hydro,
                                      inflow_boundary, dt, true);
    test_grid2.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_Z_P, hydro,
                                      inflow_boundary, dt, true);

    // conserved variable update
    test_grid1.update_conserved_variables(dt, true);
    test_grid2.update_conserved_variables(dt, true);

    // primitive variable update
    test_grid1.update_primitive_variables(hydro);
    test_grid2.update_primitive_variables(hydro);
  }

  // verify that ionization state was indeed advected!
  bool gained_neutral_h = false;
  for (auto cellit = test_grid2.hydro_begin(); cellit != test_grid2.hydro_end(); ++cellit) {
    if (cellit.get_ionization_variables().get_ionic_fraction(ION_H_n) > 0.11) {
      gained_neutral_h = true;
      break;
    }
  }
  assert_condition(gained_neutral_h);

  /// write a restart file
  {
    RestartWriter writer("test_hydrodensitysubgrid.restart");
    test_grid1.write_restart_file(writer);
  }

  /// read the restart file and check that both grids are the same
  {
    RestartReader reader("test_hydrodensitysubgrid.restart");
    HydroDensitySubGrid grid2(reader);
    assert_condition(test_grid1.get_number_of_cells() ==
                     grid2.get_number_of_cells());
    auto it = test_grid1.hydro_begin();
    auto it2 = grid2.hydro_begin();
    while (it != test_grid1.hydro_end() && it2 != grid2.hydro_end()) {
      assert_condition(
          it.get_ionization_variables().get_ionic_fraction(ION_H_n) ==
          it2.get_ionization_variables().get_ionic_fraction(ION_H_n));
      assert_condition(it.get_hydro_variables().get_primitives_density() ==
                       it2.get_hydro_variables().get_primitives_density());
      ++it;
      ++it2;
    }
  }

  std::ofstream ofile("testHydroDensitySubGrid_result.txt");
  ofile << "# x\trho\tvx\tP\n";
  for (auto cellit = test_grid1.hydro_begin(); cellit != test_grid1.hydro_end();
       ++cellit) {
    const CoordinateVector<> p = cellit.get_cell_midpoint();
    const HydroVariables &hydrovars = cellit.get_hydro_variables();
    ofile << p.x() << "\t" << hydrovars.get_primitives_density() << "\t"
          << hydrovars.get_primitives_velocity().x() << "\t"
          << hydrovars.get_primitives_pressure() << "\n";
  }
  for (auto cellit = test_grid2.hydro_begin(); cellit != test_grid2.hydro_end();
       ++cellit) {
    const CoordinateVector<> p = cellit.get_cell_midpoint();
    const HydroVariables &hydrovars = cellit.get_hydro_variables();
    ofile << p.x() << "\t" << hydrovars.get_primitives_density() << "\t"
          << hydrovars.get_primitives_velocity().x() << "\t"
          << hydrovars.get_primitives_pressure() << "\n";
  }

  return 0;
}
