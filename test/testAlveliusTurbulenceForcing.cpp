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
 * @file testAlveliusTurbulenceForcing.cpp
 *
 * @brief Unit test for the AlveliusTurbulenceForcing class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "AlveliusTurbulenceForcing.hpp"
#include "Assert.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "InitialTurbulence.hpp"

typedef DensitySubGridCreator< HydroDensitySubGrid > HydroGridCreator;

/** @brief Create initialized hydro variables with an optional linear shear. */
static void initialize_rectangular_grid(HydroGridCreator &grid_creator,
                                        const Hydro &hydro,
                                        const double shear_gradient) {
  HomogeneousDensityFunction density_function;
  grid_creator.initialize(density_function);
  for (auto grid = grid_creator.begin(); grid != grid_creator.original_end();
       ++grid) {
    for (auto cell = (*grid).hydro_begin(); cell != (*grid).hydro_end();
         ++cell) {
      HydroVariables &variables = cell.get_hydro_variables();
      const double x = cell.get_cell_midpoint().x();
      variables.set_primitives_density(2. + 0.25 * x);
      variables.set_primitives_velocity(
          CoordinateVector<>(0., shear_gradient * x, 0.));
      variables.set_primitives_pressure(3.);
      hydro.set_conserved_variables(variables, cell.get_volume());
    }
  }
}

/** @brief Collect original-grid primitive velocities in iteration order. */
static std::vector< CoordinateVector<> >
get_velocities(HydroGridCreator &grid_creator) {
  std::vector< CoordinateVector<> > velocities;
  for (auto grid = grid_creator.begin(); grid != grid_creator.original_end();
       ++grid) {
    for (auto cell = (*grid).hydro_begin(); cell != (*grid).hydro_end();
         ++cell) {
      velocities.push_back(
          cell.get_hydro_variables().get_primitives_velocity());
    }
  }
  return velocities;
}

/** @brief Require two velocity fields to be bitwise identical. */
static void assert_identical(
    const std::vector< CoordinateVector<> > &first,
    const std::vector< CoordinateVector<> > &second) {
  assert_condition(first.size() == second.size());
  for (size_t i = 0; i < first.size(); ++i) {
    assert_condition(first[i].x() == second[i].x());
    assert_condition(first[i].y() == second[i].y());
    assert_condition(first[i].z() == second[i].z());
  }
}

/**
 * @brief Unit test for the AlveliusTurbulenceForcing class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  AlveliusTurbulenceForcing forcing(1, 32, Box<>(0., 1.), 2., 3., 2.5, 0.2, 1.,
                                    42, 1.e-6, 0.);

  double box[6] = {0., 0., 0., 1., 1., 1.};
  CoordinateVector< int_fast32_t > ncell(32, 32, 32);

  HydroDensitySubGrid subgrid(box, ncell);
  for (auto cellit = subgrid.hydro_begin(); cellit != subgrid.hydro_end();
       ++cellit) {
    cellit.get_hydro_variables().conserved(0) = 1.;
  }

  HydroDensitySubGrid subgrid2(box, ncell);
  for (auto cellit = subgrid2.hydro_begin(); cellit != subgrid2.hydro_end();
       ++cellit) {
    cellit.get_hydro_variables().conserved(0) = 1.;
  }

  forcing.update_turbulence(1.e-5);
  forcing.add_turbulent_forcing(0, subgrid);

  std::ofstream ofile("test_alvelius.txt");
  ofile << "#x (m)\ty (m)\tz (m)\tvx (m s^-1)\tvy (m s^-1)\tvz (m s^-1)\n";
  for (auto cellit = subgrid.hydro_begin(); cellit != subgrid.hydro_end();
       ++cellit) {
    const CoordinateVector<> x = cellit.get_cell_midpoint();
    const double m = cellit.get_hydro_variables().get_conserved_mass();
    const CoordinateVector<> v =
        cellit.get_hydro_variables().get_conserved_momentum() / m;
    ofile << x.x() << "\t" << x.y() << "\t" << x.z() << "\t" << v.x() << "\t"
          << v.y() << "\t" << v.z() << "\n";
  }

  forcing.add_turbulent_forcing(0, subgrid2);

  {
    RestartWriter writer("test_alvelius.restart");
    forcing.write_restart_file(writer);
  }

  forcing.update_turbulence(1.e-5);
  forcing.add_turbulent_forcing(0, subgrid);

  {
    RestartReader reader("test_alvelius.restart");
    AlveliusTurbulenceForcing forcing2(reader);
    forcing2.update_turbulence(1.e-5);

    auto cellit1 = subgrid.hydro_begin();
    auto cellit2 = subgrid2.hydro_begin();
    while (cellit1 != subgrid.hydro_end()) {

      const CoordinateVector<> p1 = cellit1.get_cell_midpoint();
      const CoordinateVector<> p2 = cellit2.get_cell_midpoint();
      assert_condition(p1.x() == p2.x());
      assert_condition(p1.y() == p2.y());
      assert_condition(p1.z() == p2.z());

      const CoordinateVector<> v1 =
          cellit1.get_hydro_variables().get_primitives_velocity();
      const CoordinateVector<> v2 =
          cellit2.get_hydro_variables().get_primitives_velocity();
      assert_condition(v1.x() == v2.x());
      assert_condition(v1.y() == v2.y());
      assert_condition(v1.z() == v2.z());

      ++cellit1;
      ++cellit2;
    }
  }

  // One-time turbulence: use a genuinely rectangular 2:2:1 box, with the
  // modes defined on the virtual 2:2:2 cube and sampled at the real cells.
  const Box<> rectangular_box(CoordinateVector<>(-1., -1., -0.5),
                              CoordinateVector<>(2., 2., 1.));
  const CoordinateVector< int_fast32_t > total_cells(8, 8, 4);
  const CoordinateVector< int_fast32_t > subgrids(2, 2, 1);
  const CoordinateVector< bool > periodicity(true, true, false);
  const Abundances abundances;
  const Hydro hydro(5. / 3., 100., 1.e4, 1.e99, false, abundances);
  const double shear_gradient = 3.;
  const double target_rms = 8.5e3;

  HydroGridCreator same_seed_a(rectangular_box, total_cells, subgrids,
                               periodicity);
  HydroGridCreator same_seed_b(rectangular_box, total_cells, subgrids,
                               periodicity);
  HydroGridCreator different_seed(rectangular_box, total_cells, subgrids,
                                  periodicity);
  HydroGridCreator no_shear(rectangular_box, total_cells, subgrids,
                            periodicity);
  HydroGridCreator disabled(rectangular_box, total_cells, subgrids,
                            periodicity);
  initialize_rectangular_grid(same_seed_a, hydro, shear_gradient);
  initialize_rectangular_grid(same_seed_b, hydro, shear_gradient);
  initialize_rectangular_grid(different_seed, hydro, shear_gradient);
  initialize_rectangular_grid(no_shear, hydro, 0.);
  initialize_rectangular_grid(disabled, hydro, shear_gradient);

  const std::vector< CoordinateVector<> > baseline =
      get_velocities(same_seed_a);
  const std::vector< CoordinateVector<> > disabled_before =
      get_velocities(disabled);
  const InitialTurbulenceStatistics stats_a =
      InitialTurbulence::initialize(
          true, same_seed_a, hydro, target_rms, 42);
  const InitialTurbulenceStatistics stats_b =
      InitialTurbulence::initialize(
          true, same_seed_b, hydro, target_rms, 42);
  InitialTurbulence::initialize(
      true, different_seed, hydro, target_rms, 43);
  InitialTurbulence::initialize(
      true, no_shear, hydro, target_rms, 42);
  InitialTurbulence::initialize(
      false, disabled, hydro, target_rms, 42);

  // Disabling the feature is a no-op, and equal seeds are exactly
  // deterministic.
  assert_identical(disabled_before, get_velocities(disabled));
  assert_identical(get_velocities(same_seed_a),
                   get_velocities(same_seed_b));
  assert_values_equal_tol(stats_a.final_mean_velocity.x(), 0., 1.e-8);
  assert_values_equal_tol(stats_a.final_mean_velocity.y(), 0., 1.e-8);
  assert_values_equal_tol(stats_a.final_mean_velocity.z(), 0., 1.e-8);
  assert_values_equal_tol(stats_a.final_rms_velocity, target_rms, 1.e-8);
  assert_condition(stats_a.raw_rms_velocity > 0.);
  assert_condition(stats_a.normalization_factor > 0.);
  assert_values_equal_tol(stats_a.final_rms_velocity,
                          stats_b.final_rms_velocity, 1.e-12);

  // A different seed must change at least one cell.
  const std::vector< CoordinateVector<> > velocity_a =
      get_velocities(same_seed_a);
  const std::vector< CoordinateVector<> > velocity_different =
      get_velocities(different_seed);
  bool found_difference = false;
  for (size_t i = 0; i < velocity_a.size(); ++i) {
    found_difference =
        found_difference || velocity_a[i].x() != velocity_different[i].x() ||
        velocity_a[i].y() != velocity_different[i].y() ||
        velocity_a[i].z() != velocity_different[i].z();
  }
  assert_condition(found_difference);

  // Subtracting the same-seed no-shear realization must recover the original
  // background shear, proving that turbulence was added on top of it.
  const std::vector< CoordinateVector<> > velocity_no_shear =
      get_velocities(no_shear);
  for (size_t i = 0; i < velocity_a.size(); ++i) {
    assert_values_equal_tol(velocity_a[i].x() - velocity_no_shear[i].x(),
                            baseline[i].x(), 1.e-8);
    assert_values_equal_tol(velocity_a[i].y() - velocity_no_shear[i].y(),
                            baseline[i].y(), 1.e-8);
    assert_values_equal_tol(velocity_a[i].z() - velocity_no_shear[i].z(),
                            baseline[i].z(), 1.e-8);
  }

  // Independently recover the mass-weighted statistics from the velocity
  // increment and verify that momentum and total energy were reconstructed.
  double total_mass = 0.;
  double turbulent_second_moment = 0.;
  CoordinateVector<> turbulent_weighted_sum;
  size_t velocity_index = 0;
  for (auto grid = same_seed_a.begin(); grid != same_seed_a.original_end();
       ++grid) {
    for (auto cell = (*grid).hydro_begin(); cell != (*grid).hydro_end();
         ++cell, ++velocity_index) {
      const HydroVariables &variables = cell.get_hydro_variables();
      const double mass = variables.get_conserved_mass();
      const CoordinateVector<> velocity = variables.get_primitives_velocity();
      const CoordinateVector<> turbulent_velocity =
          velocity - baseline[velocity_index];
      total_mass += mass;
      turbulent_weighted_sum += mass * turbulent_velocity;
      turbulent_second_moment += mass * turbulent_velocity.norm2();
      const CoordinateVector<> expected_momentum = mass * velocity;
      assert_values_equal_tol(variables.get_conserved_momentum().x(),
                              expected_momentum.x(), 1.e-10);
      assert_values_equal_tol(variables.get_conserved_momentum().y(),
                              expected_momentum.y(), 1.e-10);
      assert_values_equal_tol(variables.get_conserved_momentum().z(),
                              expected_momentum.z(), 1.e-10);
      const double expected_energy =
          variables.get_primitives_pressure() * cell.get_volume() /
              (5. / 3. - 1.) +
          0.5 * mass * velocity.norm2();
      assert_values_equal_tol(variables.get_conserved_total_energy(),
                              expected_energy, 1.e-10);
    }
  }
  const CoordinateVector<> measured_mean =
      turbulent_weighted_sum / total_mass;
  const double measured_rms = std::sqrt(std::max(
      0., turbulent_second_moment / total_mass - measured_mean.norm2()));
  assert_values_equal_tol(measured_mean.x(), 0., 1.e-8);
  assert_values_equal_tol(measured_mean.y(), 0., 1.e-8);
  assert_values_equal_tol(measured_mean.z(), 0., 1.e-8);
  assert_values_equal_tol(measured_rms, target_rms, 1.e-8);

  return 0;
}
