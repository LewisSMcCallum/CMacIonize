/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2026 Lewis McCallum
 *
 * CMacIonize is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 ******************************************************************************/

#include "Assert.hpp"
#include "GalacticShearingBox.hpp"
#include "HomogeneousDensityFunction.hpp"

int main(int argc, char **argv) {
  ParameterFile params("test_galactic_shearing_box.param");
  GalacticShearingBox shearing_box(params);

  const Box<> box(CoordinateVector<>(-2., -4., -0.5),
                  CoordinateVector<>(4., 8., 1.));
  DensitySubGridCreator< HydroDensitySubGrid > grid_creator(
      box, CoordinateVector< int_fast32_t >(4, 8, 1),
      CoordinateVector< int_fast32_t >(2, 2, 1),
      CoordinateVector< bool >(true, true, false));
  HomogeneousDensityFunction density_function;
  grid_creator.initialize(density_function);
  shearing_box.validate_boundaries(CoordinateVector< bool >(true, true, false),
                                   grid_creator);

  const Abundances abundances;
  const Hydro hydro(5. / 3., 100., 1.e4, 1.e99, false, abundances);

  // A uniform state on top of v_y = Omega*x must not see the edge-to-edge
  // tangential velocity jump after the frame transformation.
  for (auto grid = grid_creator.begin(); grid != grid_creator.original_end();
       ++grid) {
    for (auto cell = (*grid).hydro_begin(); cell != (*grid).hydro_end();
         ++cell) {
      HydroVariables &variables = cell.get_hydro_variables();
      variables.set_primitives_density(1.);
      variables.set_primitives_velocity(
          CoordinateVector<>(0., cell.get_cell_midpoint().x(), 0.));
      variables.set_primitives_pressure(1.);
      cell.get_ionization_variables().set_ionic_fraction(ION_H_n, 0.5);
    }
    (*grid).initialize_hydrodynamic_variables(hydro, false);
  }
  shearing_box.add_boundary_fluxes(grid_creator, hydro, 0.375, 0.01, true);
  for (auto grid = grid_creator.begin(); grid != grid_creator.original_end();
       ++grid) {
    for (auto cell = (*grid).hydro_begin(); cell != (*grid).hydro_end();
         ++cell) {
      const HydroVariables &variables = cell.get_hydro_variables();
      assert_values_equal_tol(variables.delta_conserved(0), 0., 1.e-12);
      assert_values_equal_tol(variables.delta_conserved(2), 0., 1.e-12);
      assert_values_equal_tol(variables.delta_conserved(4), 0., 1.e-12);
      assert_values_equal_tol(
          cell.get_ionization_variables().get_delta_ionic_fraction(ION_H_n),
          0., 1.e-12);
    }
  }

  // Reset and use a non-uniform state with radial flow. The remapped face
  // exchange must conserve mass, x/z momentum and advected ion mass.
  int_fast32_t row = 0;
  for (auto grid = grid_creator.begin(); grid != grid_creator.original_end();
       ++grid) {
    for (auto cell = (*grid).hydro_begin(); cell != (*grid).hydro_end();
         ++cell) {
      HydroVariables &variables = cell.get_hydro_variables();
      for (uint_fast8_t i = 0; i < 5; ++i) {
        variables.delta_conserved(i) = 0.;
        variables.primitive_gradients(i) = CoordinateVector<>(0.);
      }
      cell.get_ionization_variables().reset_delta_ionic_fractions();
      const double y = cell.get_cell_midpoint().y();
      variables.set_primitives_density(1.2 + 0.03 * y);
      variables.set_primitives_velocity(
          CoordinateVector<>(0.15 + 0.01 * y,
                             cell.get_cell_midpoint().x(), 0.04 * y));
      variables.set_primitives_pressure(2. + 0.02 * y);
      cell.get_ionization_variables().set_ionic_fraction(
          ION_H_n, 0.4 + 0.01 * row);
      ++row;
    }
    (*grid).initialize_hydrodynamic_variables(hydro, false);
  }
  shearing_box.add_boundary_fluxes(grid_creator, hydro, 0.375, 0.01, true);

  double total_mass = 0.;
  double total_x_momentum = 0.;
  double total_z_momentum = 0.;
  double total_ion_mass = 0.;
  for (auto grid = grid_creator.begin(); grid != grid_creator.original_end();
       ++grid) {
    for (auto cell = (*grid).hydro_begin(); cell != (*grid).hydro_end();
         ++cell) {
      const HydroVariables &variables = cell.get_hydro_variables();
      total_mass += variables.delta_conserved(0);
      total_x_momentum += variables.delta_conserved(1);
      total_z_momentum += variables.delta_conserved(3);
      total_ion_mass +=
          cell.get_ionization_variables().get_delta_ionic_fraction(ION_H_n);
    }
  }
  assert_values_equal_tol(total_mass, 0., 1.e-12);
  assert_values_equal_tol(total_x_momentum, 0., 1.e-12);
  assert_values_equal_tol(total_z_momentum, 0., 1.e-12);
  assert_values_equal_tol(total_ion_mass, 0., 1.e-12);

  return 0;
}
