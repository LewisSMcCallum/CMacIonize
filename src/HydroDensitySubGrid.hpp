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
 * @file HydroDensitySubGrid.hpp
 *
 * @brief Extension of DensitySubGrid that adds hydro variables.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HYDRODENSITYSUBGRID_HPP
#define HYDRODENSITYSUBGRID_HPP

#include "DensitySubGrid.hpp"
#include "DensityValues.hpp"
#include "Hydro.hpp"
#include "HydroVariables.hpp"
#include "IonizationAdvection.hpp"

/**
 * @brief Extension of DensitySubGrid that adds hydro variables.
 */
class HydroDensitySubGrid : public DensitySubGrid {
private:
  /**
   * @brief Build limited ionic slopes and local bounds for one cell.
   *
   * The ordinary symmetric MC stencil is used whenever both neighbours are
   * available. At a local subgrid edge a bounded one-sided three-cell slope is
   * used instead. For the direction normal to a subgrid interface, the cell on
   * the neighbouring subgrid can be supplied as an override, restoring the
   * centred stencil across that interface without storing ion gradients.
   */
  inline void get_ionic_reconstruction_data(
      const int_fast32_t index, const int_fast32_t normal_direction,
      const IonizationVariables *normal_minus_override,
      const IonizationVariables *normal_plus_override,
      double slopes[3][NUMBER_OF_IONNAMES],
      double lower[NUMBER_OF_IONNAMES],
      double upper[NUMBER_OF_IONNAMES]) const {

    const int_fast32_t ix = index / _number_of_cells[3];
    const int_fast32_t remainder = index - ix * _number_of_cells[3];
    const int_fast32_t iy = remainder / _number_of_cells[2];
    const int_fast32_t iz = remainder - iy * _number_of_cells[2];
    const int_fast32_t coordinate[3] = {ix, iy, iz};
    const int_fast32_t stride[3] = {_number_of_cells[3],
                                    _number_of_cells[2], 1};

    const IonizationVariables &cell = _ionization_variables[index];
    for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
      const double x0 = cell.get_ionic_fraction(ion);
      lower[ion] = x0;
      upper[ion] = x0;
    }

    for (int_fast32_t direction = 0; direction < 3; ++direction) {
      const IonizationVariables *minus = nullptr;
      const IonizationVariables *plus = nullptr;
      const IonizationVariables *minus_two = nullptr;
      const IonizationVariables *plus_two = nullptr;

      if (coordinate[direction] > 0) {
        minus = &_ionization_variables[index - stride[direction]];
      } else if (direction == normal_direction) {
        minus = normal_minus_override;
      }
      if (coordinate[direction] + 1 < _number_of_cells[direction]) {
        plus = &_ionization_variables[index + stride[direction]];
      } else if (direction == normal_direction) {
        plus = normal_plus_override;
      }
      if (coordinate[direction] > 1) {
        minus_two = &_ionization_variables[index - 2 * stride[direction]];
      }
      if (coordinate[direction] + 2 < _number_of_cells[direction]) {
        plus_two = &_ionization_variables[index + 2 * stride[direction]];
      }

      for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
        slopes[direction][ion] = IonizationAdvection::limited_cell_slope(
            cell, minus, plus, minus_two, plus_two, ion);
        if (minus != nullptr) {
          lower[ion] = std::min(lower[ion], minus->get_ionic_fraction(ion));
          upper[ion] = std::max(upper[ion], minus->get_ionic_fraction(ion));
        }
        if (plus != nullptr) {
          lower[ion] = std::min(lower[ion], plus->get_ionic_fraction(ion));
          upper[ion] = std::max(upper[ion], plus->get_ionic_fraction(ion));
        }
        if (minus_two != nullptr) {
          lower[ion] =
              std::min(lower[ion], minus_two->get_ionic_fraction(ion));
          upper[ion] =
              std::max(upper[ion], minus_two->get_ionic_fraction(ion));
        }
        if (plus_two != nullptr) {
          lower[ion] =
              std::min(lower[ion], plus_two->get_ionic_fraction(ion));
          upper[ion] =
              std::max(upper[ion], plus_two->get_ionic_fraction(ion));
        }
      }
    }
  }

  /**
   * @brief Run the ordinary hydro Riemann solve and advect ion mass with the
   * exact same final hydro mass flux.
   *
   * Only the upwind donor composition is reconstructed. Its face state uses a
   * full three-dimensional MUSCL-Hancock predictor, and is then projected onto
   * the positivity-preserving ionic flux budget before the equal/opposite
   * conservative tracer flux is deposited.
   */
  inline static void do_reconstructed_flux(
      const Hydro &hydro, const uint_fast8_t direction,
      HydroDensitySubGrid &left_grid, const int_fast32_t index_left,
      HydroDensitySubGrid &right_grid, const int_fast32_t index_right,
      const double dx, const double area, const double dt,
      const bool advect_ionization) {

    HydroVariables &left_hydro = left_grid._hydro_variables[index_left];
    HydroVariables &right_hydro = right_grid._hydro_variables[index_right];
    IonizationVariables &left_ionization =
        left_grid._ionization_variables[index_left];
    IonizationVariables &right_ionization =
        right_grid._ionization_variables[index_right];

    if (!advect_ionization) {
      hydro.do_flux_calculation(direction, left_hydro, left_ionization,
                                right_hydro, right_ionization, dx, area, dt,
                                false);
      return;
    }

    const double left_mass = left_hydro.get_conserved_mass();
    const double right_mass = right_hydro.get_conserved_mass();
    const double left_mass_delta_before = left_hydro.delta_conserved(0);
    const double right_mass_delta_before = right_hydro.delta_conserved(0);

    // Let hydro decide exactly how much gas crosses the face, including the
    // HLLC solve, second-order hydro reconstruction and all flux limiting.
    hydro.do_flux_calculation(direction, left_hydro, left_ionization,
                              right_hydro, right_ionization, dx, area, dt,
                              false);
    const double mflux =
        left_mass_delta_before - left_hydro.delta_conserved(0);
    if (mflux == 0.) {
      return;
    }

    HydroDensitySubGrid &donor_grid = mflux > 0. ? left_grid : right_grid;
    const int_fast32_t donor_index = mflux > 0. ? index_left : index_right;
    IonizationVariables &donor = mflux > 0. ? left_ionization : right_ionization;
    HydroVariables &donor_hydro = mflux > 0. ? left_hydro : right_hydro;
    const double donor_mass = mflux > 0. ? left_mass : right_mass;
    const double donor_mass_delta_before =
        mflux > 0. ? left_mass_delta_before : right_mass_delta_before;

    double slopes[3][NUMBER_OF_IONNAMES];
    double lower[NUMBER_OF_IONNAMES];
    double upper[NUMBER_OF_IONNAMES];
    donor_grid.get_ionic_reconstruction_data(
        donor_index, direction,
        mflux > 0. ? nullptr : &left_ionization,
        mflux > 0. ? &right_ionization : nullptr, slopes, lower, upper);

    const double cell_size[3] = {donor_grid._cell_size[0],
                                 donor_grid._cell_size[1],
                                 donor_grid._cell_size[2]};
    double donor_face[NUMBER_OF_IONNAMES];
    IonizationAdvection::reconstruct_face_3d(
        donor, slopes, lower, upper, donor_hydro.get_primitives_velocity(),
        cell_size, direction, mflux > 0., dt, donor_face);

    IonizationAdvection::limit_outflow_face(
        donor, donor_mass, donor_mass_delta_before, std::abs(mflux), dt,
        donor_face);
    IonizationAdvection::add_donor_interface_flux(
        left_ionization, right_ionization, mflux, donor_face);
  }

  /**
   * @brief Hydro flux plus multidimensional ionic reconstruction at a physical
   * domain boundary.
   */
  inline static void do_reconstructed_ghost_flux(
      const Hydro &hydro, const uint_fast8_t direction,
      const CoordinateVector<> ghost_position, HydroDensitySubGrid &grid,
      const int_fast32_t index, const HydroBoundary &boundary, const double dx,
      const double area, const double dt, const bool advect_ionization) {

    HydroVariables &cell_hydro = grid._hydro_variables[index];
    IonizationVariables &cell_ionization = grid._ionization_variables[index];
    if (!advect_ionization) {
      hydro.do_ghost_flux_calculation(direction, ghost_position, cell_hydro,
                                      cell_ionization, boundary, dx, area, dt,
                                      false);
      return;
    }

    const double old_mass = cell_hydro.get_conserved_mass();
    const double old_mass_delta = cell_hydro.delta_conserved(0);
    hydro.do_ghost_flux_calculation(direction, ghost_position, cell_hydro,
                                    cell_ionization, boundary, dx, area, dt,
                                    false);
    const double mflux = old_mass_delta - cell_hydro.delta_conserved(0);
    if (mflux == 0.) {
      return;
    }

    double slopes[3][NUMBER_OF_IONNAMES];
    double lower[NUMBER_OF_IONNAMES];
    double upper[NUMBER_OF_IONNAMES];
    grid.get_ionic_reconstruction_data(index, direction, nullptr, nullptr,
                                       slopes, lower, upper);
    const double cell_size[3] = {grid._cell_size[0], grid._cell_size[1],
                                 grid._cell_size[2]};
    double face[NUMBER_OF_IONNAMES];
    IonizationAdvection::reconstruct_face_3d(
        cell_ionization, slopes, lower, upper,
        cell_hydro.get_primitives_velocity(), cell_size, direction, dx > 0., dt,
        face);

    // Positive mflux is out of the real cell and therefore needs the donor
    // positivity constraint. Negative mflux enters from the outflow ghost and
    // retains the extrapolated interior composition used historically.
    if (mflux > 0.) {
      IonizationAdvection::limit_outflow_face(
          cell_ionization, old_mass, old_mass_delta, mflux, dt, face);
    }
    IonizationAdvection::add_boundary_flux(cell_ionization, mflux, face);
  }

  /*! @brief Volume of a single cell (in m^3). */
  double _cell_volume;

  /*! @brief Inverse volume of a single cell (in m^-3). */
  double _inverse_cell_volume;

  /*! @brief Surface areas of a single cell (in m^2). */
  double _cell_areas[3];

  /*! @brief Hydrodynamical variables. */
  HydroVariables *_hydro_variables;

  /*! @brief Gradient limiters for the primitive hydrodynamical variables. */
  double *_primitive_variable_limiters;

  /*! @brief Indices of the hydro tasks associated with this subgrid. */
  size_t _hydro_tasks[18];

public:
  /**
   * @brief Constructor.
   *
   * @param box Dimensions of the box that contains the grid (in m; first 3
   * elements are the anchor of the box, 3 last elements are the side lengths
   * of the box).
   * @param ncell Number of cells in each dimension.
   */
  inline HydroDensitySubGrid(const double *box,
                             const CoordinateVector< int_fast32_t > ncell)
      : DensitySubGrid(box, ncell),
        _cell_volume(_cell_size[0] * _cell_size[1] * _cell_size[2]),
        _inverse_cell_volume(1. / _cell_volume),
        _cell_areas{_cell_size[1] * _cell_size[2],
                    _cell_size[0] * _cell_size[2],
                    _cell_size[0] * _cell_size[1]} {

    const int_fast32_t tot_ncell = _number_of_cells[3] * ncell[0];
    _hydro_variables = new HydroVariables[tot_ncell];
    _primitive_variable_limiters = new double[tot_ncell * 10];

    for (int_fast32_t i = 0; i < 5 * tot_ncell; ++i) {
      _primitive_variable_limiters[2 * i] = DBL_MAX;
      _primitive_variable_limiters[2 * i + 1] = -DBL_MAX;
    }
  }

  /** @brief Copy constructor. */
  inline HydroDensitySubGrid(const HydroDensitySubGrid &original)
      : DensitySubGrid(original), _cell_volume(original._cell_volume),
        _inverse_cell_volume(original._inverse_cell_volume),
        _cell_areas{original._cell_areas[0], original._cell_areas[1],
                    original._cell_areas[2]} {

    const int_fast32_t tot_ncell = _number_of_cells[3] * _number_of_cells[0];
    _hydro_variables = new HydroVariables[tot_ncell];
    _primitive_variable_limiters = new double[tot_ncell * 10];

    for (int_fast32_t i = 0; i < tot_ncell; ++i) {
      _hydro_variables[i].copy_all(original._hydro_variables[i]);
    }
    for (int_fast32_t i = 0; i < 10 * tot_ncell; ++i) {
      _primitive_variable_limiters[i] =
          original._primitive_variable_limiters[i];
    }
  }

  /** @brief Destructor. */
  virtual ~HydroDensitySubGrid() {
    delete[] _hydro_variables;
    delete[] _primitive_variable_limiters;
  }

  inline double initialize_hydrodynamic_variables(const Hydro &hydro,
                                                  const bool do_primitives) {
    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    double timestep = DBL_MAX;
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      if (do_primitives) {
        hydro.ionization_to_hydro(_ionization_variables[i],
                                  _hydro_variables[i]);
      }
      hydro.set_conserved_variables(_hydro_variables[i], _cell_volume);
      timestep = std::min(timestep, hydro.get_timestep(_hydro_variables[i],
                                                       _ionization_variables[i],
                                                       _cell_volume));
    }
    return timestep;
  }

  inline void update_conserved_variables(const double timestep,
                                         const bool advect_ionization = false) {
    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      const double old_mass = _hydro_variables[i].get_conserved_mass();
      const CoordinateVector<> a =
          _hydro_variables[i].get_gravitational_acceleration();
      const CoordinateVector<> p =
          _hydro_variables[i].get_conserved_momentum();
      const double mdt = _hydro_variables[i].get_conserved_mass() * timestep;
      _hydro_variables[i].conserved(1) += mdt * a.x();
      _hydro_variables[i].conserved(2) += mdt * a.y();
      _hydro_variables[i].conserved(3) += mdt * a.z();
      _hydro_variables[i].conserved(4) +=
          timestep * CoordinateVector<>::dot_product(p, a);
      _hydro_variables[i].conserved(4) +=
          0.5 * mdt * timestep * a.norm2();
      _hydro_variables[i].conserved(4) +=
          _hydro_variables[i].get_energy_term();
      _hydro_variables[i].set_energy_term(0.);
      for (int_fast8_t j = 0; j < 5; ++j) {
        _hydro_variables[i].conserved(j) +=
            _hydro_variables[i].delta_conserved(j) * timestep;
        _hydro_variables[i].delta_conserved(j) = 0;
        _hydro_variables[i].primitive_gradients(j) = CoordinateVector<>(0.);
        _primitive_variable_limiters[10 * i + 2 * j] = DBL_MAX;
        _primitive_variable_limiters[10 * i + 2 * j + 1] = -DBL_MAX;
      }

      if (advect_ionization) {
        const double new_mass = _hydro_variables[i].get_conserved_mass();
        if (new_mass > 0.) {
          for (int_fast32_t j = 0; j < NUMBER_OF_IONNAMES; ++j) {
            const double old_fraction =
                _ionization_variables[i].get_ionic_fraction(j);
            const double delta_fraction =
                _ionization_variables[i].get_delta_ionic_fraction(j);
            double new_fraction =
                (old_mass * old_fraction + delta_fraction * timestep) /
                new_mass;
            // The face limiter should make these clamps roundoff-only. Keep
            // them as a last-resort guard for specialised boundary paths.
            if (new_fraction < 0.) {
              new_fraction = 0.;
            } else if (new_fraction > 1.) {
              new_fraction = 1.;
            }
            _ionization_variables[i].set_ionic_fraction(j, new_fraction);
          }
          IonizationAdvection::enforce_ionic_simplex(
              _ionization_variables[i]);
        }
        _ionization_variables[i].reset_delta_ionic_fractions();
      }

      cmac_assert(_hydro_variables[i].get_conserved_mass() ==
                  _hydro_variables[i].get_conserved_mass());
      cmac_assert(_hydro_variables[i].get_conserved_momentum().x() ==
                  _hydro_variables[i].get_conserved_momentum().x());
      cmac_assert(_hydro_variables[i].get_conserved_momentum().y() ==
                  _hydro_variables[i].get_conserved_momentum().y());
      cmac_assert(_hydro_variables[i].get_conserved_momentum().z() ==
                  _hydro_variables[i].get_conserved_momentum().z());
      cmac_assert(_hydro_variables[i].get_conserved_total_energy() ==
                  _hydro_variables[i].get_conserved_total_energy());

#ifdef SAFE_HYDRO_VARIABLES
      _hydro_variables[i].conserved(0) =
          std::max(_hydro_variables[i].get_conserved_mass(), 0.);
      _hydro_variables[i].conserved(4) =
          std::max(_hydro_variables[i].get_conserved_total_energy(), 0.);
#else
      cmac_assert(_hydro_variables[i].get_conserved_mass() >= 0.);
      cmac_assert(_hydro_variables[i].get_conserved_total_energy() >= 0.);
#endif
      cmac_assert_message(_hydro_variables[i].get_conserved_mass() > 0.0,
                          "about to set mass = 0");
      cmac_assert_message(
          _hydro_variables[i].get_conserved_total_energy() >
              0.5 * CoordinateVector<>::dot_product(
                        _hydro_variables[i].get_conserved_momentum(),
                        _hydro_variables[i].get_conserved_momentum()) /
                  _hydro_variables[i].get_conserved_mass(),
          "about to set kinetic greater than total energy....");
    }
  }

  inline void update_primitive_variables(const Hydro &hydro) {
    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      hydro.set_primitive_variables(
          _hydro_variables[i], _ionization_variables[i], _inverse_cell_volume);
    }
  }

  inline void
  update_ionization_variables(const Hydro &hydro,
                              const double maximum_neutral_fraction) {
    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      hydro.hydro_to_ionization(_hydro_variables[i], _ionization_variables[i]);
      const double previous_h0 =
          _ionization_variables[i].get_prev_ionic_fraction(ION_H_n);
      const bool has_physical_previous_state =
          previous_h0 >= 0. && previous_h0 <= 1.;
      if (!has_physical_previous_state && maximum_neutral_fraction > 0. &&
          _ionization_variables[i].get_ionic_fraction(ION_H_n) >
              maximum_neutral_fraction) {
        _ionization_variables[i].set_ionic_fraction(
            ION_H_n, maximum_neutral_fraction);
      }
    }
  }

  inline void add_ionization_energy(const Hydro &hydro, const double timestep) {
    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      hydro.add_ionization_energy(_ionization_variables[i],
                                  _hydro_variables[i], _inverse_cell_volume,
                                  timestep);
    }
  }

  inline void predict_primitive_variables(const Hydro &hydro,
                                          const double timestep) {
    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      hydro.predict_primitive_variables(_hydro_variables[i], timestep);
    }
  }

  inline void apply_slope_limiter(const Hydro &hydro) {
    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      hydro.apply_slope_limiter(_hydro_variables[i],
                                &_primitive_variable_limiters[10 * i],
                                _cell_size);
    }
  }

  inline void inner_flux_sweep(const Hydro &hydro, const double dt,
                               const bool advect_ionization = false) {
    for (int_fast32_t ix = 0; ix < _number_of_cells[0] - 1; ++ix) {
      for (int_fast32_t iy = 0; iy < _number_of_cells[1]; ++iy) {
        for (int_fast32_t iz = 0; iz < _number_of_cells[2]; ++iz) {
          const int_fast32_t index000 =
              ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          const int_fast32_t index100 =
              (ix + 1) * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          do_reconstructed_flux(hydro, 0, *this, index000, *this, index100,
                                _cell_size[0], _cell_areas[0], dt,
                                advect_ionization);
        }
      }
    }
    for (int_fast32_t ix = 0; ix < _number_of_cells[0]; ++ix) {
      for (int_fast32_t iy = 0; iy < _number_of_cells[1] - 1; ++iy) {
        for (int_fast32_t iz = 0; iz < _number_of_cells[2]; ++iz) {
          const int_fast32_t index000 =
              ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          const int_fast32_t index010 =
              ix * _number_of_cells[3] + (iy + 1) * _number_of_cells[2] + iz;
          do_reconstructed_flux(hydro, 1, *this, index000, *this, index010,
                                _cell_size[1], _cell_areas[1], dt,
                                advect_ionization);
        }
      }
    }
    for (int_fast32_t ix = 0; ix < _number_of_cells[0]; ++ix) {
      for (int_fast32_t iy = 0; iy < _number_of_cells[1]; ++iy) {
        for (int_fast32_t iz = 0; iz < _number_of_cells[2] - 1; ++iz) {
          const int_fast32_t index000 =
              ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          const int_fast32_t index001 = index000 + 1;
          do_reconstructed_flux(hydro, 2, *this, index000, *this, index001,
                                _cell_size[2], _cell_areas[2], dt,
                                advect_ionization);
        }
      }
    }
  }

  inline void outer_flux_sweep(const int_fast32_t direction, const Hydro &hydro,
                               HydroDensitySubGrid &neighbour,
                               const double dt,
                               const bool advect_ionization = false) {
    int_fast32_t i, start_index_left, start_index_right, row_increment,
        row_length, column_increment, column_length;
    double dx, A;
    HydroDensitySubGrid *left_grid, *right_grid;
    switch (direction) {
    case TRAVELDIRECTION_FACE_X_P:
      i = 0;
      left_grid = this;
      right_grid = &neighbour;
      start_index_left = (_number_of_cells[0] - 1) * _number_of_cells[3];
      start_index_right = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[2];
      column_length = _number_of_cells[1];
      dx = _cell_size[0];
      A = _cell_areas[0];
      break;
    case TRAVELDIRECTION_FACE_X_N:
      i = 0;
      left_grid = &neighbour;
      right_grid = this;
      start_index_left = (_number_of_cells[0] - 1) * _number_of_cells[3];
      start_index_right = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[2];
      column_length = _number_of_cells[1];
      dx = _cell_size[0];
      A = _cell_areas[0];
      break;
    case TRAVELDIRECTION_FACE_Y_P:
      i = 1;
      left_grid = this;
      right_grid = &neighbour;
      start_index_left = (_number_of_cells[1] - 1) * _number_of_cells[2];
      start_index_right = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dx = _cell_size[1];
      A = _cell_areas[1];
      break;
    case TRAVELDIRECTION_FACE_Y_N:
      i = 1;
      left_grid = &neighbour;
      right_grid = this;
      start_index_left = (_number_of_cells[1] - 1) * _number_of_cells[2];
      start_index_right = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dx = _cell_size[1];
      A = _cell_areas[1];
      break;
    case TRAVELDIRECTION_FACE_Z_P:
      i = 2;
      left_grid = this;
      right_grid = &neighbour;
      start_index_left = _number_of_cells[2] - 1;
      start_index_right = 0;
      row_increment = _number_of_cells[2];
      row_length = _number_of_cells[1];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dx = _cell_size[2];
      A = _cell_areas[2];
      break;
    case TRAVELDIRECTION_FACE_Z_N:
      i = 2;
      left_grid = &neighbour;
      right_grid = this;
      start_index_left = _number_of_cells[2] - 1;
      start_index_right = 0;
      row_increment = _number_of_cells[2];
      row_length = _number_of_cells[1];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dx = _cell_size[2];
      A = _cell_areas[2];
      break;
    default:
      cmac_error("Unknown hydro neighbour: %" PRIiFAST32, direction);
      break;
    }

    for (int_fast32_t ic = 0; ic < column_length; ++ic) {
      for (int_fast32_t ir = 0; ir < row_length; ++ir) {
        const int_fast32_t index_left =
            start_index_left + ic * column_increment + ir * row_increment;
        const int_fast32_t index_right =
            start_index_right + ic * column_increment + ir * row_increment;
        do_reconstructed_flux(hydro, i, *left_grid, index_left, *right_grid,
                              index_right, dx, A, dt, advect_ionization);
      }
    }
  }

  inline void outer_ghost_flux_sweep(
      const int_fast32_t direction, const Hydro &hydro,
      const HydroBoundary &boundary, const double dt,
      const bool advect_ionization = false) {
    int_fast32_t i, start_index_left, row_increment, row_length,
        column_increment, column_length;
    double dx, A;
    CoordinateVector<> offset;
    switch (direction) {
    case TRAVELDIRECTION_FACE_X_P:
      i = 0;
      start_index_left = (_number_of_cells[0] - 1) * _number_of_cells[3];
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[2];
      column_length = _number_of_cells[1];
      dx = _cell_size[0];
      A = _cell_areas[0];
      offset = CoordinateVector<>(_cell_size[0], 0., 0.);
      break;
    case TRAVELDIRECTION_FACE_X_N:
      i = 0;
      start_index_left = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[2];
      column_length = _number_of_cells[1];
      dx = -_cell_size[0];
      A = _cell_areas[0];
      offset = CoordinateVector<>(-_cell_size[0], 0., 0.);
      break;
    case TRAVELDIRECTION_FACE_Y_P:
      i = 1;
      start_index_left = (_number_of_cells[1] - 1) * _number_of_cells[2];
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dx = _cell_size[1];
      A = _cell_areas[1];
      offset = CoordinateVector<>(0., _cell_size[1], 0.);
      break;
    case TRAVELDIRECTION_FACE_Y_N:
      i = 1;
      start_index_left = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dx = -_cell_size[1];
      A = _cell_areas[1];
      offset = CoordinateVector<>(0., -_cell_size[1], 0.);
      break;
    case TRAVELDIRECTION_FACE_Z_P:
      i = 2;
      start_index_left = _number_of_cells[2] - 1;
      row_increment = _number_of_cells[2];
      row_length = _number_of_cells[1];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dx = _cell_size[2];
      A = _cell_areas[2];
      offset = CoordinateVector<>(0., 0., _cell_size[2]);
      break;
    case TRAVELDIRECTION_FACE_Z_N:
      i = 2;
      start_index_left = 0;
      row_increment = _number_of_cells[2];
      row_length = _number_of_cells[1];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dx = -_cell_size[2];
      A = _cell_areas[2];
      offset = CoordinateVector<>(0., 0., -_cell_size[2]);
      break;
    default:
      cmac_error("Unknown hydro neighbour: %" PRIiFAST32, direction);
      break;
    }

    for (int_fast32_t ic = 0; ic < column_length; ++ic) {
      for (int_fast32_t ir = 0; ir < row_length; ++ir) {
        const int_fast32_t index_left =
            start_index_left + ic * column_increment + ir * row_increment;
        do_reconstructed_ghost_flux(
            hydro, i, get_cell_midpoint(index_left) + offset, *this, index_left,
            boundary, dx, A, dt, advect_ionization);
      }
    }
  }

  inline void inner_gradient_sweep(const Hydro &hydro) {
    for (int_fast32_t ix = 0; ix < _number_of_cells[0] - 1; ++ix) {
      for (int_fast32_t iy = 0; iy < _number_of_cells[1]; ++iy) {
        for (int_fast32_t iz = 0; iz < _number_of_cells[2]; ++iz) {
          const int_fast32_t index000 =
              ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          const int_fast32_t index100 =
              (ix + 1) * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          hydro.do_gradient_calculation(
              0, _hydro_variables[index000], _hydro_variables[index100],
              _inv_cell_size[0], &_primitive_variable_limiters[10 * index000],
              &_primitive_variable_limiters[10 * index100]);
        }
      }
    }
    for (int_fast32_t ix = 0; ix < _number_of_cells[0]; ++ix) {
      for (int_fast32_t iy = 0; iy < _number_of_cells[1] - 1; ++iy) {
        for (int_fast32_t iz = 0; iz < _number_of_cells[2]; ++iz) {
          const int_fast32_t index000 =
              ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          const int_fast32_t index010 =
              ix * _number_of_cells[3] + (iy + 1) * _number_of_cells[2] + iz;
          hydro.do_gradient_calculation(
              1, _hydro_variables[index000], _hydro_variables[index010],
              _inv_cell_size[1], &_primitive_variable_limiters[10 * index000],
              &_primitive_variable_limiters[10 * index010]);
        }
      }
    }
    for (int_fast32_t ix = 0; ix < _number_of_cells[0]; ++ix) {
      for (int_fast32_t iy = 0; iy < _number_of_cells[1]; ++iy) {
        for (int_fast32_t iz = 0; iz < _number_of_cells[2] - 1; ++iz) {
          const int_fast32_t index000 =
              ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          const int_fast32_t index001 =
              ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz + 1;
          hydro.do_gradient_calculation(
              2, _hydro_variables[index000], _hydro_variables[index001],
              _inv_cell_size[2], &_primitive_variable_limiters[10 * index000],
              &_primitive_variable_limiters[10 * index001]);
        }
      }
    }
  }

  inline void outer_gradient_sweep(const int_fast32_t direction,
                                   const Hydro &hydro,
                                   HydroDensitySubGrid &neighbour) {
    int_fast32_t i, start_index_left, start_index_right, row_increment,
        row_length, column_increment, column_length;
    double dxinv;
    HydroDensitySubGrid *left_grid, *right_grid;
    switch (direction) {
    case TRAVELDIRECTION_FACE_X_P:
      i = 0;
      left_grid = this;
      right_grid = &neighbour;
      start_index_left = (_number_of_cells[0] - 1) * _number_of_cells[3];
      start_index_right = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[2];
      column_length = _number_of_cells[1];
      dxinv = _inv_cell_size[0];
      break;
    case TRAVELDIRECTION_FACE_X_N:
      i = 0;
      left_grid = &neighbour;
      right_grid = this;
      start_index_left = (_number_of_cells[0] - 1) * _number_of_cells[3];
      start_index_right = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[2];
      column_length = _number_of_cells[1];
      dxinv = _inv_cell_size[0];
      break;
    case TRAVELDIRECTION_FACE_Y_P:
      i = 1;
      left_grid = this;
      right_grid = &neighbour;
      start_index_left = (_number_of_cells[1] - 1) * _number_of_cells[2];
      start_index_right = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dxinv = _inv_cell_size[1];
      break;
    case TRAVELDIRECTION_FACE_Y_N:
      i = 1;
      left_grid = &neighbour;
      right_grid = this;
      start_index_left = (_number_of_cells[1] - 1) * _number_of_cells[2];
      start_index_right = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dxinv = _inv_cell_size[1];
      break;
    case TRAVELDIRECTION_FACE_Z_P:
      i = 2;
      left_grid = this;
      right_grid = &neighbour;
      start_index_left = _number_of_cells[2] - 1;
      start_index_right = 0;
      row_increment = _number_of_cells[2];
      row_length = _number_of_cells[1];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dxinv = _inv_cell_size[2];
      break;
    case TRAVELDIRECTION_FACE_Z_N:
      i = 2;
      left_grid = &neighbour;
      right_grid = this;
      start_index_left = _number_of_cells[2] - 1;
      start_index_right = 0;
      row_increment = _number_of_cells[2];
      row_length = _number_of_cells[1];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dxinv = _inv_cell_size[2];
      break;
    default:
      cmac_error("Unknown hydro neighbour: %" PRIiFAST32, direction);
      break;
    }

    for (int_fast32_t ic = 0; ic < column_length; ++ic) {
      for (int_fast32_t ir = 0; ir < row_length; ++ir) {
        const int_fast32_t index_left =
            start_index_left + ic * column_increment + ir * row_increment;
        const int_fast32_t index_right =
            start_index_right + ic * column_increment + ir * row_increment;
        hydro.do_gradient_calculation(
            i, left_grid->_hydro_variables[index_left],
            right_grid->_hydro_variables[index_right], dxinv,
            &left_grid->_primitive_variable_limiters[10 * index_left],
            &right_grid->_primitive_variable_limiters[10 * index_right]);
      }
    }
  }

  inline void add_x_remapped_ghost_gradient(const uint_fast32_t index,
                                            const Hydro &hydro,
                                            HydroVariables ghost,
                                            const bool ghost_on_right) {
    double ghost_limiters[10];
    for (uint_fast8_t i = 0; i < 5; ++i) {
      ghost_limiters[2 * i] = DBL_MAX;
      ghost_limiters[2 * i + 1] = -DBL_MAX;
    }
    if (ghost_on_right) {
      hydro.do_gradient_calculation(
          0, _hydro_variables[index], ghost, _inv_cell_size[0],
          &_primitive_variable_limiters[10 * index], ghost_limiters);
    } else {
      hydro.do_gradient_calculation(
          0, ghost, _hydro_variables[index], _inv_cell_size[0],
          ghost_limiters, &_primitive_variable_limiters[10 * index]);
    }
  }

  inline void outer_ghost_gradient_sweep(const int_fast32_t direction,
                                         const Hydro &hydro,
                                         const HydroBoundary &boundary) {
    int_fast32_t i, start_index_left, row_increment, row_length,
        column_increment, column_length;
    double dxinv;
    HydroDensitySubGrid *left_grid;
    CoordinateVector<> offset;
    switch (direction) {
    case TRAVELDIRECTION_FACE_X_P:
      i = 0;
      left_grid = this;
      start_index_left = (_number_of_cells[0] - 1) * _number_of_cells[3];
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[2];
      column_length = _number_of_cells[1];
      dxinv = _inv_cell_size[0];
      offset = CoordinateVector<>(_cell_size[0], 0., 0.);
      break;
    case TRAVELDIRECTION_FACE_X_N:
      i = 0;
      left_grid = this;
      start_index_left = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[2];
      column_length = _number_of_cells[1];
      dxinv = -_inv_cell_size[0];
      offset = CoordinateVector<>(-_cell_size[0], 0., 0.);
      break;
    case TRAVELDIRECTION_FACE_Y_P:
      i = 1;
      left_grid = this;
      start_index_left = (_number_of_cells[1] - 1) * _number_of_cells[2];
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dxinv = _inv_cell_size[1];
      offset = CoordinateVector<>(0., _cell_size[1], 0.);
      break;
    case TRAVELDIRECTION_FACE_Y_N:
      i = 1;
      left_grid = this;
      start_index_left = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dxinv = -_inv_cell_size[1];
      offset = CoordinateVector<>(0., -_cell_size[1], 0.);
      break;
    case TRAVELDIRECTION_FACE_Z_P:
      i = 2;
      left_grid = this;
      start_index_left = _number_of_cells[2] - 1;
      row_increment = _number_of_cells[2];
      row_length = _number_of_cells[1];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dxinv = _inv_cell_size[2];
      offset = CoordinateVector<>(0., 0., _cell_size[2]);
      break;
    case TRAVELDIRECTION_FACE_Z_N:
      i = 2;
      left_grid = this;
      start_index_left = 0;
      row_increment = _number_of_cells[2];
      row_length = _number_of_cells[1];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dxinv = -_inv_cell_size[2];
      offset = CoordinateVector<>(0., 0., -_cell_size[2]);
      break;
    default:
      cmac_error("Unknown hydro neighbour: %" PRIiFAST32, direction);
      break;
    }

    for (int_fast32_t ic = 0; ic < column_length; ++ic) {
      for (int_fast32_t ir = 0; ir < row_length; ++ir) {
        const int_fast32_t index_left =
            start_index_left + ic * column_increment + ir * row_increment;
        hydro.do_ghost_gradient_calculation(
            i, get_cell_midpoint(index_left) + offset,
            left_grid->_hydro_variables[index_left], boundary, dxinv,
            &left_grid->_primitive_variable_limiters[10 * index_left]);
      }
    }
  }

  inline void set_hydro_task(const int_fast32_t i, const size_t task) {
    _hydro_tasks[i] = task;
  }

  inline size_t get_hydro_task(const int_fast32_t i) const {
    return _hydro_tasks[i];
  }

  virtual void initialize_hydro(const uint_fast32_t index,
                                const DensityValues &values) {
    _hydro_variables[index].set_primitives_velocity(values.get_velocity());
  }

  class hydroiterator : public Cell {
  private:
    uint_fast32_t _index;
    HydroDensitySubGrid *_subgrid;

  public:
    inline hydroiterator(const uint_fast32_t index,
                         HydroDensitySubGrid &subgrid)
        : _index(index), _subgrid(&subgrid) {}

    virtual CoordinateVector<> get_cell_midpoint() const {
      return _subgrid->get_cell_midpoint(_index);
    }

    virtual double get_volume() const { return _subgrid->_cell_volume; }

    virtual std::vector< Face > get_faces() const {
      return std::vector< Face >();
    }

    inline const HydroVariables &get_hydro_variables() const {
      return _subgrid->_hydro_variables[_index];
    }

    inline HydroVariables &get_hydro_variables() {
      return _subgrid->_hydro_variables[_index];
    }

    inline const IonizationVariables &get_ionization_variables() const {
      return _subgrid->_ionization_variables[_index];
    }

    inline IonizationVariables &get_ionization_variables() {
      return _subgrid->_ionization_variables[_index];
    }

    inline hydroiterator &operator++() {
      ++_index;
      return *this;
    }

    inline hydroiterator &operator+=(const uint_fast32_t increment) {
      _index += increment;
      return *this;
    }

    inline hydroiterator operator+(const uint_fast32_t increment) const {
      hydroiterator it(*this);
      it += increment;
      return it;
    }

    inline uint_fast32_t get_index() const { return _index; }

    inline bool operator==(hydroiterator it) const {
      return (_subgrid == it._subgrid && _index == it._index);
    }

    inline bool operator!=(hydroiterator it) const { return !(*this == it); }
  };

  inline hydroiterator hydro_begin() { return hydroiterator(0, *this); }

  inline hydroiterator hydro_end() {
    return hydroiterator(
        _number_of_cells[0] * _number_of_cells[1] * _number_of_cells[2], *this);
  }

  inline hydroiterator get_hydro_cell(const CoordinateVector<> position) {
    CoordinateVector< int_fast32_t > three_index;
    return hydroiterator(get_start_index(position - _anchor,
                                         TRAVELDIRECTION_INSIDE, three_index),
                         *this);
  }

  virtual void write_restart_file(RestartWriter &restart_writer) const {
    DensitySubGrid::write_restart_file(restart_writer);
    restart_writer.write(_cell_volume);
    restart_writer.write(_cell_areas[0]);
    restart_writer.write(_cell_areas[1]);
    restart_writer.write(_cell_areas[2]);
    const int_fast32_t number_of_cells =
        _number_of_cells[0] * _number_of_cells[1] * _number_of_cells[2];
    for (int_fast32_t i = 0; i < number_of_cells; ++i) {
      _hydro_variables[i].write_restart_file(restart_writer);
      _ionization_variables[i].write_restart_file(restart_writer);
    }
  }

  inline HydroDensitySubGrid(RestartReader &restart_reader)
      : DensitySubGrid(restart_reader) {
    _cell_volume = restart_reader.read< double >();
    _inverse_cell_volume = 1. / _cell_volume;
    _cell_areas[0] = restart_reader.read< double >();
    _cell_areas[1] = restart_reader.read< double >();
    _cell_areas[2] = restart_reader.read< double >();
    const int_fast32_t number_of_cells =
        _number_of_cells[0] * _number_of_cells[1] * _number_of_cells[2];
    _hydro_variables = new HydroVariables[number_of_cells];
    for (int_fast32_t i = 0; i < number_of_cells; ++i) {
      _hydro_variables[i] = HydroVariables(restart_reader);
      _ionization_variables[i] = IonizationVariables(restart_reader);
    }
    _primitive_variable_limiters = new double[10 * number_of_cells];
    for (int_fast32_t i = 0; i < 5 * number_of_cells; ++i) {
      _primitive_variable_limiters[2 * i] = DBL_MAX;
      _primitive_variable_limiters[2 * i + 1] = -DBL_MAX;
    }
  }
};

#endif // HYDRODENSITYSUBGRID_HPP
