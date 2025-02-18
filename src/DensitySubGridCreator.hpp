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
 * @file DensitySubGridCreator.hpp
 *
 * @brief Class responsible for creating DensitySubGrid instances that make up
 * a larger grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYSUBGRIDCREATOR_HPP
#define DENSITYSUBGRIDCREATOR_HPP

#include "Box.hpp"
#include "DensityFunction.hpp"
#include "DensitySubGrid.hpp"
#include "Error.hpp"
#include "OpenMP.hpp"
#include "ParameterFile.hpp"
#include "HydroDensitySubGrid.hpp"

#include <cinttypes>
#include <vector>
#include <cmath>

/**
 * @brief Class responsible for creating DensitySubGrid instances that make up
 * a larger grid.
 */
template < class _subgrid_type_ > class DensitySubGridCreator {
private:
  /*! @brief Subgrids that make up the grid. */
  std::vector< _subgrid_type_ * > _subgrids;

  /*! @brief Indices of the original subgrid corresponding to each copy. */
  std::vector< size_t > _originals;

  /*! @brief Indices of the first copy of each subgrid. */
  std::vector< size_t > _copies;

  /*! @brief Dimensions of the simulation box (in m). */
  const Box<> _box;

  /*! @brief Side length of a single subgrid (in m). */
  const CoordinateVector<> _subgrid_sides;

  /*! @brief Number of subgrids in each coordinate direction. */
  const CoordinateVector< int_fast32_t > _number_of_subgrids;

  /*! @brief Number of cells in each coordinate direction for a single
   *  subgrid. */
  const CoordinateVector< int_fast32_t > _subgrid_number_of_cells;

  /*! @brief Periodicity flags. */
  const CoordinateVector< bool > _periodicity;

public:
  /**
   * @brief Constructor.
   *
   * @param box Dimensions of the simulation box (in m).
   * @param number_of_cells Number of cells in each coordinate direction.
   * @param number_of_subgrids Number of subgrids in each coordinate direction.
   * @param periodicity Periodicity flags.
   */
  inline DensitySubGridCreator(
      const Box<> box, const CoordinateVector< int_fast32_t > number_of_cells,
      const CoordinateVector< int_fast32_t > number_of_subgrids,
      const CoordinateVector< bool > periodicity)
      : _box(box), _subgrid_sides(box.get_sides()[0] / number_of_subgrids[0],
                                  box.get_sides()[1] / number_of_subgrids[1],
                                  box.get_sides()[2] / number_of_subgrids[2]),
        _number_of_subgrids(number_of_subgrids[0], number_of_subgrids[1],
                            number_of_subgrids[2]),
        _subgrid_number_of_cells(number_of_cells[0] / number_of_subgrids[0],
                                 number_of_cells[1] / number_of_subgrids[1],
                                 number_of_cells[2] / number_of_subgrids[2]),
        _periodicity(periodicity) {

    for (uint_fast8_t i = 0; i < 3; ++i) {
      if (number_of_cells[i] % number_of_subgrids[i] != 0) {
        cmac_error("Number of subgrids not compatible with number of cells!");
      }
    }

    _subgrids.resize(_number_of_subgrids[0] * _number_of_subgrids[1] *
                         _number_of_subgrids[2],
                     nullptr);
    _copies.resize(_subgrids.size(), 0xffffffff);
  }

  /**
   * @brief ParameterFile constructor.
   *
   * This method will read the following parameters from the parameter file:
   *  - (DensityGrid:)number of cells: number of cells in each coordinate
   *    direction (default: [64, 64, 64])
   *  - number of subgrids: number of subgrids in each coordinate direction
   *    (default: [8, 8, 8])
   *
   * @param box Dimensions of the simulation box (in m).
   * @param params ParameterFile to read from.
   */
  inline DensitySubGridCreator(const Box<> box, ParameterFile &params)
      : DensitySubGridCreator(
            box,
            params.get_value< CoordinateVector< int_fast32_t > >(
                "DensityGrid:number of cells",
                CoordinateVector< int_fast32_t >(64)),
            params.get_value< CoordinateVector< int_fast32_t > >(
                "DensitySubGridCreator:number of subgrids",
                CoordinateVector< int_fast32_t >(8)),
            params.get_value< CoordinateVector< bool > >(
                "DensitySubGridCreator:periodicity",
                CoordinateVector< bool >(false))) {}

  /**
   * @brief Destructor.
   */
  inline ~DensitySubGridCreator() {
    for (uint_fast32_t igrid = 0; igrid < _subgrids.size(); ++igrid) {
      delete _subgrids[igrid];
    }
  }

  /**
   * @brief Get the number of subgrids created by the creator.
   *
   * @return Total number of original subgrids.
   */
  inline uint_fast32_t number_of_original_subgrids() const {
    return _number_of_subgrids.x() * _number_of_subgrids.y() *
           _number_of_subgrids.z();
  }

  /**
   * @brief Get the actual number of subgrids including all copies.
   *
   * @return Total number of subgrids.
   */
  inline uint_fast32_t number_of_actual_subgrids() const {
    return _subgrids.size();
  }

  /**
   * @brief Total number of cells in the grid.
   *
   * @return Total number of cells.
   */
  inline uint_fast64_t number_of_cells() const {
    return number_of_original_subgrids() * _subgrid_number_of_cells.x() *
           _subgrid_number_of_cells.y() * _subgrid_number_of_cells.z();
  }

  /**
   * @brief Get the number of subgrids in each coordinate direction.
   *
   * @return Number of subgrids in each coordinate direction.
   */
  inline CoordinateVector< int_fast32_t > get_subgrid_layout() const {
    return _number_of_subgrids;
  }

  /**
   * @brief Get the number of cells in each coordinate direction per subgrid.
   *
   * @return Number of cells in each coordinate direction per subgrid.
   */
  inline CoordinateVector< int_fast32_t > get_subgrid_cell_layout() const {
    return _subgrid_number_of_cells;
  }

  /**
   * @brief Get the dimensions of the box containing the grid.
   *
   * @return Dimensions of the box containing the grid (in m).
   */
  inline Box<> get_box() const { return _box; }

  /**
   * @brief Get the 3D integer coordinates of the given subgrid index within
   * the subgrid grid layout.
   *
   * @param index Subgrid index (needs to be smaller than number_of_subgrids).
   * @return 3D integer coordinates of the subgrid within the subgrid grid
   * layout.
   */
  inline CoordinateVector< int_fast32_t >
  get_grid_position(const size_t index) const {

    const int_fast32_t ix =
        index / (_number_of_subgrids[1] * _number_of_subgrids[2]);
    const int_fast32_t iy =
        (index - ix * _number_of_subgrids[1] * _number_of_subgrids[2]) /
        _number_of_subgrids[2];
    const int_fast32_t iz =
        index - ix * _number_of_subgrids[1] * _number_of_subgrids[2] -
        iy * _number_of_subgrids[2];
    return CoordinateVector< int_fast32_t >(ix, iy, iz);
  }

  /**
   * @brief Get the indices for the neighbouring subgrids of the given subgrid.
   *
   * @param index Subgrid index (needs to be smaller than number_of_subgrids).
   * @param neighbours Return array containing the indices of the neighbouring
   * subgrids.
   * @return Number of existing neighbours.
   */
  inline uint_fast8_t get_neighbours(const size_t index,
                                     size_t neighbours[6]) const {

    const CoordinateVector< int_fast32_t > p = get_grid_position(index);
    uint_fast8_t number_of_neighbours = 0;
    if (p.x() > 0) {
      const size_t ngbi =
          (p.x() - 1) * _number_of_subgrids[1] * _number_of_subgrids[2] +
          p.y() * _number_of_subgrids[2] + p.z();
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    } else if (_periodicity.x()) {
      const size_t ngbi = (_number_of_subgrids[0] - 1) *
                              _number_of_subgrids[1] * _number_of_subgrids[2] +
                          p.y() * _number_of_subgrids[2] + p.z();
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    }
    if (p.x() < _number_of_subgrids[0] - 1) {
      const size_t ngbi =
          (p.x() + 1) * _number_of_subgrids[1] * _number_of_subgrids[2] +
          p.y() * _number_of_subgrids[2] + p.z();
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    } else if (_periodicity.x()) {
      const size_t ngbi = p.y() * _number_of_subgrids[2] + p.z();
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    }
    if (p.y() > 0) {
      const size_t ngbi =
          p.x() * _number_of_subgrids[1] * _number_of_subgrids[2] +
          (p.y() - 1) * _number_of_subgrids[2] + p.z();
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    } else if (_periodicity.y()) {
      const size_t ngbi =
          p.x() * _number_of_subgrids[1] * _number_of_subgrids[2] +
          (_number_of_subgrids[1] - 1) * _number_of_subgrids[2] + p.z();
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    }
    if (p.y() < _number_of_subgrids[1] - 1) {
      const size_t ngbi =
          p.x() * _number_of_subgrids[1] * _number_of_subgrids[2] +
          (p.y() + 1) * _number_of_subgrids[2] + p.z();
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    } else if (_periodicity.y()) {
      const size_t ngbi =
          p.x() * _number_of_subgrids[1] * _number_of_subgrids[2] + p.z();
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    }
    if (p.z() > 0) {
      const size_t ngbi =
          p.x() * _number_of_subgrids[1] * _number_of_subgrids[2] +
          p.y() * _number_of_subgrids[2] + p.z() - 1;
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    } else if (_periodicity.z()) {
      const size_t ngbi =
          p.x() * _number_of_subgrids[1] * _number_of_subgrids[2] +
          p.y() * _number_of_subgrids[2] + _number_of_subgrids[2] - 1;
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    }
    if (p.z() < _number_of_subgrids[2] - 1) {
      const size_t ngbi =
          p.x() * _number_of_subgrids[1] * _number_of_subgrids[2] +
          p.y() * _number_of_subgrids[2] + p.z() + 1;
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    } else if (_periodicity.z()) {
      const size_t ngbi =
          p.x() * _number_of_subgrids[1] * _number_of_subgrids[2] +
          p.y() * _number_of_subgrids[2];
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    }
    return number_of_neighbours;
  }

  /**
   * @brief Create the DensitySubGrid with the given index.
   *
   * @param index Index of the subgrid.
   * @return Pointer to a newly created DensitySubGrid instance. Memory
   * management for the pointer is transferred to the caller.
   */
  inline _subgrid_type_ *create_subgrid(const uint_fast32_t index) const {

    const int_fast32_t ix =
        index / (_number_of_subgrids[1] * _number_of_subgrids[2]);
    const int_fast32_t iy =
        (index - ix * _number_of_subgrids[1] * _number_of_subgrids[2]) /
        _number_of_subgrids[2];
    const int_fast32_t iz =
        index - ix * _number_of_subgrids[1] * _number_of_subgrids[2] -
        iy * _number_of_subgrids[2];
    const double subgrid_box[6] = {
        _box.get_anchor()[0] + ix * _subgrid_sides[0],
        _box.get_anchor()[1] + iy * _subgrid_sides[1],
        _box.get_anchor()[2] + iz * _subgrid_sides[2],
        _subgrid_sides[0],
        _subgrid_sides[1],
        _subgrid_sides[2]};
    _subgrid_type_ *this_grid =
        new _subgrid_type_(subgrid_box, _subgrid_number_of_cells);
    for (int_fast32_t i = 0; i < TRAVELDIRECTION_NUMBER; ++i) {
      this_grid->set_neighbour(i, NEIGHBOUR_OUTSIDE);
      this_grid->set_active_buffer(i, NEIGHBOUR_OUTSIDE);
    }
    for (int_fast32_t nix = -1; nix < 2; ++nix) {
      for (int_fast32_t niy = -1; niy < 2; ++niy) {
        for (int_fast32_t niz = -1; niz < 2; ++niz) {
          // get neighbour corrected indices
          int_fast32_t cix = ix + nix;
          int_fast32_t ciy = iy + niy;
          int_fast32_t ciz = iz + niz;
          if (_periodicity.x()) {
            if (cix < 0) {
              cix = _number_of_subgrids[0] - 1;
            }
            if (cix >= _number_of_subgrids[0]) {
              cix = 0;
            }
          }
          if (_periodicity.y()) {
            if (ciy < 0) {
              ciy = _number_of_subgrids[1] - 1;
            }
            if (ciy >= _number_of_subgrids[1]) {
              ciy = 0;
            }
          }
          if (_periodicity.z()) {
            if (ciz < 0) {
              ciz = _number_of_subgrids[2] - 1;
            }
            if (ciz >= _number_of_subgrids[2]) {
              ciz = 0;
            }
          }
          // if the indices above point to a real subgrid: set up the
          // neighbour relations
          if ((cix >= 0 && cix < _number_of_subgrids[0]) &&
              (ciy >= 0 && ciy < _number_of_subgrids[1]) &&
              (ciz >= 0 && ciz < _number_of_subgrids[2])) {
            // we use get_output_direction() to get the correct index
            // for the neighbour
            // the three_index components will either be
            //  - -ncell --> negative --> lower limit
            //  - 0 --> in range --> inside
            //  - ncell --> upper limit
            const CoordinateVector< int_fast32_t > three_index(
                nix * _subgrid_number_of_cells[0],
                niy * _subgrid_number_of_cells[1],
                niz * _subgrid_number_of_cells[2]);
            const int_fast32_t ngbi =
                this_grid->get_output_direction(three_index);
            // now get the actual ngb index
            const uint_fast32_t ngb_index =
                cix * _number_of_subgrids[1] * _number_of_subgrids[2] +
                ciy * _number_of_subgrids[2] + ciz;
            this_grid->set_neighbour(ngbi, ngb_index);
          } // if ci
        }   // for niz
      }     // for niy
    }       // for nix

    return this_grid;
  }

  /**
   * @brief Initialize the subgrids that make up the grid.
   *
   * @param density_function DensityFunction to use to initialize the cell
   * variables.
   */
  inline void initialize(DensityFunction &density_function) {
    AtomicValue< size_t > igrid(0);
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
    while (igrid.value() < _subgrids.size()) {
      const size_t this_igrid = igrid.post_increment();
      if (this_igrid < _subgrids.size()) {
        _subgrids[this_igrid] = create_subgrid(this_igrid);
        _subgrids[this_igrid]->set_owning_thread(get_thread_index());
        for (auto it = _subgrids[this_igrid]->begin();
             it != _subgrids[this_igrid]->end(); ++it) {
          DensityValues values = density_function(it);
          it.get_ionization_variables().set_number_density(
              values.get_number_density());
          it.get_ionization_variables().set_dust_density(
              values.get_dust_gas_ratio()*values.get_number_density()*1.67e-27);
          // it.get_ionization_variables().set_fraction_silicon(
          //     values.get_fraction_silicates());

          for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
            it.get_ionization_variables().set_ionic_fraction(
                ion, values.get_ionic_fraction(ion));
          }
          it.get_ionization_variables().set_temperature(
              values.get_temperature());
          _subgrids[this_igrid]->initialize_hydro(it.get_index(), values);
        }
      }
    }
  }

  /**
   * @brief Create copies for the given subgrids according to the given copy
   * level specification.
   *
   * @param copy_levels Desired copy level for each subgrid.
   */
  inline void create_copies(std::vector< uint_fast8_t > &copy_levels) {

    // we need to store the original number of subgrids for reference
    const int_fast32_t number_of_unique_subgrids =
        number_of_original_subgrids();
    cmac_assert_message(number_of_unique_subgrids ==
                            _number_of_subgrids[0] * _number_of_subgrids[1] *
                                _number_of_subgrids[2],
                        "Number of subgrids does not match expectation!");
    // we need to do 2 loops:
    //  - one loop to create the copies and store the offset of the first copy
    //    for each subgrid
    //  - a second loop that sets the neighbours (and has access to all
    //  necessary
    //    copies to set inter-copy neighbour relations)

    // array to store the offsets of new copies in
    for (int_fast32_t i = 0; i < number_of_unique_subgrids; ++i) {
      const uint_fast8_t level = copy_levels[i];
      const uint_fast32_t number_of_copies = 1 << level;
      // create the copies
      if (number_of_copies > 1) {
        _copies[i] = _subgrids.size();
      }
      for (uint_fast32_t j = 1; j < number_of_copies; ++j) {
        _subgrids.push_back(new _subgrid_type_(*_subgrids[i]));
        _originals.push_back(i);
      }
    }

    // neighbour setting
    for (int_fast32_t i = 0; i < number_of_unique_subgrids; ++i) {
      const uint_fast8_t level = copy_levels[i];
      const uint_fast32_t number_of_copies = 1 << level;
      // first do the self-reference for each copy (if there are copies)
      for (uint_fast32_t j = 1; j < number_of_copies; ++j) {
        const uint_fast32_t copy = _copies[i] + j - 1;
        _subgrids[copy]->set_neighbour(0, copy);
      }
      // now do the actual neighbours
      for (int_fast32_t j = 1; j < TRAVELDIRECTION_NUMBER; ++j) {
        const uint_fast32_t original_ngb = _subgrids[i]->get_neighbour(j);
        if (original_ngb != NEIGHBOUR_OUTSIDE) {
          const uint_fast8_t ngb_level = copy_levels[original_ngb];
          // check how the neighbour level compares to the subgrid level
          if (ngb_level == level) {
            // same, easy: just make copies mutual neighbours
            // and leave the original grid as is
            for (uint_fast32_t k = 1; k < number_of_copies; ++k) {
              const uint_fast32_t copy = _copies[i] + k - 1;
              const uint_fast32_t ngb_copy = _copies[original_ngb] + k - 1;
              _subgrids[copy]->set_neighbour(j, ngb_copy);
            }
          } else {
            // not the same: there are 2 options
            if (level > ngb_level) {
              // we have less neighbour copies, so some of our copies need to
              // share the same neighbour
              // some of our copies might also need to share the original
              // neighbour
              const uint_fast32_t number_of_ngb_copies = 1
                                                         << (level - ngb_level);
              for (uint_fast32_t k = 1; k < number_of_copies; ++k) {
                const uint_fast32_t copy = _copies[i] + k - 1;
                // this term will round down, which is what we want
                const uint_fast32_t ngb_index = k / number_of_ngb_copies;
                const uint_fast32_t ngb_copy =
                    (ngb_index > 0) ? _copies[original_ngb] + ngb_index - 1
                                    : original_ngb;
                _subgrids[copy]->set_neighbour(j, ngb_copy);
              }
            } else {
              // we have more neighbour copies: pick a subset
              const uint_fast32_t number_of_own_copies = 1
                                                         << (ngb_level - level);
              for (uint_fast32_t k = 1; k < number_of_copies; ++k) {
                const uint_fast32_t copy = _copies[i] + k - 1;
                // the second term will skip some neighbour copies, which is
                // what we want
                const uint_fast32_t ngb_copy =
                    _copies[original_ngb] + (k - 1) * number_of_own_copies;
                _subgrids[copy]->set_neighbour(j, ngb_copy);
              }
            }
          }
        } else {
          // flag this neighbour as NEIGHBOUR_OUTSIDE for all copies
          for (uint_fast32_t k = 1; k < number_of_copies; ++k) {
            const uint_fast32_t copy = _copies[i] + k - 1;
            _subgrids[copy]->set_neighbour(j, NEIGHBOUR_OUTSIDE);
          }
        }
      }
    }
  }

  /**
   * @brief Create copies for the given subgrids according to the given copy
   * level specification.
   *
   * @param copy_levels Desired copy level for each subgrid.
   */
  inline void update_copies(std::vector< uint_fast8_t > &copy_levels) {

    const uint_fast32_t original_number = number_of_original_subgrids();
    for (uint_fast32_t igrid = original_number; igrid < _subgrids.size();
         ++igrid) {
      delete _subgrids[igrid];
    }
    _subgrids.resize(original_number);
    _originals.clear();

    create_copies(copy_levels);
  }

  /**
   * @brief Update the counters of all original subgrids with the contributions
   * from their copies.
   */
  inline void update_original_counters() {
    AtomicValue< size_t > ioriginal(0);
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
    while (ioriginal.value() < _copies.size()) {
      const size_t this_ioriginal = ioriginal.post_increment();
      if (this_ioriginal < _copies.size() &&
          _copies[this_ioriginal] != 0xffffffff) {
        size_t copy_index = _copies[this_ioriginal] - _copies.size();
        while (copy_index < _originals.size() &&
               _originals[copy_index] == this_ioriginal) {
          _subgrids[this_ioriginal]->update_intensities(
              *_subgrids[copy_index + _copies.size()]);
          ++copy_index;
        }
      }
    }
  }

  /**
   * @brief Update the properties of subgrid copies with the changed properties
   * of their original.
   */
  inline void update_copy_properties() {
    AtomicValue< size_t > ioriginal(0);
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
    while (ioriginal.value() < _copies.size()) {
      const size_t this_ioriginal = ioriginal.post_increment();
      if (this_ioriginal < _copies.size() &&
          _copies[this_ioriginal] != 0xffffffff) {
        size_t copy_index = _copies[this_ioriginal] - _copies.size();
        while (copy_index < _originals.size() &&
               _originals[copy_index] == this_ioriginal) {
          _subgrids[copy_index + _copies.size()]->update_neutral_fractions(
              *_subgrids[this_ioriginal]);
          ++copy_index;
        }
      }
    }
  }

  /**
   * @brief iterator to loop over subgrids.
   */
  class iterator {
  private:
    /*! @brief Index of the subgrid the iterator is currently pointing to. */
    size_t _index;

    /*! @brief Pointer to the underlying grid creator (we cannot use a
     *  reference, since then things like it = it would not work). */
    DensitySubGridCreator *_grid_creator;

  public:
    /**
     * @brief Constructor.
     *
     * @param index Index of the cell the iterator is currently pointing to.
     * @param grid_creator DensitySubGridCreator over which we iterate.
     */
    inline iterator(const uint_fast32_t index,
                    DensitySubGridCreator &grid_creator)
        : _index(index), _grid_creator(&grid_creator) {}

    /**
     * @brief Dereference operator.
     *
     * @return Reference to the subgrid the iterator is currently pointing to.
     */
    inline _subgrid_type_ &operator*() {
      return *_grid_creator->_subgrids[_index];
    }

    /**
     * @brief Get an iterator to the first and beyond last copy of the subgrid
     * the iterator is currently pointing to.
     *
     * @return Pair containing the first and last copy of the subgrid the
     * iterator is currently pointing to.
     */
    inline std::pair< iterator, iterator > get_copies() {
      const size_t first_copy = _grid_creator->_copies[_index];
      if (first_copy == 0xffffffff) {
        return std::make_pair(
            iterator(_grid_creator->_subgrids.size(), *_grid_creator),
            iterator(_grid_creator->_subgrids.size(), *_grid_creator));
      }
      size_t last_copy = first_copy;
      while (
          last_copy < _grid_creator->_subgrids.size() &&
          _grid_creator
                  ->_originals[last_copy -
                               _grid_creator->number_of_original_subgrids()] ==
              _index) {
        ++last_copy;
      }
      return std::make_pair(iterator(first_copy, *_grid_creator),
                            iterator(last_copy, *_grid_creator));
    }

    // Iterator functionality

    /**
     * @brief Increment operator.
     *
     * We only implemented the pre-increment version, since the post-increment
     * version creates a new object and is computationally more expensive.
     *
     * @return Reference to the incremented iterator.
     */
    inline iterator &operator++() {
      ++_index;
      return *this;
    }

    /**
     * @brief Increment operator.
     *
     * @param increment Increment to add.
     * @return Reference to the incremented iterator.
     */
    inline iterator &operator+=(const uint_fast32_t increment) {
      _index += increment;
      return *this;
    }

    /**
     * @brief Free addition operator.
     *
     * @param increment Increment to add to the iterator.
     * @return Incremented iterator.
     */
    inline iterator operator+(const uint_fast32_t increment) const {
      iterator it(*this);
      it += increment;
      return it;
    }

    /**
     * @brief Get the index of the subgrid the iterator is currently pointing
     * to.
     *
     * @return Index of the current cell.
     */
    inline size_t get_index() const { return _index; }

    /**
     * @brief Compare iterators.
     *
     * @param it Iterator to compare with.
     * @return True if the iterators point to the same subgrid of the same grid.
     */
    inline bool operator==(iterator it) const {
      return (_grid_creator == it._grid_creator && _index == it._index);
    }

    /**
     * @brief Compare iterators.
     *
     * @param it Iterator to compare with.
     * @return True if the iterators do not point to the same subgrid of the
     * same grid.
     */
    inline bool operator!=(iterator it) const { return !(*this == it); }
  };

  /**
   * @brief Get an iterator to the beginning of the grid.
   *
   * @return Iterator to the first subgrid.
   */
  inline iterator begin() { return iterator(0, *this); }

  /**
   * @brief Get an iterator to the end of the original subgrids.
   *
   * @return Iterator to the end of the original subgrids.
   */
  inline iterator original_end() {
    return iterator(number_of_original_subgrids(), *this);
  }

  /**
   * @brief Get an iterator to the end of all subgrids.
   *
   * @return Iterator to the end of all subgrids.
   */
  inline iterator all_end() { return iterator(_subgrids.size(), *this); }

  /**
   * @brief Get an iterator to the subgrid that contains the given position.
   *
   * @param position Position (in m).
   * @return Iterator to the subgrid that contains that position.
   */
  inline iterator get_subgrid(const CoordinateVector<> position) {
    const int_fast32_t ix =
        std::floor((position.x() - _box.get_anchor().x()) / _subgrid_sides.x());
    const int_fast32_t iy =
        std::floor((position.y() - _box.get_anchor().y()) / _subgrid_sides.y());
    const int_fast32_t iz =
        std::floor((position.z() - _box.get_anchor().z()) / _subgrid_sides.z());
    return iterator(ix * _number_of_subgrids[1] * _number_of_subgrids[2] +
                        iy * _number_of_subgrids[2] + iz,
                    *this);
  }

  /**
   * @brief Get the subgrid with the given index.
   *
   * @param index Index of a subgrid.
   * @return Corresponding subgrid.
   */
  inline iterator get_subgrid(const size_t index) {
    return iterator(index, *this);
  }





std::vector<std::pair<uint_fast32_t, uint_fast32_t>> cells_within_radius(CoordinateVector<double> midpoint, double radius) {
    std::vector<std::pair<uint_fast32_t, uint_fast32_t>> result;

    // Get the grid spacing
    double grid_spacingx = (_box.get_sides()[0]) / (_subgrid_number_of_cells[0] * _number_of_subgrids[0]);
    double grid_spacingy = (_box.get_sides()[1]) / (_subgrid_number_of_cells[1] * _number_of_subgrids[1]);
    double grid_spacingz = (_box.get_sides()[2]) / (_subgrid_number_of_cells[2] * _number_of_subgrids[2]);

    CoordinateVector<double> anchor = _box.get_anchor();

    CoordinateVector<double> grid_spacing(grid_spacingx, grid_spacingy, grid_spacingz);

    // Box dimensions for periodic boundaries
    CoordinateVector<double> box_sides = _box.get_sides();

    // Calculate the minimum and maximum indices for each axis
    CoordinateVector<int> min_idx;
    CoordinateVector<int> max_idx;

    for (uint_fast32_t i = 0; i < 3; i++) {
        min_idx[i] = floor((midpoint[i] - radius - anchor[i]) / grid_spacing[i]);
        max_idx[i] = ceil((midpoint[i] + radius - anchor[i]) / grid_spacing[i]);
        // // Clamp indices for non-periodic boundaries
        // if (!_periodicity[i]) {
        //     if (min_idx[i] < 0) {
        //       std::cout << "SNe leaking from axis " << i << " in the negative direction" << std::endl;
        //         min_idx[i] = 0;
        //     }
        //     if (max_idx[i] >= _subgrid_number_of_cells[i] * _number_of_subgrids[i]) {
        //       std::cout << "SNe leaking from axis " << i << " in the positive direction" << std::endl;
        //         max_idx[i] = _subgrid_number_of_cells[i] * _number_of_subgrids[i] - 1;
        //     }
        // }
    }


    // Loop through the subgrids within the bounds, including periodic wrapping
    for (int x_idx = min_idx.x(); x_idx <= max_idx.x(); ++x_idx) {
        for (int y_idx = min_idx.y(); y_idx <= max_idx.y(); ++y_idx) {
            for (int z_idx = min_idx.z(); z_idx <= max_idx.z(); ++z_idx) {
                // Wrap indices for periodic dimensions
                int wrapped_x_idx = x_idx;
                int wrapped_y_idx = y_idx;
                int wrapped_z_idx = z_idx;

                double x_offset = 0.0, y_offset = 0.0, z_offset = 0.0;

                if (_periodicity[0]) {
                    if (x_idx < 0) {
                        wrapped_x_idx = x_idx + _subgrid_number_of_cells[0] * _number_of_subgrids[0];
                        x_offset = -box_sides[0];
                    } else if (x_idx >= _subgrid_number_of_cells[0] * _number_of_subgrids[0]) {
                        wrapped_x_idx = x_idx - _subgrid_number_of_cells[0] * _number_of_subgrids[0];
                        x_offset = box_sides[0];
                    }
                } else {
                  if (x_idx < 0) {
                    continue;
                  } else if (x_idx >= _subgrid_number_of_cells[0] * _number_of_subgrids[0]) {
                    continue;
                  }
                }
                if (_periodicity[1]) {
                    if (y_idx < 0) {
                        wrapped_y_idx = y_idx + _subgrid_number_of_cells[1] * _number_of_subgrids[1];
                        y_offset = -box_sides[1];
                    } else if (y_idx >= _subgrid_number_of_cells[1] * _number_of_subgrids[1]) {
                        wrapped_y_idx = y_idx - _subgrid_number_of_cells[1] * _number_of_subgrids[1];
                        y_offset = box_sides[1];
                    }
                } else {
                  if (y_idx < 0) {
                    continue;
                  } else if (y_idx >= _subgrid_number_of_cells[1] * _number_of_subgrids[1]) {
                    continue;
                  }
                }
                if (_periodicity[2]) {
                    if (z_idx < 0) {
                        wrapped_z_idx = z_idx + _subgrid_number_of_cells[2] * _number_of_subgrids[2];
                        z_offset = -box_sides[2];
                    } else if (z_idx >= _subgrid_number_of_cells[2] * _number_of_subgrids[2]) {
                        wrapped_z_idx = z_idx - _subgrid_number_of_cells[2] * _number_of_subgrids[2];
                        z_offset = box_sides[2];
                    }
                } else {
                  if (z_idx < 0) {
                    continue;
                  } else if (z_idx >= _subgrid_number_of_cells[2] * _number_of_subgrids[2]) {
                    continue;
                  }
                }

                // Wrapped position for subgrid lookup
                CoordinateVector<double> wrapped_midpoint(wrapped_x_idx * grid_spacing[0],
                                                          wrapped_y_idx * grid_spacing[1],
                                                          wrapped_z_idx * grid_spacing[2]);
                wrapped_midpoint += anchor;
                wrapped_midpoint[0] += grid_spacing[0] / 2.;
                wrapped_midpoint[1] += grid_spacing[1] / 2.;
                wrapped_midpoint[2] += grid_spacing[2] / 2.;

                // Unwrapped position for distance calculation
                CoordinateVector<double> unwrapped_midpoint = wrapped_midpoint;
                unwrapped_midpoint[0] += x_offset;
                unwrapped_midpoint[1] += y_offset;
                unwrapped_midpoint[2] += z_offset;

                // Check if the cell is within the specified radius
                double distance = (midpoint - unwrapped_midpoint).norm();
                if (distance <= radius) {
                    uint_fast32_t subgrid_idx = this->get_subgrid(wrapped_midpoint).get_index();
                    HydroDensitySubGrid& subgrid = *this->get_subgrid(subgrid_idx);
                    uint_fast32_t cell_idx = subgrid.get_cell(wrapped_midpoint).get_index();

                    // Add the pair of indices to the result vector
                    result.emplace_back(subgrid_idx, cell_idx);
                }
            }
        }
    }

    return result;
}


  // std::vector<std::pair<uint_fast32_t, uint_fast32_t>> cells_within_radius(CoordinateVector<double> midpoint, double radius) {
  //   std::vector<std::pair<uint_fast32_t, uint_fast32_t>> result;

  //   // Get the grid spacing
  //   double grid_spacingx = (_box.get_sides()[0])/(_subgrid_number_of_cells[0]*_number_of_subgrids[0]);
  //   double grid_spacingy = (_box.get_sides()[1])/(_subgrid_number_of_cells[1]*_number_of_subgrids[1]);
  //   double grid_spacingz = (_box.get_sides()[2])/(_subgrid_number_of_cells[2]*_number_of_subgrids[2]);

  //   CoordinateVector<double> anchor = _box.get_anchor();

  //   CoordinateVector<double> grid_spacing(grid_spacingx,grid_spacingy,grid_spacingz);




  //   // Calculate the minimum and maximum indices for each axis
  //   CoordinateVector<int> min_idx;
  //   CoordinateVector<int> max_idx;

  //   for (uint_fast32_t i = 0; i < 3; i++) {
  //     min_idx[i] = floor((midpoint[i] - radius - anchor[i]) / grid_spacing[i]);
  //     if (min_idx[i] < 0) {
  //       min_idx[i] = 0;
  //     }
  //     max_idx[i] = ceil((midpoint[i] + radius - anchor[i]) / grid_spacing[i]);
  //     if (max_idx[i] > _subgrid_number_of_cells[i]*_number_of_subgrids[i] -1) {
  //       max_idx[i] = _subgrid_number_of_cells[i]*_number_of_subgrids[i] -1;
  //     }

  //   }

  //   // Loop through the subgrids within the bounds
  //   for (int x_idx = min_idx.x(); x_idx <= max_idx.x(); ++x_idx) {
  //     for (int y_idx = min_idx.y(); y_idx <= max_idx.y(); ++y_idx) {
  //       for (int z_idx = min_idx.z(); z_idx <= max_idx.z(); ++z_idx) {
  //         CoordinateVector<double> cell_midpoint(x_idx * grid_spacing[0], y_idx * grid_spacing[1], z_idx * grid_spacing[2]);

  //         cell_midpoint += anchor;
  //         cell_midpoint[0] += grid_spacing[0]/2.;
  //         cell_midpoint[1] += grid_spacing[1]/2.;
  //         cell_midpoint[2] += grid_spacing[2]/2.;



  //         // Check if the cell is within the specified radius
  //         double distance = (midpoint - cell_midpoint).norm();
  //         if (distance <= radius) {
  //           uint_fast32_t subgrid_idx = this->get_subgrid(cell_midpoint).get_index();
  //           HydroDensitySubGrid& subgrid = *this->get_subgrid(subgrid_idx);
  //           uint_fast32_t cell_idx = subgrid.get_cell(cell_midpoint).get_index();

  //           // Add the pair of indices to the result vector
  //           result.emplace_back(subgrid_idx, cell_idx);
  //         }
  //       }
  //     }
  //   }

  //   return result;
  // }


  /**
   * @brief Dump the subgrids to the given restart file.
   *
   * @param restart_writer RestartWriter to write to.
   */
  inline void write_restart_file(RestartWriter &restart_writer) const {

    // const members
    _box.write_restart_file(restart_writer);
    _subgrid_sides.write_restart_file(restart_writer);
    _number_of_subgrids.write_restart_file(restart_writer);
    _subgrid_number_of_cells.write_restart_file(restart_writer);
    _periodicity.write_restart_file(restart_writer);

    const size_t number_of_subgrids = _subgrids.size();
    restart_writer.write(number_of_subgrids);
    for (size_t i = 0; i < number_of_subgrids; ++i) {
      _subgrids[i]->write_restart_file(restart_writer);
    }
    const size_t number_of_copies = _originals.size();
    restart_writer.write(number_of_copies);
    for (size_t i = 0; i < number_of_copies; ++i) {
      restart_writer.write(_originals[i]);
    }
    const size_t number_of_originals = _copies.size();
    restart_writer.write(number_of_originals);
    for (size_t i = 0; i < number_of_originals; ++i) {
      restart_writer.write(_copies[i]);
    }
  }

  /**
   * @brief Restart constructor.
   *
   * @param restart_reader Restart file to read from.
   */
  inline DensitySubGridCreator(RestartReader &restart_reader)
      : _box(restart_reader), _subgrid_sides(restart_reader),
        _number_of_subgrids(restart_reader),
        _subgrid_number_of_cells(restart_reader), _periodicity(restart_reader) {

    const size_t number_of_subgrids = restart_reader.read< size_t >();
    _subgrids.resize(number_of_subgrids, nullptr);
    for (size_t i = 0; i < number_of_subgrids; ++i) {
      _subgrids[i] = new _subgrid_type_(restart_reader);
    }
    const size_t number_of_copies = restart_reader.read< size_t >();
    _originals.resize(number_of_copies, 0);
    for (size_t i = 0; i < number_of_copies; ++i) {
      _originals[i] = restart_reader.read< size_t >();
    }
    const size_t number_of_originals = restart_reader.read< size_t >();
    _copies.resize(number_of_originals, 0);
    for (size_t i = 0; i < number_of_originals; ++i) {
      _copies[i] = restart_reader.read< size_t >();
    }
  }
};

#endif // DENSITYSUBGRIDCREATOR_HPP
