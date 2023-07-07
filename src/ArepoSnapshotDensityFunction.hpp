/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file GadgetSnapshotDensityFunction.hpp
 *
 * @brief DensityFunction that reads a density field from a Gadget HDF5 snapshot
 * file: header.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef AREPOSNAPSHOTDENSITYFUNCTION_HPP
#define AREPOSNAPSHOTDENSITYFUNCTION_HPP

#include "Box.hpp"
#include "DensityFunction.hpp"
#include "Octree.hpp"
#include <string>
#include <vector>

class Log;
class ParameterFile;

/**
 * @brief DensityFunction that reads a density field from an Arepo snapshot.
 */
class ArepoSnapshotDensityFunction : public DensityFunction {
private:
  /*! @brief Simulation box, only initialized if the box is periodic (if the box
   *  is not periodic, the components of the Box will all be zero). */
  Box<> _box;

  /*! @brief Positions of the Voronoi mesh generating points in the snapshot (in m). */
  std::vector< CoordinateVector<> > _positions;

  std::vector< CoordinateVector<> > _velocities;


  std::vector< CoordinateVector<> > _chemabundances;


  /*! @brief Densities of the Voronoi cells in the snapshot (in kg m^-3).*/

  std::vector< double > _densities;

  /*! @brief Temperatures of the Voronoi cells in the snapshot (in K). */
  std::vector< double > _temperatures;

  /*! @brief Neutral fractions of the Voronoi cells in the snapshot (if
   *  present). */
  std::vector< double > _neutral_fractions;

  /* @brief Masses of the Voronoi cells in the snapshot */

  std::vector< double > _masses;

  std::vector< double > _internal_energy;

  double _dust_to_gas;
  double _fraction_silicates;

  bool _read_temp_nfrac;



  /*! @brief Octree used to speed up neighbour searching. */
  Octree *_octree;

  /*! @brief Log to write logging info to. */
  Log *_log;

public:
  ArepoSnapshotDensityFunction(std::string name,
                                bool fallback_periodic = false,
                                double fallback_unit_length_in_SI = 1.,
                                double fallback_unit_mass_in_SI = 1.,
                                double fallback_unit_temperature_in_SI = 1.,
                                bool use_neutral_fraction = false,
                                double fallback_temperature = 8000.,
                                bool comoving_integration = false,
                                double hubble_parameter = 0.7,
                                double dust_to_gas = 0.0,
                                double fraction_silicates=0.6,
                                bool read_temp_nfrac=true,
                                Log *log = nullptr);

  ArepoSnapshotDensityFunction(ParameterFile &params, Log *log = nullptr);

  virtual ~ArepoSnapshotDensityFunction();

  virtual DensityValues operator()(const Cell &cell);

  double get_total_hydrogen_number() const;
};

#endif // AREPOSNAPSHOTDENSITYFUNCTION_HPP
