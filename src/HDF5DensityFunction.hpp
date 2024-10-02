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
#ifndef HDF5DENSITYFUNCTION_HPP
#define HDF5DENSITYFUNCTION_HPP

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
class HDF5DensityFunction : public DensityFunction {
private:
  /*! @brief Simulation box, only initialized if the box is periodic (if the box
   *  is not periodic, the components of the Box will all be zero). */


   /*! @brief Densities of the Voronoi cells in the snapshot (in kg m^-3).*/

   std::vector<  double >  _densities;
   std::vector<  double >  _temperatures;
  

   double _cell_vol;
  CoordinateVector< uint_fast32_t > _ncell;

  bool _read_temps;

  double _dust_gas_ratio;

  double _fraction_silicates;

  Box<> _box;


  /*! @brief Log to write logging info to. */
  Log *_log;

  /*! @brief Positions of the Voronoi mesh generating points in the snapshot (in m). */

  //std::vector< CoordinateVector<> > _positions;







public:
  HDF5DensityFunction(std::string name,
                           CoordinateVector< uint_fast32_t > ncell,
                           bool read_temps,
                             const double dust_gas_ratio,
                             const double fraction_silicates,
                                const Box<> &simulation_box,
                                Log *log = nullptr);

  HDF5DensityFunction(ParameterFile &params, Log *log = nullptr);

  virtual ~HDF5DensityFunction();

  virtual DensityValues operator()(const Cell &cell);

  double get_total_hydrogen_number() const;
};

#endif // HDF5DENSITYFUNCTION_HPP
