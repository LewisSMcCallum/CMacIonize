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
 * @file HDF5DensityFunction.cpp
 *
 * @brief DensityFunction that reads a density field from an Arepo HDF5 snapshot
 * file: implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "HDF5DensityFunction.hpp"
#include "HDF5Tools.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "UnitConverter.hpp"
#include <cfloat>
#include <fstream>
#include <string>
#include <iostream>



/**
 * @brief Constructor.
 *
 * @param name Name of the snapshot file to read.
 * @param fallback_periodic Periodicity flag used in case RuntimePars group is
 * not found in the snapshot file.
 * @param fallback_unit_length_in_SI Length unit to use if the Units group is
 * not found in the snapshot file.
 * @param fallback_unit_mass_in_SI Mass unit to use if the Units group is not
 * found in the snapshot file.
 * @param fallback_unit_temperature_in_SI Temperature unit to use if the Units
 * group is not found in the snapshot file.
 * @param use_neutral_fraction Whether or not to use the neutral fraction from
 * the snapshot file as initial guess for the neutral fraction (if it is
 * present in the snapshot).
 * @param fallback_temperature Initial temperature to use if no temperature
 * block was found in the snapshot file.
 * @param comoving_integration Comoving integration flag indicating whether
 * comoving integration was used in the snapshot.
 * @param hubble_parameter Hubble parameter used to convert from comoving to
 * physical coordinates. This is a dimensionless parameter, defined as the
 * actual assumed Hubble constant divided by 100 km/s/Mpc.
 * @param log Log to write logging information to.
 */
HDF5DensityFunction::HDF5DensityFunction(
    std::string name, CoordinateVector< uint_fast32_t > ncell,
                const double dust_gas_ratio, const double fraction_silicates,
                const Box<> &simulation_box, Log *log)
    : _ncell(ncell), _dust_gas_ratio(dust_gas_ratio),
         _fraction_silicates(fraction_silicates),_log(log)  {


   _box = simulation_box;


   _cell_vol = (_box.get_sides()[0]/_ncell[0])*(_box.get_sides()[1]/_ncell[1])*(_box.get_sides()[2]/_ncell[2]);


  // turn off default HDF5 error handling: we catch errors ourselves
  HDF5Tools::initialize();

  // open the file for reading
  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(name, HDF5Tools::HDF5FILEMODE_READ);
  // open the Parameters group



  // open the group containing the cell data
  HDF5Tools::HDF5Group maingroup = HDF5Tools::open_group(file, "/PartType0");
  // read the positions, masses
  _densities = HDF5Tools::read_dataset< double >(maingroup,
                                                             "NumberDensity");







  HDF5Tools::close_group(maingroup);

  // close the file
  HDF5Tools::close_file(file);



  if (_log) {
    _log->write_status("Successfully read densities from file \"", name, "\".");
  }

}




/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - filename: Name of the snapshot file (required)
 *  - fallback periodic flag: Periodicity flag to use when no periodicity flag
 *    can be found in the snapshot file (default: false)
 *  - fallback unit length: Length unit to use if no units can be found in the
 *    snapshot file (default: 1. m, with warning)
 *  - fallback unit mass: Mass unit to use if no units can be found in the
 *    snapshot file (default: 1. kg, with warning)
 *  - fallback unit temperature: Temperature unit to use if no units can be
 *    found in the snapshot file (default: 1. K, with warning)
 *  - use neutral fraction: Use initial neutral fractions from the snapshot file
 *    (if present; default: false)?
 *  - fallback initial temperature: Initial temperature to use if no temperature
 *    values can be found in the snapshot file (default: 8000. K, with warning)
 *  - comoving integration flag: Was comoving integration used in the original
 *    simulation (default: false)?
 *  - hubble parameter: Reduced Hubble parameter used for the original
 *    simulation (default: 0.7)
 *
 * @param params ParameterFile to read.
 * @param log Log to write logging information to.
 */
 HDF5DensityFunction::HDF5DensityFunction(ParameterFile &params,
                                                     Log *log)
      : HDF5DensityFunction(
            params.get_filename("DensityFunction:filename"),
            params.get_value< CoordinateVector< uint_fast32_t > >(
                "DensityGrid:number of cells",
                CoordinateVector< uint_fast32_t >(64)),
            params.get_value<double>("DensityFunction:dust to gas",0.01),
            params.get_value<double>("DensityFunction:fraction silicates",1.0),
            Box<>(params.get_physical_vector< QUANTITY_LENGTH >(
                      "SimulationBox:anchor", "[-5. pc, -5. pc, -5. pc]"),
                  params.get_physical_vector< QUANTITY_LENGTH >(
                      "SimulationBox:sides", "[10. pc, 10. pc, 10. pc]")),
            log) {}

/**
 * @brief Destructor.
 *
 * Deletes the internal Octree.
 */
HDF5DensityFunction::~HDF5DensityFunction() {

}

/**
 * @brief Function that gives the density for a given cell.
 *
 * @param cell Geometrical information about the cell.
 * @return Initial physical field values for that cell.
 */
DensityValues HDF5DensityFunction::operator()(const Cell &cell) {

  DensityValues values;

  const CoordinateVector<> position = cell.get_cell_midpoint();
  uint_fast32_t ix, iy, iz;
   ix = (position.x() - _box.get_anchor().x()) / _box.get_sides().x() *
        _ncell.x();
   iy = (position.y() - _box.get_anchor().y()) / _box.get_sides().y() *
        _ncell.y();
   iz = (position.z() - _box.get_anchor().z()) / _box.get_sides().z() *
        _ncell.z();




   int ind1d = (ix*_ncell.z()*_ncell.y()) + (iy*_ncell.z()) + iz;


   values.set_number_density(_densities[ind1d]);

  

   values.set_temperature(8000);
   values.set_ionic_fraction(ION_H_n,1e-6);


  return values;
}

/**
 * @brief Get the total number of hydrogen atoms in the snapshot.
 *
 * @return Sum of the hydrogen number of all Voronoi cells in the snapshot.
 */
double HDF5DensityFunction::get_total_hydrogen_number() const {
  double mtot = 0.;
  for (size_t i = 0; i < _densities.size(); ++i) {
    mtot += _densities[i]*_cell_vol;
  }
  return mtot / 1.6737236e-27;
}
