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
 * @file ArepoSnapshotDensityFunction.cpp
 *
 * @brief DensityFunction that reads a density field from an Arepo HDF5 snapshot
 * file: implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "ArepoSnapshotDensityFunction.hpp"
#include "HDF5Tools.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "UnitConverter.hpp"
#include <cfloat>
#include <fstream>
#include <string>



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
ArepoSnapshotDensityFunction::ArepoSnapshotDensityFunction(
    std::string name, bool fallback_periodic, double fallback_unit_length_in_SI,
    double fallback_unit_mass_in_SI, double fallback_unit_temperature_in_SI,
    bool use_neutral_fraction, double fallback_temperature,
    bool comoving_integration, double hubble_parameter,
    double dust_to_gas, double fraction_silicates, bool read_temp_nfrac,
    Log *log)
    : _log(log) {

  _dust_to_gas = dust_to_gas;
  _fraction_silicates = fraction_silicates;

  _read_temp_nfrac = read_temp_nfrac;
  // turn off default HDF5 error handling: we catch errors ourselves
  HDF5Tools::initialize();

  // open the file for reading
  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(name, HDF5Tools::HDF5FILEMODE_READ);
  // open the Parameters group
  bool periodic = fallback_periodic;
  double unit_length_in_SI = fallback_unit_length_in_SI;
  double unit_mass_in_SI = fallback_unit_mass_in_SI;

  double unit_vel_in_SI = 1;
  if (HDF5Tools::group_exists(file, "/Parameters")) {
    HDF5Tools::HDF5Group runtimepars =
        HDF5Tools::open_group(file, "/Parameters");
    // read the PeriodicBoundariesOn flag
    periodic = HDF5Tools::read_attribute< int32_t >(
                   runtimepars, "PeriodicBoundariesOn") != 0;
    const double unit_length_in_cgs =
                   HDF5Tools::read_attribute< double >(runtimepars, "UnitLength_in_cm");
    const double unit_mass_in_cgs =
                  HDF5Tools::read_attribute< double >(runtimepars, "UnitMass_in_g");
    const double unit_velocity_in_cgs =
                  HDF5Tools::read_attribute< double >(runtimepars, "UnitVelocity_in_cm_per_s");
    unit_length_in_SI =
                  UnitConverter::to_SI< QUANTITY_LENGTH >(unit_length_in_cgs, "cm");
    unit_mass_in_SI =
                  UnitConverter::to_SI< QUANTITY_MASS >(unit_mass_in_cgs, "g");
    unit_vel_in_SI = unit_velocity_in_cgs/100.0;
    if (fallback_unit_temperature_in_SI == 0.) {
      if (_log) {
        _log->write_warning(
            "No fallback temperature unit found in parameter file, using K!");
      }

    }
    // close the group
    HDF5Tools::close_group(runtimepars);
  } else {
    _log->write_warning("No Parameters found!");
    if (periodic) {
      _log->write_warning("Assuming a periodic box.");
    } else {
      _log->write_warning("Assuming a non-periodic box.");
    }
  }
  if (periodic) {
    // open the Header group
    HDF5Tools::HDF5Group header = HDF5Tools::open_group(file, "/Header");
    // Read the box size
    CoordinateVector<> sides =
        HDF5Tools::read_attribute< CoordinateVector<> >(header, "BoxSize");
    // in this case, the anchor is just (0., 0., 0.)
    CoordinateVector<> anchor;
    _box = Box<>(anchor, sides);
    // close the Header group
    HDF5Tools::close_group(header);
  }


  if (comoving_integration) {
    // code values are in comoving units
    unit_length_in_SI /= hubble_parameter;
    unit_mass_in_SI /= hubble_parameter;
  }

  const double unit_length_in_SI_squared =
      unit_length_in_SI * unit_length_in_SI;
  const double unit_density_in_SI =
      unit_mass_in_SI / unit_length_in_SI / unit_length_in_SI_squared;

  // open the group containing the cell data
  HDF5Tools::HDF5Group voronoicells = HDF5Tools::open_group(file, "/PartType0");
  // read the positions, masses
  _positions = HDF5Tools::read_dataset< CoordinateVector<> >(voronoicells,
                                                             "Coordinates");

  _velocities = HDF5Tools::read_dataset< CoordinateVector<> >(voronoicells,
                                                            "Velocities");
  _densities = HDF5Tools::read_dataset< double >(voronoicells, "Density");
  _masses = HDF5Tools::read_dataset< double >(voronoicells, "Masses");
  if (HDF5Tools::group_exists(voronoicells, "InternalEnergy")) {
    _internal_energy =
        HDF5Tools::read_dataset< double >(voronoicells, "InternalEnergy");
  } else {
    if (_log) {
      _log->write_warning("No internal energy block found to calculate temperature, using fallback initial "
                          "temperature value.");
    }
    if (fallback_temperature == 0.) {
      fallback_temperature = 8000.;
      if (_log) {
        _log->write_warning(
            "No fallback initial temperature provided either, using 8000. K.");
      }
    }
    _temperatures.resize(_densities.size(), fallback_temperature);

  }
  // close the group



  _chemabundances = HDF5Tools::read_dataset< CoordinateVector<> >(voronoicells,
                                                             "ChemicalAbundances");


  HDF5Tools::close_group(voronoicells);

  // close the file
  HDF5Tools::close_file(file);

  // unit conversion + treebox data collection
  CoordinateVector<> minpos(DBL_MAX);
  CoordinateVector<> maxpos(-DBL_MAX);
  for (size_t i = 0; i < _positions.size(); ++i) {
    _positions[i][0] *= unit_length_in_SI;
    _positions[i][1] *= unit_length_in_SI;
    _positions[i][2] *= unit_length_in_SI;
    _masses[i] *= unit_mass_in_SI;
    _densities[i] *= unit_density_in_SI;
    _internal_energy[i] *= unit_vel_in_SI*unit_vel_in_SI;
    minpos = CoordinateVector<>::min(minpos, _positions[i]);
    maxpos = CoordinateVector<>::max(maxpos, _positions[i]);
  }
  if (periodic) {
    _box.get_anchor()[0] *= unit_length_in_SI;
    _box.get_anchor()[1] *= unit_length_in_SI;
    _box.get_anchor()[2] *= unit_length_in_SI;
    _box.get_sides()[0] *= unit_length_in_SI;
    _box.get_sides()[1] *= unit_length_in_SI;
    _box.get_sides()[2] *= unit_length_in_SI;
  }


  if (_log) {
    _log->write_status("Successfully read densities from file \"", name, "\".");
  }


  Box<> box(_box);
  if (!periodic) {
    // set box to particle extents + small margin
    CoordinateVector<> sides = maxpos - minpos;
    const CoordinateVector<> anchor = minpos - 0.005 * sides;
    sides *= 1.01;
    box = Box<>(anchor, sides);
  }
  if (_log) {
    std::string pstring;
    if (periodic) {
      pstring = "periodic ";
    }
    _log->write_status("Creating octree in ", pstring, "box with anchor [",
                       box.get_anchor().x(), " m, ", box.get_anchor().y(),
                       " m, ", box.get_anchor().z(), " m] and sides [",
                       box.get_sides().x(), " m, ", box.get_sides().y(), " m, ",
                       box.get_sides().z(), " m]...");
  }

  _octree = new Octree(_positions, box, periodic);

  if (_log) {
    _log->write_status("Done creating octree.");
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
ArepoSnapshotDensityFunction::ArepoSnapshotDensityFunction(
    ParameterFile &params, Log *log)
    : ArepoSnapshotDensityFunction(
          params.get_filename("DensityFunction:filename"),
          params.get_value< bool >("DensityFunction:fallback periodic flag",
                                   false),
          params.get_physical_value< QUANTITY_LENGTH >(
              "DensityFunction:fallback unit length", "0. m"),
          params.get_physical_value< QUANTITY_MASS >(
              "DensityFunction:fallback unit mass", "0. kg"),
          params.get_physical_value< QUANTITY_TEMPERATURE >(
              "DensityFunction:fallback unit temperature", "0. K"),
          params.get_value< bool >("DensityFunction:use neutral fraction",
                                   false),
          params.get_physical_value< QUANTITY_TEMPERATURE >(
              "DensityFunction:fallback initial temperature", "0. K"),
          params.get_value< bool >("DensityFunction:comoving integration flag",
                                   false),
          params.get_value< double >("DensityFunction:hubble parameter", 0.7),
          params.get_value<double>("DensityFunction:dust to gas",0.0),
          params.get_value<double>("DensityFunction:fraction silicates",0.6),
          params.get_value<bool>("DensityFunction:read temperature",true),
          log) {}

/**
 * @brief Destructor.
 *
 * Deletes the internal Octree.
 */
ArepoSnapshotDensityFunction::~ArepoSnapshotDensityFunction() {
  delete _octree;
}

/**
 * @brief Function that gives the density for a given cell.
 *
 * @param cell Geometrical information about the cell.
 * @return Initial physical field values for that cell.
 */
DensityValues ArepoSnapshotDensityFunction::operator()(const Cell &cell) {

  DensityValues values;

  const CoordinateVector<> position = cell.get_cell_midpoint();
  double density = 0.;
  double temperature = 0.;
  double neutral_fraction = -1.;
  std::vector<double> chemistry(5,0.0);

  const uint_fast32_t ngb = _octree->get_closest_ngb(position);
  CoordinateVector<> chems = _chemabundances[ngb];
  density = _densities[ngb];
  // nHtot = nHI + nHp + 2*nH2
  // nTOT = nHI + nH2 + nHp + ne + nHe
  double nHtot = density/((1. + 4.0 * 0.1) * 1.6737236e-27);
  double ntot = nHtot*(1.0 + chems[1] - chems[0] + 0.1);

  double nCO = chems[2]*nHtot;
  double nHp = chems[1]*nHtot;
  double nH2 = chems[0]*nHtot;
  double nHI = (1.0 - chems[1] - 2.0*chems[0])*nHtot;
  double mu = density/(ntot*1.6726231e-27);

  chemistry[0] = nHI;
  chemistry[1] = nH2;
  chemistry[2] = nHp;
  chemistry[3] = nCO;
  chemistry[4] = mu;



  //values.set_chemistry(chemistry);

  values.set_velocity(_velocities[ngb]);





  temperature = (2.0/3.0)*_internal_energy[ngb]*mu*1.6726231e-27/1.3806485e-23;


  neutral_fraction = chems[1];
  if (!_read_temp_nfrac) {
    temperature = 10000;
    neutral_fraction = 1.e-4;
  }



  values.set_number_density(nHtot);
  values.set_temperature(temperature);
  values.set_dust_gas_ratio(_dust_to_gas);
  values.set_fraction_silicates(_fraction_silicates);


  if (neutral_fraction >= 0.) {
    values.set_ionic_fraction(ION_H_n, neutral_fraction);
  } else {
    values.set_ionic_fraction(ION_H_n, 1.e-6);
  }
#ifdef HAS_HELIUM
  values.set_ionic_fraction(ION_He_n, 1.e-6);
#endif
  return values;
}

/**
 * @brief Get the total number of hydrogen atoms in the snapshot.
 *
 * @return Sum of the hydrogen number of all Voronoi cells in the snapshot.
 */
double ArepoSnapshotDensityFunction::get_total_hydrogen_number() const {
  double mtot = 0.;
  for (size_t i = 0; i < _masses.size(); ++i) {
    mtot += _masses[i];
  }
  return mtot / 1.6737236e-27;
}
