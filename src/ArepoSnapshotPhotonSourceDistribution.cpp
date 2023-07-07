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
 * @file ArepoSnapshotPhotonSourceDistribution.cpp
 *
 * @brief ArepoSnapshotPhotonSourceDistribution implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "ArepoSnapshotPhotonSourceDistribution.hpp"
#include "WMBasicPhotonSourceSpectrum.hpp"
#include "HDF5Tools.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "UVLuminosityFunctionFactory.hpp"
#include "UnitConverter.hpp"

/**
 * @brief Constructor.
 *
 * Reads in the sources from the file and stores them in internal arrays.
 *
 * @param filename Name of the snapshot file to read.
 * @param formation_time_name Name of the formation time data set in the
 * snapshot file.
 * @param box Dimensions of the simulation box (in m).
 * @param luminosity_function UVLuminosityFunction to use to assign UV
 * luminosities to star particles.
 * @param fallback_unit_length_in_SI Length unit to use if the units group is
 * not found in the snapshot file.
 * @param fallback_unit_time_in_SI Time unit to use if the units group is not
 * found in the snapshot file.
 * @param fallback_unit_mass_in_SI Mass unit to use if the units group is not
 * found in the snapshot file.
 * @param cutoff_age Upper age limit for stellar populations that emit UV
 * radiation (in s).
 * @param use_gas Do the gas particles contain stars?
 * @param SFR_unit Unit used for star formation rate values in the snapshot (if
 * set to 0., mass_unit/time_unit is assumed).
 * @param comoving_integration Comoving integration flag indicating whether
 * comoving integration was used in the snapshot.
 * @param hubble_parameter  Hubble parameter used to convert from comoving to
 * physical coordinates. This is a dimensionless parameter, defined as the
 * actual assumed Hubble constant divided by 100 km/s/Mpc.
 * @param log Log to write logging information to.



 */


 int pop_from_temp(double teff){
    if (teff < 35500){
          return 0;
    } else if (teff >= 35000 and teff < 38500){
         return 1;
    } else if(teff >= 38500 and teff < 41500){
          return 2;
    }  else if (teff >= 41500 and teff < 44500){
          return 3;
    }  else if (teff >= 44500 and teff < 47500){
          return 4;
    }  else if (teff >= 47500){
          return 5;
    } else {
      std::cout << "SOMETHING GONE WRONG IN SETTING UP SPECTRA TYPE" << std::endl;
      return -1;
    }
  }

ArepoSnapshotPhotonSourceDistribution::ArepoSnapshotPhotonSourceDistribution(
    std::string filename, std::string formation_time_name, const Box<> box,
    UVLuminosityFunction *luminosity_function,
    double fallback_unit_length_in_SI, double fallback_unit_time_in_SI,
    double fallback_unit_mass_in_SI, double cutoff_age, bool use_gas,
    double SFR_unit, bool comoving_integration, double hubble_parameter,
    Log *log)
    : _log(log) {





  _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(34000,60,log));
  _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(37000,60,log));
  _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(40000,60,log));
  _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(43000,60,log));
  _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(46000,60,log));
  _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(49000,60,log));

  // turn off default HDF5 error handling: we catch errors ourselves
  HDF5Tools::initialize();

  // open the file for reading
  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_READ);


  // units
  double unit_length_in_SI = fallback_unit_length_in_SI;
  double unit_time_in_SI = fallback_unit_time_in_SI;
  double unit_mass_in_SI = fallback_unit_mass_in_SI;
  if (HDF5Tools::group_exists(file, "/Parameters")) {
    HDF5Tools::HDF5Group units = HDF5Tools::open_group(file, "/Parameters");
    const double unit_length_in_cgs =
        HDF5Tools::read_attribute< double >(units, "UnitLength_in_cm");
    const double unit_vel_in_cgs =
        HDF5Tools::read_attribute< double >(units, "UnitVelocity_in_cm_per_s");
    const double unit_mass_in_cgs =
        HDF5Tools::read_attribute< double >(units, "UnitMass_in_g");
    const double unit_time_in_cgs =
        unit_length_in_cgs/unit_vel_in_cgs;
    unit_length_in_SI =
        UnitConverter::to_SI< QUANTITY_LENGTH >(unit_length_in_cgs, "cm");
    // seconds are seconds
    unit_time_in_SI = unit_time_in_cgs;
    unit_mass_in_SI =
        UnitConverter::to_SI< QUANTITY_MASS >(unit_mass_in_cgs, "g");
    HDF5Tools::close_group(units);
  } else {
    if (_log) {
      _log->write_warning("No Units group found!");
      _log->write_warning("Using fallback units.");
    }

    if (fallback_unit_length_in_SI == 0.) {
      if (_log) {
        _log->write_warning(
            "No fallback length unit found in parameter file, using m!");
      }
      unit_length_in_SI = 1.;
    }

    if (fallback_unit_time_in_SI == 0.) {
      if (_log) {
        _log->write_warning(
            "No fallback time unit found in parameter file, using s!");
      }
      unit_time_in_SI = 1.;
    }

    if (fallback_unit_mass_in_SI == 0.) {
      if (_log) {
        _log->write_warning(
            "No fallback mass unit found in parameter file, using kg!");
      }
      unit_mass_in_SI = 1.;
    }
  }

  if (comoving_integration) {
    // code values are in comoving units
    unit_length_in_SI /= hubble_parameter;
    unit_mass_in_SI /= hubble_parameter;
    unit_time_in_SI /= hubble_parameter;
  }

  _total_luminosity = 0.;

  // open the group containing the star particle data
  HDF5Tools::HDF5Group clusterparticles =
      HDF5Tools::open_group(file, "/SinkSources");
  // read the positions
  std::vector< CoordinateVector<> > positions =
      HDF5Tools::read_dataset< CoordinateVector<> >(clusterparticles,
                                                      "Coordinates");
  // read the luminosites
  std::vector< double > luminosity =
      HDF5Tools::read_dataset< double >(clusterparticles, "Luminosity");


  std::vector<double> teff =
      HDF5Tools::read_dataset<double >(clusterparticles,"EffectiveTemp");
  // close the group
  HDF5Tools::close_group(clusterparticles);
  // close the file
  HDF5Tools::close_file(file);

  for (size_t i = 0; i < luminosity.size(); ++i) {
    const CoordinateVector<> position = positions[i] * unit_length_in_SI;
    if (box.inside(position)) {
      const double UV_luminosity = luminosity[i];
      if (UV_luminosity > 0.) {
        _positions.push_back(position);
        _luminosities.push_back(UV_luminosity);
        _total_luminosity += UV_luminosity;
        _spectrum_index.push_back(pop_from_temp(teff[i]));
        if (pop_from_temp(teff[i]) == -1) {
          std::cout << "UH OH NOT IMPLEMENTED RIGHT" << std::endl;
        }
      }
    }
  }


  if (_log) {
    _log->write_status("Succesfully read in photon sources from \"", filename,
                       "\".");
    _log->write_status("Found ", _positions.size(),
                       " active sources, with a total luminosity of ",
                       _total_luminosity, " s^-1.");
    if (_total_luminosity == 0.) {
      _log->write_warning("Total luminosity of sources in snapshot is zero!");
    }
  }

  // we are done with the UVLuminosityFunction, delete it
  delete luminosity_function;
}

/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - filename: Name of the snapshot file to read (required)
 *  - formation time name: Name of the data field that contains the stellar
 *    formation time values (default: FormationTime)
 *  - fallback unit length: Length unit to use if no units can be found in the
 *    snapshot file (default: 1. m, with warning)
 *  - fallback unit time: Time unit to use if no units can be found in the
 *    snapshot file (default: 1. s, with warning)
 *  - fallback unit mass: Mass unit to use if no units can be found in the
 *    snapshot file (default: 1. kg, with warning)
 *  - cutoff age: Upper limit for the age of star particles that emit UV
 *    radiation (default: 5. Myr)
 *  - use gas: Do gas particles contain stars (default: false)?
 *  - SFR unit: Unit for SFR values in the snapshot (default: (mass unit)/(time
 *    unit))
 *  - comoving integration flag: Was comoving integration active in the original
 *    simulation (default: false)?
 *  - hubble parameter: Reduced Hubble parameter used for the original
 *    simulation (default: 0.7)
 *
 * @param params ParameterFile to read from.
 * @param log Log to write logging information to.
 */
ArepoSnapshotPhotonSourceDistribution::ArepoSnapshotPhotonSourceDistribution(
    ParameterFile &params, Log *log)
    : ArepoSnapshotPhotonSourceDistribution(
          params.get_filename("PhotonSourceDistribution:filename"),
          params.get_value< std::string >(
              "PhotonSourceDistribution:formation time name", "FormationTime"),
          Box<>(params.get_physical_vector< QUANTITY_LENGTH >(
                    "SimulationBox:anchor"),
                params.get_physical_vector< QUANTITY_LENGTH >(
                    "SimulationBox:sides")),
          UVLuminosityFunctionFactory::generate(params, log),
          params.get_physical_value< QUANTITY_LENGTH >(
              "PhotonSourceDistribution:fallback unit length", "0. m"),
          params.get_physical_value< QUANTITY_TIME >(
              "PhotonSourceDistribution:fallback unit time", "0. s"),
          params.get_physical_value< QUANTITY_MASS >(
              "PhotonSourceDistribution:fallback unit mass", "0. kg"),
          params.get_physical_value< QUANTITY_TIME >(
              "PhotonSourceDistribution:cutoff age", "5. Myr"),
          params.get_value< bool >("PhotonSourceDistribution:use gas", false),
          params.get_physical_value< QUANTITY_MASS_RATE >(
              "PhotonSourceDistribution:SFR unit", "0. kg s^-1"),
          params.get_value< bool >(
              "PhotonSourceDistribution:comoving integration flag", false),
          params.get_value< double >(
              "PhotonSourceDistribution:hubble parameter", 0.7),
          log) {}

/**
 * @brief Get the number of sources in the snapshot file.
 *
 * @return Number of sources.
 */
photonsourcenumber_t
ArepoSnapshotPhotonSourceDistribution::get_number_of_sources() const {
  return _positions.size();
}

/**
 * @brief Get the position of one of the sources.
 *
 * Note that this function can alter the internal state of the
 * PhotonSourceDistribution, as for some implementations the positions are
 * decided randomly based on a RandomGenerator.
 *
 * @param index Valid index of a source, must be an integer in between 0 and
 * get_number_of_sources().
 * @return Position of the given source (in m).
 */
CoordinateVector<> ArepoSnapshotPhotonSourceDistribution::get_position(
    photonsourcenumber_t index) {
  return _positions[index];
}




/**
 * @brief Get the weight of one of the sources.
 *
 * At the moment, all sources have an equal weight.
 *
 * @param index Valid index of a source, must be an integer in between 0 and
 * get_number_of_sources().
 * @return Weight of the given source.
 */
double ArepoSnapshotPhotonSourceDistribution::get_weight(
    photonsourcenumber_t index) const {
  return _luminosities[index] / _total_luminosity;
}

/**
 * @brief Get the total luminosity of all sources in the snapshot file.
 *
 * @return Total luminosity (in s^-1).
 */
double ArepoSnapshotPhotonSourceDistribution::get_total_luminosity() const {
  return _total_luminosity;
}



double ArepoSnapshotPhotonSourceDistribution::get_photon_frequency(RandomGenerator &random_generator,
  photonsourcenumber_t index) {




  return _all_spectra[_spectrum_index[index]]->get_random_frequency(random_generator,0.0);

}
