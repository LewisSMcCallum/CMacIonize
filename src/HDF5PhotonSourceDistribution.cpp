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
 * @file HDF5PhotonSourceDistribution.cpp
 *
 * @brief HDF5PhotonSourceDistribution implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "HDF5PhotonSourceDistribution.hpp"
#include "HDF5Tools.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "UVLuminosityFunctionFactory.hpp"
#include "UnitConverter.hpp"
#include "WMBasicPhotonSourceSpectrum.hpp"
#include "Pegase3PhotonSourceSpectrum.hpp"
#include "RestartReader.hpp"
#include "RestartWriter.hpp"
#include "DensitySubGridCreator.hpp"
#include "SupernovaHandler.hpp"

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
HDF5PhotonSourceDistribution::HDF5PhotonSourceDistribution(
    std::string filename, const Box<> box, const double update_interval, const bool has_lifetimes,
    Log *log)
    : _log(log),_update_interval(update_interval),_has_lifetimes(has_lifetimes) {



   _has_exploded=false;
   _number_of_updates = 0;


    _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(32000,25,log));
    _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(34000,25,log));
    _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(34000,25,log));
    _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(35000,40,log));
    _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(36000,25,log));
    _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(37000,25,log));
    _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(39000,25,log));
    _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(39000,25,log));
    _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(40000,25,log));
    _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(41000,40,log));
    _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(42000,40,log));
    _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(43000,40,log));
    _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(44000,40,log));
    _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(45000,40,log));
    _all_spectra.push_back(new Pegase3PhotonSourceSpectrum(1e10,0.02,log));

  // turn off default HDF5 error handling: we catch errors ourselves
  HDF5Tools::initialize();

  // open the file for reading
  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_READ);



  _total_luminosity = 0.;

  // open the group containing the star particle data
  HDF5Tools::HDF5Group maingroup =
      HDF5Tools::open_group(file, "/PartType0");
  // read the positions
  std::vector< CoordinateVector<> > positions =
      HDF5Tools::read_dataset< CoordinateVector<> >(maingroup,
                                                      "Sources");





   std::vector< double > read_lums =
       HDF5Tools::read_dataset< double > (maingroup, "SourceLuminosities");



  _spectrum_index = HDF5Tools::read_dataset<int> (maingroup, "spec_index");

  if (_has_lifetimes) {
    _source_lifetimes = HDF5Tools::read_dataset<double> (maingroup, "Lifetimes");
  }
  


  // close the group
  HDF5Tools::close_group(maingroup);
  // close the file
  HDF5Tools::close_file(file);


  novahandler = new SupernovaHandler(1.e44);


  for (size_t i = 0; i < positions.size(); ++i) {
    CoordinateVector<> position = positions[i];
    if (box.inside(position)) {
      double UV_luminosity = read_lums[i];
      
      _positions.push_back(position);
      _luminosities.push_back(UV_luminosity);
      _total_luminosity += UV_luminosity;
      
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
HDF5PhotonSourceDistribution::HDF5PhotonSourceDistribution(
    ParameterFile &params, Log *log)
    : HDF5PhotonSourceDistribution(
          params.get_filename("PhotonSourceDistribution:filename"),
          Box<>(params.get_physical_vector< QUANTITY_LENGTH >(
                    "SimulationBox:anchor"),
                params.get_physical_vector< QUANTITY_LENGTH >(
                    "SimulationBox:sides")),
          params.get_physical_value<QUANTITY_TIME>("PhotonSourceDistribution:update interval", "0.05 Myr"),
          params.get_value<bool>("PhotonSourceDistribution:has lifetimes", false),
          log) {
          }

/**
 * @brief Get the number of sources in the snapshot file.
 *
 * @return Number of sources.
 */
photonsourcenumber_t
HDF5PhotonSourceDistribution::get_number_of_sources() const {
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
CoordinateVector<> HDF5PhotonSourceDistribution::get_position(
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
double HDF5PhotonSourceDistribution::get_weight(
    photonsourcenumber_t index) const {
  return _luminosities[index] / _total_luminosity;
}

/**
 * @brief Get the total luminosity of all sources in the snapshot file.
 *
 * @return Total luminosity (in s^-1).
 */
double HDF5PhotonSourceDistribution::get_total_luminosity() const {
  return _total_luminosity;
}



double HDF5PhotonSourceDistribution::get_photon_frequency(RandomGenerator &random_generator,
  photonsourcenumber_t index) {

  return _all_spectra[_spectrum_index[index]]->get_random_frequency(random_generator,0.0);

}


  /**
   * @brief Write the distribution to the given restart file.
   *
   * @param restart_writer RestartWriter to use.
   */
  void HDF5PhotonSourceDistribution::write_restart_file(RestartWriter &restart_writer) const {

    restart_writer.write(_update_interval);
    restart_writer.write(_has_lifetimes);
    const size_t number_of_sources = _positions.size();
    restart_writer.write(number_of_sources);
    for (size_t i = 0; i < number_of_sources; ++i) {
      _positions[i].write_restart_file(restart_writer);
      restart_writer.write(_luminosities[i]);
      restart_writer.write(_spectrum_index[i]);
    }
    restart_writer.write(_total_luminosity);
  }

  /**
   * @brief Restart constructor.
   *
   * @param restart_reader Restart file to read from.
   */
  HDF5PhotonSourceDistribution::HDF5PhotonSourceDistribution(RestartReader &restart_reader):
  _update_interval(restart_reader.read< double >()),_has_lifetimes(restart_reader.read< bool >()) {

    const size_t number_of_sources = restart_reader.read< size_t >();
    _positions.resize(number_of_sources);
    _luminosities.resize(number_of_sources);
    for (size_t i = 0; i < number_of_sources; ++i) {
      _positions[i] = CoordinateVector<>(restart_reader);
      _luminosities[i] = restart_reader.read< double >();
      _spectrum_index[i] = restart_reader.read<int>();
    }
    _total_luminosity = restart_reader.read< double >();
    
   

    _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(40000,25,nullptr));
    _all_spectra.push_back(new Pegase3PhotonSourceSpectrum(1e10,0.02,nullptr));

    novahandler = new SupernovaHandler(1.e44);
  }





 bool HDF5PhotonSourceDistribution::do_stellar_feedback(const double current_time) const {
    return (_to_do_feedback.size() > 0);
  }



   void HDF5PhotonSourceDistribution::get_sne_radii(DensitySubGridCreator< HydroDensitySubGrid > &grid_creator) {

       for (uint_fast32_t i = 0; i < _to_do_feedback.size(); ++i) {

        double r_inj,r_st,nbar,num_inj;


        std::tie(r_inj,r_st,nbar,num_inj) = novahandler->get_r_inj(&grid_creator,_to_do_feedback[i]);

         _r_inj.push_back(r_inj);
         _r_st.push_back(r_st);
         _nbar.push_back(nbar);
         _num_cells_injected.push_back(num_inj);


       }
  }



  void HDF5PhotonSourceDistribution::add_stellar_feedback(HydroDensitySubGrid &subgrid, Hydro &hydro) {



    for (uint_fast32_t i = 0; i < _to_do_feedback.size(); ++i) {



      novahandler->inject_sne(subgrid, hydro, _to_do_feedback[i], _r_inj[i],_r_st[i],_nbar[i],_num_cells_injected[i]);


    }
  }

  void HDF5PhotonSourceDistribution::done_stellar_feedback() {

    for (uint_fast32_t i=0; i<_to_do_feedback.size();i++) {

    std::cout << "\n SNe INJECTION HERE: R_inj = " << _r_inj[i] << " R_st = " <<  _r_st[i]
       << " num_cells = " <<  _num_cells_injected[i] << " nbar = "  << _nbar[i] << "\n";
    }


    _has_exploded = true;
    _to_do_feedback.clear();
    _r_inj.clear();
    _r_st.clear();
    _num_cells_injected.clear();
    _nbar.clear();

  }


  bool HDF5PhotonSourceDistribution::update(DensitySubGridCreator< HydroDensitySubGrid > *grid_creator, double actual_timestep) {

    bool changed = false;


    // clear out sources which no longer exist and add them to SNe todo list
    size_t i = 0;
    while (i < _source_lifetimes.size()) {
      _source_lifetimes[i] -= actual_timestep;
      if (_source_lifetimes[i] <= 0.) {

        // remove the element
        _to_do_feedback.push_back(_positions[i]);
        _positions.erase(_positions.begin() + i);
        _source_lifetimes.erase(_source_lifetimes.begin() + i);
        _luminosities.erase(_luminosities.begin() + i);
        _spectrum_index.erase(_spectrum_index.begin() + i);
        changed = true;

      } else {
        // check the next element
        ++i;
      }
    }

    _number_of_updates += 1;


    return changed;
  }

