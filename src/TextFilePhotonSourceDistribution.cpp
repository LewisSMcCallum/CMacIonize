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
 * @file TextFilePhotonSourceDistribution.cpp
 *
 * @brief TextFilePhotonSourceDistribution implementation.
 *
 * @author Lewis McCallum (lm261@st-andrews.ac.uk)
 */
#include "TextFilePhotonSourceDistribution.hpp"
#include "WMBasicPhotonSourceSpectrum.hpp"
#include "PowerLawPhotonSourceSpectrum.hpp"
#include "CoordinateVector.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
//#include "UnitConverter.hpp"

#include <cinttypes>
#include <fstream>

/**
 * @brief Constructor.
 *
 * @param filename Name of the text file to read.
 * @param log Log to write logging info to.
 */
TextFilePhotonSourceDistribution::TextFilePhotonSourceDistribution(
    std::string filename, double time, Log *log)
    : _log(log) {





  _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(40000,25,log));
  _all_spectra.push_back(new PowerLawPhotonSourceSpectrum(4,log));



  std::ifstream file;
  file.open(filename);
  if (!file.is_open()) {
    cmac_error("Could not open file \"%s\"!", filename.c_str());
  } else {
    std::cout << "Opened file - " << filename << std::endl;
  }


  double time_val,posx,posy,posz,luminosity,mass;
  int event,index;


  std::string dummyLine,star_type;

  std::getline(file, dummyLine);

  time_val = 0.0;

  while (!file.eof() && time_val <= time) {
    file >> time_val >> posx >> posy >> posz >> event >> index >> luminosity >> mass >> star_type;

    if (event == 2) {
      _to_delete.push_back(index);
    }
  }



  file.close();

  file.open(filename);


  std::getline(file, dummyLine);


  time_val = 0.0;
  while (!file.eof() && time_val <= time) {
    file >> time_val >> posx >> posy >> posz >> event >> index >> luminosity >> mass >> star_type;
    if (event == 1) {
      if (std::find(_to_delete.begin(), _to_delete.end(), index) == _to_delete.end() && luminosity > 0.0) {
        _positions.push_back(CoordinateVector<double>(posx,posy,posz));

        _luminosities.push_back(luminosity);
        if (star_type == "HOLMES") {
          _spectrum_index.push_back(1);
        } else {
          _spectrum_index.push_back(0);
        }
      }

    }

  }

  int number_of_positions = _luminosities.size();


 _total_luminosity = 0.0;

 for (int i=0;i<number_of_positions;i++) {
   _total_luminosity += _luminosities[i];
 }








  file.close();

  if (_log) {
    _log->write_status("TextFilePhotonSourceDistribution with ",
                       number_of_positions, " sources and total luminosity of ",
                     _total_luminosity, " s^-1.");
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - filename: Name of the file (default: sinks.txt)
 *
 * @param params ParameterFile to read from.
 * @param log Log to write logging info to.
 */
TextFilePhotonSourceDistribution::TextFilePhotonSourceDistribution(
    ParameterFile &params, Log *log)
    : TextFilePhotonSourceDistribution(
          params.get_value< std::string >("PhotonSourceDistribution:filename",
                                          "sinks.txt"),
          params.get_physical_value <QUANTITY_TIME>("PhotonSourceDistribution:time", "100 Myr"),
          log) {}

/**
 * @brief Get the number of sources in the ASCII file.
 *
 * @return Number of sources.
 */
photonsourcenumber_t
TextFilePhotonSourceDistribution::get_number_of_sources() const {
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
CoordinateVector<> TextFilePhotonSourceDistribution::get_position(
    photonsourcenumber_t index) {
  return _positions[index];
}

/**
 * @brief Get the weight of one of the sources.
 *
 * @param index Valid index of a source, must be an integer in between 0 and
 * get_number_of_sources().
 * @return Weight of the given source.
 */
double TextFilePhotonSourceDistribution::get_weight(
    photonsourcenumber_t index) const {
  return _luminosities[index]/_total_luminosity;
}

double TextFilePhotonSourceDistribution::get_photon_weighting(photonsourcenumber_t index) const {


  return _luminosities[index]/_total_luminosity;
}

/**
 * @brief Get the total luminosity of all sources in the ASCII file.
 *
 * @return Total luminosity (in s^-1).
 */
double TextFilePhotonSourceDistribution::get_total_luminosity() const {
  return _total_luminosity;
}


double TextFilePhotonSourceDistribution::get_photon_frequency(RandomGenerator &random_generator,
  photonsourcenumber_t index) {




  return _all_spectra[_spectrum_index[index]]->get_random_frequency(random_generator,0.0);

}
