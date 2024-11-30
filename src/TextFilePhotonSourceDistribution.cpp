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
#include "Pegase3PhotonSourceSpectrum.hpp"
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


// Function to perform linear interpolation for descending xVals
double interpolate(double x, const std::vector<double>& xVals, const std::vector<double>& yVals) {
    // Ensure inputs are valid
    if (xVals.size() != yVals.size() || xVals.empty()) {
        throw std::invalid_argument("Invalid input: xVals and yVals must have the same size and cannot be empty.");
    }


    if (x > xVals.front()){
      return yVals.front();
    } else if (x < xVals.back()) {
      return yVals.back();
    }

    // Find the interval containing x
    for (size_t i = 0; i < xVals.size() - 1; ++i) {
        if (xVals[i] >= x && x >= xVals[i + 1]) {
            // Perform linear interpolation
            double x1 = xVals[i];
            double x2 = xVals[i + 1];
            double y1 = yVals[i];
            double y2 = yVals[i + 1];
            return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
        }
    }

    // If we reach here, x was not in any valid interval
    throw std::logic_error("Could not interpolate: x is not within any interval.");
}


size_t findClosestIndex(double value, const std::vector<double>& values) {
    if (values.empty()) {
        throw std::invalid_argument("The values vector cannot be empty.");
    }

    size_t closestIndex = 0;
    double minDifference = std::abs(value - values[0]);

    for (size_t i = 1; i < values.size(); ++i) {
        double difference = std::abs(value - values[i]);
        if (difference < minDifference) {
            minDifference = difference;
            closestIndex = i;
        }
    }

    return closestIndex;
}

TextFilePhotonSourceDistribution::TextFilePhotonSourceDistribution(
    std::string filename, double time, Log *log)
    : _log(log) {





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

    std::vector<double> stellarMasses = {57.95, 46.94, 38.08, 34.39, 30.98, 28.0, 25.29, 22.90, 20.76, 18.80, 17.08, 15.55};
    std::vector<double> temperatures = {44852, 42857, 40862, 39865, 38867, 37870, 36872, 35874, 34877, 33879, 32882, 31884};

    std::vector<double> avail_temps = {32000,34000,34000,35000,36000,37000,39000,39000,40000,41000,42000,43000,44000,45000};

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
          _spectrum_index.push_back(14);
        } else {
          double interpolatedTemp = interpolate(mass, stellarMasses, temperatures);
          size_t closestIndex = findClosestIndex(interpolatedTemp, avail_temps);
          std::cout << "Adding star of mass " << mass << " temp of " << interpolatedTemp << " for spec index " << closestIndex << std::endl;
          _spectrum_index.push_back(closestIndex);
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
