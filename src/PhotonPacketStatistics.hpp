/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file PhotonPacketStatistics.hpp
 *
 * @brief Statistical information about the photon packets and re-emission
 * events.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONPACKETSTATISTICS_HPP
#define PHOTONPACKETSTATISTICS_HPP

#include "AtomicValue.hpp"
#include "ParameterFile.hpp"
#include "PhotonPacket.hpp"
#include <fstream>
#include <vector>

/**
 * @brief Statistical information about the photon packets and re-emission
 * events.
 */
class PhotonPacketStatistics {
private:
  /**
   * @brief variable to store the number of photons on bins according to number
   * of scatterings.
   */
  std::vector< AtomicValue< uint_fast64_t > > _scatter_histogram;


  AtomicValue<uint_fast32_t> _num_abs;
  AtomicValue<uint_fast32_t> _num_escape;

  AtomicValue<uint_fast32_t> _num_abs_dust;
  


  AtomicValue<uint_fast32_t> _num_abs_dens;
  AtomicValue<uint_fast32_t> _num_abs_dif;

  const uint_fast32_t _track_source_num;

  const uint_fast32_t numbins = 1000;

  std::vector<double> _frequencies;
  std::vector<AtomicValue<uint_fast32_t>> _ingoing_spectrum;
  std::vector<AtomicValue<uint_fast32_t>> _outgoing_spectrum;
  const double min_frequency = 3.289e15;
  const double max_frequency = 4. * min_frequency;
  const double _bin_width = (max_frequency - min_frequency) / numbins;


public:
  /**
   * @brief constructor
   *
   * @param max_scatter maximum number of scatters recorded by the statistics
   */
  inline PhotonPacketStatistics(uint_fast32_t max_scatter, uint_fast32_t track_source_num)
      : _scatter_histogram(max_scatter + 2), _track_source_num(track_source_num),
       _frequencies(numbins),_ingoing_spectrum(numbins),
      _outgoing_spectrum(numbins){
        _num_abs.set(0);
        _num_escape.set(0);
        _num_abs_dens.set(0);
        _num_abs_dif.set(0);

        for (uint_fast32_t i = 0; i < numbins; ++i) {
          _frequencies[i] =
            min_frequency + i * (max_frequency - min_frequency) /
                            (numbins - 1.);
          _ingoing_spectrum[i].set(0.0);
          _outgoing_spectrum[i].set(0.0);
        }
      }
  /**
   * @brief parameter file constructor
   *
   * These are the parameters that are used by this function:
   *   - maximum number of scatters; (default: 5)
   * @param params reference to the parameter file
   */
  inline PhotonPacketStatistics(ParameterFile &params)
      : PhotonPacketStatistics(params.get_value< uint_fast32_t >(
            "PhotonPacketStatistics:maximum number of scatters", 5),
            params.get_value< uint_fast32_t >(
            "PhotonPacketStatistics:track source index", 0)) {}
  /**
   * @brief function that implements photon absorption termination
   *
   * @param packet photon packet retrieved to be terminated
   */
  inline void absorb_photon(bool dense) {
    _num_abs.pre_increment();
    if (dense){
      _num_abs_dens.pre_increment();
    } else {
      _num_abs_dif.pre_increment();
    }
  }

  inline void absorb_photon_dust() {
      _num_abs_dust.pre_increment();
  }
  inline void absorb_photon(const PhotonPacket &packet) {
    _num_abs.pre_increment();
    size_t scatter_counter = packet.get_scatter_counter();
    _scatter_histogram[std::min(scatter_counter, _scatter_histogram.size() - 1)]
        .pre_increment();
  }

  inline void injected_photon(const PhotonPacket &packet) {
    double photon_energy = packet.get_energy();
    if (photon_energy >= min_frequency && photon_energy < max_frequency) {
      // Determine the bin index for this energy
      uint_fast32_t bin_index = static_cast<uint_fast32_t>((photon_energy - min_frequency) / _bin_width);
      // Increment the count in the corresponding bin
      _ingoing_spectrum[bin_index].pre_increment();
    }
  }

  inline uint_fast32_t get_tracked_source() {return _track_source_num;}


  /**
   * @brief Function that implements photon escape termination
   *
   * @param packet Photon packet.
   */
  inline void escape_photon(const PhotonPacket &packet) {
    _num_escape.pre_increment();
    size_t scatter_counter = packet.get_scatter_counter();
    _scatter_histogram[std::min(scatter_counter, _scatter_histogram.size() - 1)]
        .pre_increment();
    double photon_energy = packet.get_energy();
    if (photon_energy >= min_frequency && photon_energy < max_frequency) {
      // Determine the bin index for this energy
      uint_fast32_t bin_index = static_cast<uint_fast32_t>((photon_energy - min_frequency) / _bin_width);
      // Increment the count in the corresponding bin
      _outgoing_spectrum[bin_index].pre_increment();
    }
  }



  inline uint_fast32_t get_num_escaped() {
    return _num_escape.value();
  }

  inline uint_fast32_t get_num_absorbed() {
    return _num_abs.value();
  }

  inline uint_fast32_t get_num_abs_dens() {
    return _num_abs_dens.value();
  }

  inline uint_fast32_t get_num_abs_dif() {
    return _num_abs_dif.value();
  }

  inline uint_fast32_t get_num_abs_dust() {
    return _num_abs_dust.value();
  }


  inline void reset_counters() {
    _num_abs.set(0);
    _num_escape.set(0);
    _num_abs_dust.set(0);
    _num_abs_dens.set(0);
    _num_abs_dif.set(0);
    for (uint_fast32_t i = 0; i < numbins; ++i) {
      _ingoing_spectrum[i].set(0.0);
      _outgoing_spectrum[i].set(0.0);
    }
  }
  /**
   * @brief function that outputs re-emission statistics of photons
   */
  inline void print_stats() {
// histogram thing
  {
    std::ofstream output_stats("photon_statistics.txt");
    output_stats << "# Scattering statistics for photons\n";
    output_stats << "# Nscatter\t BinCount  \n";
    for (uint_fast32_t i = 0; i < _scatter_histogram.size(); i++) {
      output_stats << i << "\t" << _scatter_histogram[i].value() << "\n";
    }
  }
// outgoing spectrum
  {
    std::ofstream output_stats("leaking_spectrum.txt");
    output_stats << "# Frequency\t Spectrum \n";
    for (uint_fast32_t i = 0; i < _frequencies.size(); i++) {
      output_stats << _frequencies[i] << "\t" << _outgoing_spectrum[i].value() << "\n";
    }
  }  
//ingoing spectrum
  {
    std::ofstream output_stats("input_spectrum.txt");
    output_stats << "# Frequency\t Spectrum \n";
    for (uint_fast32_t i = 0; i < _frequencies.size(); i++) {
      output_stats << _frequencies[i] << "\t" << _ingoing_spectrum[i].value() << "\n";
    }
  }

}
};

#endif // PHOTONPACKETSTATISTICS_HPP
