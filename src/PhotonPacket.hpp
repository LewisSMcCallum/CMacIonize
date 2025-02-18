/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017, 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file PhotonPacket.hpp
 *
 * @brief Photon packet.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONPACKET_HPP
#define PHOTONPACKET_HPP

/*! @brief Size of the MPI buffer necessary to store a single Photon. */
#define PHOTON_MPI_SIZE ((9 + NUMBER_OF_IONNAMES) * sizeof(double))

#include "Configuration.hpp"
#include "CoordinateVector.hpp"
#include "ElementNames.hpp"
#include "PhotonType.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/**
 * @brief Photon packet.
 */
class PhotonPacket {
private:
  /*! @brief Current position of the photon packet (in m). */
  CoordinateVector<> _position;

  /*! @brief Propagation direction of the photon packet. */
  CoordinateVector<> _direction;

  /*! @brief Target optical depth for the photon packet. */
  double _target_optical_depth;

  /*! @brief Photoionization cross section of the photons in the photon packet
   *  (in m^2). Note that for ions other than hydrogen, these values contain an
   *  additional abundance factor. */
  double _photoionization_cross_section[NUMBER_OF_IONNAMES];

  double _dust_opacity;


  /*! @brief Energy of the photon packet (in Hz). */
  double _energy;

  /*! @brief Weight of the photon packet. */
  double _weight;

  /*! @brief Type of the photon. All photons start off as PHOTONTYPE_PRIMARY,
   *  but their type can change during reemission events. */
  PhotonType _type;

  /*! @brief Tracer for the number of scatterings the photon experiences */
  uint_fast32_t _scatter_counter;

  double _distance_travelled = 0.0;

public:
#ifdef HAVE_MPI
  /**
   * @brief Store the contents of the PhotonPacket in the given MPI
   * communication buffer.
   *
   * @param buffer Buffer to use (should be preallocated and have at least size
   * PHOTON_MPI_SIZE).
   */
  inline void pack(char buffer[PHOTON_MPI_SIZE]) {
    int_least32_t buffer_position = 0;
    MPI_Pack(&_position[0], 3, MPI_DOUBLE, buffer, PHOTON_MPI_SIZE,
             &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&_direction[0], 3, MPI_DOUBLE, buffer, PHOTON_MPI_SIZE,
             &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&_target_optical_depth, 1, MPI_DOUBLE, buffer, PHOTON_MPI_SIZE,
             &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(_photoionization_cross_section, NUMBER_OF_IONNAMES, MPI_DOUBLE,
             buffer, PHOTON_MPI_SIZE, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&_energy, 1, MPI_DOUBLE, buffer, PHOTON_MPI_SIZE, &buffer_position,
             MPI_COMM_WORLD);
    MPI_Pack(&_weight, 1, MPI_DOUBLE, buffer, PHOTON_MPI_SIZE, &buffer_position,
             MPI_COMM_WORLD);
  }

  /**
   * @brief Copy the contents of the given MPI communication buffer into this
   * PhotonPacket.
   *
   * @param buffer MPI communication buffer (should have at least size
   * PHOTON_MPI_SIZE).
   */
  inline void unpack(char buffer[PHOTON_MPI_SIZE]) {
    int_least32_t buffer_position = 0;
    double temp[3];
    MPI_Unpack(buffer, PHOTON_MPI_SIZE, &buffer_position, temp, 3, MPI_DOUBLE,
               MPI_COMM_WORLD);
    _position[0] = temp[0];
    _position[1] = temp[1];
    _position[2] = temp[2];
    MPI_Unpack(buffer, PHOTON_MPI_SIZE, &buffer_position, temp, 3, MPI_DOUBLE,
               MPI_COMM_WORLD);
    _direction[0] = temp[0];
    _direction[1] = temp[1];
    _direction[2] = temp[2];
    MPI_Unpack(buffer, PHOTON_MPI_SIZE, &buffer_position,
               &_target_optical_depth, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer, PHOTON_MPI_SIZE, &buffer_position,
               _photoionization_cross_section, NUMBER_OF_IONNAMES, MPI_DOUBLE,
               MPI_COMM_WORLD);
    MPI_Unpack(buffer, PHOTON_MPI_SIZE, &buffer_position, &_energy, 1,
               MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer, PHOTON_MPI_SIZE, &buffer_position, &_weight, 1,
               MPI_DOUBLE, MPI_COMM_WORLD);
  }
#endif

  /**
   * @brief Check that the given PhotonPacket is equal to this one.
   *
   * @param other Other PhotonPacket.
   */
  inline void check_equal(const PhotonPacket &other) {

    cmac_assert_message(_position[0] == other._position[0],
                        "Positions do not match!");
    cmac_assert_message(_position[1] == other._position[1],
                        "Positions do not match!");
    cmac_assert_message(_position[2] == other._position[2],
                        "Positions do not match!");
    cmac_assert_message(_direction[0] == other._direction[0],
                        "Directions do not match!");
    cmac_assert_message(_direction[1] == other._direction[1],
                        "Directions do not match!");
    cmac_assert_message(_direction[2] == other._direction[2],
                        "Directions do not match!");
    cmac_assert_message(_target_optical_depth == other._target_optical_depth,
                        "Target optical depths do not match!");
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      cmac_assert_message(_photoionization_cross_section[i] ==
                              other._photoionization_cross_section[i],
                          "Cross sections do not match!");
    }
    cmac_assert_message(_energy == other._energy, "Energies do not match!");
    cmac_assert_message(_weight == other._weight, "Weights do not match!");
  }

  /**
   * @brief Get a constant reference to the direction of the photon packet.
   *
   * @return Constant reference to the direction.
   */
  inline const CoordinateVector<> &get_direction() const { return _direction; }

  /**
   * @brief Get a reference to the direction of the photon packet.
   *
   * @return Reference to the direction.
   */
  inline CoordinateVector<> &get_direction() { return _direction; }

  /**
   * @brief Set the direction vector directly.
   *
   * @param direction Direction for the photon packet.
   */
  inline void set_direction(const CoordinateVector<> direction) {
    _direction = direction;
    // make sure the direction is properly normalised
    _direction /= _direction.norm();
  }

  /**
   * @brief Get a constant reference to the position.
   *
   * @return Constant reference to the position (in m).
   */
  inline const CoordinateVector<> &get_position() const { return _position; }

  /**
   * @brief Update the position of the photon packet.
   *
   * @param position Position for the photon packet (in m).
   */
  inline void set_position(const CoordinateVector<> position) {
    _position = position;
  }

  /**
   * @brief Get the target optical depth for the photon packet.
   *
   * @return Target optical depth.
   */
  inline double get_target_optical_depth() const {
    return _target_optical_depth;
  }

  /**
   * @brief Set the target optical depth for the photon packet.
   *
   * @param target_optical_depth New target optical depth.
   */
  inline void set_target_optical_depth(const double target_optical_depth) {
    _target_optical_depth = target_optical_depth;
  }

  /**
   * @brief Get the photoionization cross section for the given ion for the
   * photon packet.
   *
   * @param ion Ion for which we want the photoionization cross section.
   * @return Corresponding photoionization cross section for the photon packet
   * (in m^2).
   */
  inline double
  get_photoionization_cross_section(const int_fast32_t ion) const {
    return _photoionization_cross_section[ion];
  }


  inline double get_dust_opacity() const {return _dust_opacity;}


  inline void set_distance_travelled(const double value) {
    _distance_travelled = value;
  }
  inline double get_distance_travelled() const {return _distance_travelled;}

  inline void add_distance_travelled(const double distance) {
    _distance_travelled+= distance;
  }


  /**
   * @brief Set the photoionization cross section for the given ion for the
   * photon packet.
   *
   * @param ion Ion for which we want to set the photoionization cross section.
   * @param photoionization_cross_section New photoionization cross section
   * (in m^2).
   */
  inline void set_photoionization_cross_section(
      const int_fast32_t ion, const double photoionization_cross_section) {
    _photoionization_cross_section[ion] = photoionization_cross_section;
  }

  inline void set_dust_opacity(const double opacity) {
    _dust_opacity = opacity;
  }


  /**
   * @brief Get the energy of the photon packet.
   *
   * @return Energy of the photon packet (in Hz).
   */
  inline double get_energy() const { return _energy; }

  /**
   * @brief Set the energy of the photon packet.
   *
   * @param energy Energy of the photon packet (in Hz).
   */
  inline void set_energy(const double energy) { _energy = energy; }

  /**
   * @brief Get the weight for the photon packet.
   *
   * @return Weight for the photon packet.
   */
  inline double get_weight() const { return _weight; }

  /**
   * @brief Set the weight for the photon packet.
   *
   * @param weight New weight for the photon packet.
   */
  inline void set_weight(const double weight) { _weight = weight; }

  /**
   * @brief Get the type of the photon.
   *
   * @return PhotonType type identifier.
   */
  inline PhotonType get_type() const { return _type; }

  /**
   * @brief Set the photon type.
   *
   * @param type PhotonType type identifier.
   */
  inline void set_type(PhotonType type) { _type = type; }

  /**
   * @brief Get the photon scatter counter.
   *
   * @return scatter counter for a photon packet
   *
   */

  inline uint_fast32_t get_scatter_counter() const { return _scatter_counter; }

  /**
   * @brief Set the photon scatter counter.
   *
   * @param scatter_counter sets a scatter counter for a photon packet
   */
  inline void set_scatter_counter(uint_fast32_t scatter_counter) {
    _scatter_counter = scatter_counter;
  }
};

#endif // PHOTONPACKET_HPP
