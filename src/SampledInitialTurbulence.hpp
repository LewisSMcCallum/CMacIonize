/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2026 Lewis McCallum
 *
 * CMacIonize is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *******************************************************************************/

/**
 * @file SampledInitialTurbulence.hpp
 *
 * @brief Small sampled set of solenoidal Fourier modes for one-time initial
 * turbulence.
 */
#ifndef SAMPLEDINITIALTURBULENCE_HPP
#define SAMPLEDINITIALTURBULENCE_HPP

#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "Log.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <random>
#include <set>
#include <tuple>
#include <vector>

/** @brief Efficiently evaluated sampled solenoidal turbulent field. */
class SampledInitialTurbulence {
private:
  CoordinateVector< int_fast32_t > _number_of_subgrids;
  CoordinateVector< int_fast32_t > _number_of_cells;
  std::vector< CoordinateVector<> > _amplitudes_real;
  std::vector< CoordinateVector<> > _amplitudes_imaginary;
  std::vector< double > _sin_x;
  std::vector< double > _sin_y;
  std::vector< double > _sin_z;
  std::vector< double > _cos_x;
  std::vector< double > _cos_y;
  std::vector< double > _cos_z;

  inline void get_subgrid_offsets(const uint_fast32_t index,
                                  int_fast32_t &offset_x,
                                  int_fast32_t &offset_y,
                                  int_fast32_t &offset_z) const {
    offset_x = index / (_number_of_subgrids.y() * _number_of_subgrids.z());
    offset_y =
        (index - offset_x * _number_of_subgrids.y() *
                     _number_of_subgrids.z()) /
        _number_of_subgrids.z();
    offset_z = index -
               offset_x * _number_of_subgrids.y() *
                   _number_of_subgrids.z() -
               offset_y * _number_of_subgrids.z();
  }

  inline CoordinateVector<>
  evaluate(const uint_fast32_t oix, const uint_fast32_t oiy,
           const uint_fast32_t oiz) const {
    CoordinateVector<> field;
    const uint_fast32_t number_of_modes = _amplitudes_real.size();
    for (uint_fast32_t imode = 0; imode < number_of_modes; ++imode) {
      const double cosx = _cos_x[oix + imode];
      const double cosy = _cos_y[oiy + imode];
      const double cosz = _cos_z[oiz + imode];
      const double sinx = _sin_x[oix + imode];
      const double siny = _sin_y[oiy + imode];
      const double sinz = _sin_z[oiz + imode];
      const double cosyz = cosy * cosz - siny * sinz;
      const double sinyz = siny * cosz + cosy * sinz;
      const double cosxyz = cosx * cosyz - sinx * sinyz;
      const double sinxyz = sinx * cosyz + cosx * sinyz;
      field += _amplitudes_real[imode] * cosxyz -
               _amplitudes_imaginary[imode] * sinxyz;
    }
    return field;
  }

public:
  SampledInitialTurbulence(
      const CoordinateVector< int_fast32_t > number_of_subgrids,
      const CoordinateVector< int_fast32_t > number_of_cells,
      const Box<> box, const double kmin, const double kmax,
      const double kpeak, const double concentration_factor,
      const uint_fast32_t number_of_modes, const int_fast32_t seed,
      Log *log = nullptr)
      : _number_of_subgrids(number_of_subgrids),
        _number_of_cells(number_of_cells) {
    if (!(kmin > 0.) || kmax < kmin) {
      cmac_error("Invalid initial turbulence wave-number interval: [%g, %g].",
                 kmin, kmax);
    }
    if (!(concentration_factor > 0.)) {
      cmac_error("Initial turbulence concentration factor must be positive.");
    }
    if (number_of_modes == 0) {
      cmac_error("Sampled initial turbulence requires at least one mode.");
    }

    const double mode_length =
        std::max(box.get_sides().x(),
                 std::max(box.get_sides().y(), box.get_sides().z()));
    const double inverse_mode_length = 1. / mode_length;
    const int_fast32_t ikmax = static_cast< int_fast32_t >(std::ceil(kmax));

    std::mt19937 generator(static_cast< std::mt19937::result_type >(seed));
    std::uniform_int_distribution< int_fast32_t > integer_mode(-ikmax, ikmax);
    std::uniform_real_distribution< double > uniform(0., 2. * M_PI);
    std::set< std::tuple< int_fast32_t, int_fast32_t, int_fast32_t > > used;
    std::vector< CoordinateVector<> > wavevectors;
    wavevectors.reserve(number_of_modes);
    _amplitudes_real.reserve(number_of_modes);
    _amplitudes_imaginary.reserve(number_of_modes);

    double selected_kmin = DBL_MAX;
    double selected_kmax = 0.;
    uint_fast32_t attempts = 0;
    const uint_fast32_t maximum_attempts = 1000000;
    while (wavevectors.size() < number_of_modes &&
           attempts < maximum_attempts) {
      ++attempts;
      int_fast32_t kx = integer_mode(generator);
      int_fast32_t ky = integer_mode(generator);
      int_fast32_t kz = integer_mode(generator);
      if (kx < 0 || (kx == 0 && ky < 0) ||
          (kx == 0 && ky == 0 && kz < 0)) {
        kx = -kx;
        ky = -ky;
        kz = -kz;
      }
      const double kk = static_cast< double >(kx * kx + ky * ky + kz * kz);
      if (!(kk > 0.)) {
        continue;
      }
      const double kmag = std::sqrt(kk);
      if (kmag < kmin || kmag > kmax) {
        continue;
      }
      const std::tuple< int_fast32_t, int_fast32_t, int_fast32_t > key(
          kx, ky, kz);
      if (!used.insert(key).second) {
        continue;
      }

      const CoordinateVector<> khat(kx / kmag, ky / kmag, kz / kmag);
      const CoordinateVector<> reference =
          std::abs(khat.z()) < 0.9 ? CoordinateVector<>(0., 0., 1.)
                                  : CoordinateVector<>(0., 1., 0.);
      CoordinateVector<> e1(
          khat.y() * reference.z() - khat.z() * reference.y(),
          khat.z() * reference.x() - khat.x() * reference.z(),
          khat.x() * reference.y() - khat.y() * reference.x());
      e1 /= std::sqrt(e1.norm2());
      const CoordinateVector<> e2(
          khat.y() * e1.z() - khat.z() * e1.y(),
          khat.z() * e1.x() - khat.x() * e1.z(),
          khat.x() * e1.y() - khat.y() * e1.x());

      const double kdiff = kmag - kpeak;
      const double weight =
          std::sqrt(std::exp(-kdiff * kdiff /
                             (concentration_factor * concentration_factor)) /
                    kk);
      const double phi = uniform(generator);
      const double theta1 = uniform(generator);
      const double theta2 = uniform(generator);
      const double ga = std::sin(phi);
      const double gb = std::cos(phi);
      _amplitudes_real.push_back(
          weight * (std::cos(theta1) * ga * e1 +
                    std::cos(theta2) * gb * e2));
      _amplitudes_imaginary.push_back(
          weight * (std::sin(theta1) * ga * e1 +
                    std::sin(theta2) * gb * e2));
      wavevectors.push_back(CoordinateVector<>(
          kx * inverse_mode_length, ky * inverse_mode_length,
          kz * inverse_mode_length));
      selected_kmin = std::min(selected_kmin, kmag);
      selected_kmax = std::max(selected_kmax, kmag);
    }
    if (wavevectors.size() != number_of_modes) {
      cmac_error("Could only sample %zu of %u requested initial turbulence "
                 "modes in [%g, %g].",
                 wavevectors.size(),
                 static_cast< unsigned int >(number_of_modes), kmin, kmax);
    }

    const CoordinateVector< int_fast32_t > ntot(
        number_of_subgrids.x() * number_of_cells.x(),
        number_of_subgrids.y() * number_of_cells.y(),
        number_of_subgrids.z() * number_of_cells.z());
    const CoordinateVector<> dx(box.get_sides().x() / ntot.x(),
                                box.get_sides().y() / ntot.y(),
                                box.get_sides().z() / ntot.z());
    const CoordinateVector<> anchor = box.get_anchor();
    _sin_x.resize(number_of_modes * ntot.x());
    _sin_y.resize(number_of_modes * ntot.y());
    _sin_z.resize(number_of_modes * ntot.z());
    _cos_x.resize(number_of_modes * ntot.x());
    _cos_y.resize(number_of_modes * ntot.y());
    _cos_z.resize(number_of_modes * ntot.z());
    for (uint_fast32_t imode = 0; imode < number_of_modes; ++imode) {
      for (int_fast32_t ix = 0; ix < ntot.x(); ++ix) {
        const double x = anchor.x() + (ix + 0.5) * dx.x();
        const size_t index = ix * number_of_modes + imode;
        const double angle = 2. * M_PI * wavevectors[imode].x() * x;
        _sin_x[index] = std::sin(angle);
        _cos_x[index] = std::cos(angle);
      }
      for (int_fast32_t iy = 0; iy < ntot.y(); ++iy) {
        const double y = anchor.y() + (iy + 0.5) * dx.y();
        const size_t index = iy * number_of_modes + imode;
        const double angle = 2. * M_PI * wavevectors[imode].y() * y;
        _sin_y[index] = std::sin(angle);
        _cos_y[index] = std::cos(angle);
      }
      for (int_fast32_t iz = 0; iz < ntot.z(); ++iz) {
        const double z = anchor.z() + (iz + 0.5) * dx.z();
        const size_t index = iz * number_of_modes + imode;
        const double angle = 2. * M_PI * wavevectors[imode].z() * z;
        _sin_z[index] = std::sin(angle);
        _cos_z[index] = std::cos(angle);
      }
    }

    if (log) {
      log->write_status("Sampled initial turbulence: ", number_of_modes,
                        " solenoidal modes; selected |k| range ",
                        selected_kmin, "--", selected_kmax,
                        "; requested peak ", kpeak, ".");
    }
  }

  inline void
  get_turbulent_field(const uint_fast32_t index,
                      std::vector< CoordinateVector<> > &field) const {
    int_fast32_t offset_x, offset_y, offset_z;
    get_subgrid_offsets(index, offset_x, offset_y, offset_z);
    const uint_fast32_t number_of_modes = _amplitudes_real.size();
    field.resize(static_cast< size_t >(_number_of_cells.x()) *
                 _number_of_cells.y() * _number_of_cells.z());
    size_t cell_index = 0;
    for (int_fast32_t ix = 0; ix < _number_of_cells.x(); ++ix) {
      const uint_fast32_t oix =
          (offset_x * _number_of_cells.x() + ix) * number_of_modes;
      for (int_fast32_t iy = 0; iy < _number_of_cells.y(); ++iy) {
        const uint_fast32_t oiy =
            (offset_y * _number_of_cells.y() + iy) * number_of_modes;
        for (int_fast32_t iz = 0; iz < _number_of_cells.z(); ++iz) {
          const uint_fast32_t oiz =
              (offset_z * _number_of_cells.z() + iz) * number_of_modes;
          field[cell_index] = evaluate(oix, oiy, oiz);
          ++cell_index;
        }
      }
    }
  }
};

#endif // SAMPLEDINITIALTURBULENCE_HPP
