/*******************************************************************************
 * This file is part of CMacIonize
 ******************************************************************************/

/**
 * @file SphericalDensitySubGrid.hpp
 *
 * @brief Cubed-sphere geometry for the task-based ionization solver.
 */
#ifndef SPHERICALDENSITYSUBGRID_HPP
#define SPHERICALDENSITYSUBGRID_HPP

#include "DensitySubGrid.hpp"

#include <array>
#include <limits>
#include <vector>

/**
 * @brief A radial block on one face of a gnomonic cubed sphere.
 *
 * The three logical cell coordinates are (u,v,r).  Keeping the usual six
 * face ports means the existing task queues, buffers and source-copy scheme
 * can be reused unchanged.
 */
class SphericalDensitySubGrid : public DensitySubGrid {
private:
  CoordinateVector<> _centre;
  int_fast8_t _face;
  double _u0, _u1, _v0, _v1;
  std::vector< double > _radial_edges;
  int_fast32_t _neighbour_inputs[TRAVELDIRECTION_NUMBER];

  static inline const double *dummy_box() {
    static const double box[6] = {0., 0., 0., 1., 1., 1.};
    return box;
  }

  static inline CoordinateVector<> normal(const int_fast8_t f) {
    static const double x[6][3] = {{1, 0, 0},  {-1, 0, 0}, {0, 1, 0},
                                    {0, -1, 0}, {0, 0, 1},  {0, 0, -1}};
    return CoordinateVector<>(x[f][0], x[f][1], x[f][2]);
  }

  static inline CoordinateVector<> axis_u(const int_fast8_t f) {
    static const double x[6][3] = {{0, 1, 0}, {0, -1, 0}, {-1, 0, 0},
                                    {1, 0, 0}, {0, 1, 0},  {0, 1, 0}};
    return CoordinateVector<>(x[f][0], x[f][1], x[f][2]);
  }

  static inline CoordinateVector<> axis_v(const int_fast8_t f) {
    static const double x[6][3] = {{0, 0, 1}, {0, 0, 1}, {0, 0, 1},
                                    {0, 0, 1}, {-1, 0, 0}, {1, 0, 0}};
    return CoordinateVector<>(x[f][0], x[f][1], x[f][2]);
  }

  static inline double dot(const CoordinateVector<> &a,
                           const CoordinateVector<> &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }

  inline void local_coordinates(const CoordinateVector<> &position, double &u,
                                double &v, double &r) const {
    const CoordinateVector<> x = position - _centre;
    r = x.norm();
    const double denominator = dot(normal(_face), x);
    u = dot(axis_u(_face), x) / denominator;
    v = dot(axis_v(_face), x) / denominator;
  }

  inline CoordinateVector<> direction(const double u, const double v) const {
    CoordinateVector<> q =
        normal(_face) + u * axis_u(_face) + v * axis_v(_face);
    return q / q.norm();
  }

  inline void locate(const CoordinateVector<> &position,
                     CoordinateVector< int_fast32_t > &i) const {
    double u, v, r;
    local_coordinates(position, u, v, r);
    i[0] = std::min(
        _number_of_cells[0] - 1,
        std::max< int_fast32_t >(
            0, std::floor((u - _u0) * _number_of_cells[0] / (_u1 - _u0))));
    i[1] = std::min(
        _number_of_cells[1] - 1,
        std::max< int_fast32_t >(
            0, std::floor((v - _v0) * _number_of_cells[1] / (_v1 - _v0))));
    const std::vector< double >::const_iterator it =
        std::upper_bound(_radial_edges.begin(), _radial_edges.end(), r);
    i[2] = std::min(
        _number_of_cells[2] - 1,
        std::max< int_fast32_t >(0, it - _radial_edges.begin() - 1));
  }

  inline double angular_solid_angle(const double u0, const double u1,
                                    const double v0, const double v1) const {
    // Exact solid angle of a gnomonic quadrilateral, split into two triangles.
    const CoordinateVector<> a = direction(u0, v0);
    const CoordinateVector<> b = direction(u1, v0);
    const CoordinateVector<> c = direction(u1, v1);
    const CoordinateVector<> d = direction(u0, v1);
    const auto triangle = [](const CoordinateVector<> &p,
                             const CoordinateVector<> &q,
                             const CoordinateVector<> &r) {
      const double determinant =
          p[0] * (q[1] * r[2] - q[2] * r[1]) -
          p[1] * (q[0] * r[2] - q[2] * r[0]) +
          p[2] * (q[0] * r[1] - q[1] * r[0]);
      const double denominator =
          1. + dot(p, q) + dot(q, r) + dot(r, p);
      return 2. * std::atan2(std::abs(determinant), denominator);
    };
    return triangle(a, b, c) + triangle(a, c, d);
  }

public:
  SphericalDensitySubGrid(
      const CoordinateVector<> centre, const int_fast8_t face, const double u0,
      const double u1, const double v0, const double v1,
      const std::vector< double > &radial_edges,
      const CoordinateVector< int_fast32_t > number_of_cells)
      : DensitySubGrid(dummy_box(), number_of_cells),
        _centre(centre), _face(face), _u0(u0), _u1(u1), _v0(v0), _v1(v1),
        _radial_edges(radial_edges) {
    for (int_fast32_t i = 0; i < TRAVELDIRECTION_NUMBER; ++i) {
      _neighbour_inputs[i] =
          TravelDirections::output_to_input_direction(i);
    }
  }

  SphericalDensitySubGrid(const SphericalDensitySubGrid &other)
      : DensitySubGrid(other), _centre(other._centre), _face(other._face),
        _u0(other._u0), _u1(other._u1), _v0(other._v0), _v1(other._v1),
        _radial_edges(other._radial_edges) {
    std::copy(other._neighbour_inputs,
              other._neighbour_inputs + TRAVELDIRECTION_NUMBER,
              _neighbour_inputs);
  }

  virtual DensitySubGrid *clone() const {
    return new SphericalDensitySubGrid(*this);
  }

  inline void set_neighbour_input_direction(const int_fast32_t output,
                                            const int_fast32_t input) {
    _neighbour_inputs[output] = input;
  }

  virtual int_fast32_t
  get_neighbour_input_direction(const int_fast32_t output) const {
    return _neighbour_inputs[output];
  }

  virtual bool is_in_box(const CoordinateVector<> position) const {
    double u, v, r;
    local_coordinates(position, u, v, r);
    return u >= _u0 && u <= _u1 && v >= _v0 && v <= _v1 &&
           r >= _radial_edges.front() && r <= _radial_edges.back();
  }

  virtual CoordinateVector<>
  get_cell_midpoint(const uint_fast32_t index) const {
    CoordinateVector< int_fast32_t > i;
    get_three_index(index, i);
    const double du = (_u1 - _u0) / _number_of_cells[0];
    const double dv = (_v1 - _v0) / _number_of_cells[1];
    const double u = _u0 + (i[0] + 0.5) * du;
    const double v = _v0 + (i[1] + 0.5) * dv;
    const double r =
        0.5 * (_radial_edges[i[2]] + _radial_edges[i[2] + 1]);
    return _centre + r * direction(u, v);
  }

  virtual double get_cell_volume(const uint_fast32_t index) const {
    CoordinateVector< int_fast32_t > i;
    get_three_index(index, i);
    const double du = (_u1 - _u0) / _number_of_cells[0];
    const double dv = (_v1 - _v0) / _number_of_cells[1];
    const double omega =
        angular_solid_angle(_u0 + i[0] * du, _u0 + (i[0] + 1) * du,
                            _v0 + i[1] * dv, _v0 + (i[1] + 1) * dv);
    return omega *
           (std::pow(_radial_edges[i[2] + 1], 3) -
            std::pow(_radial_edges[i[2]], 3)) /
           3.;
  }

  virtual iterator get_cell(const CoordinateVector<> position) {
    CoordinateVector< int_fast32_t > i;
    locate(position, i);
    return iterator(get_one_index(i), *this);
  }

  virtual int_fast32_t interact(PhotonPacket &photon,
                                const int_fast32_t input_direction,
                                const double maximum_distance) {
    const CoordinateVector<> d = photon.get_direction();
    CoordinateVector<> p = photon.get_position();
    const double scale = std::max(_radial_edges.back(), 1.);
    const double epsilon = 1.e-11 * scale;
    // A source may lie exactly on a u/v/r cell face.  Select the downstream
    // cell using the photon direction, just as we do at subgrid boundaries.
    p += epsilon * d;

    // Sources on a subgrid boundary are assigned to one of the two blocks.
    // A photon launched into the other block must be forwarded immediately;
    // clamping it to a cell on the wrong side creates a one-sided source.
    double start_u, start_v, start_r;
    local_coordinates(p, start_u, start_v, start_r);
    int_fast32_t initial_output = TRAVELDIRECTION_INSIDE;
    if (start_u < _u0) {
      initial_output = TRAVELDIRECTION_FACE_X_N;
    } else if (start_u > _u1) {
      initial_output = TRAVELDIRECTION_FACE_X_P;
    } else if (start_v < _v0) {
      initial_output = TRAVELDIRECTION_FACE_Y_N;
    } else if (start_v > _v1) {
      initial_output = TRAVELDIRECTION_FACE_Y_P;
    } else if (start_r < _radial_edges.front()) {
      initial_output = TRAVELDIRECTION_FACE_Z_N;
    } else if (start_r > _radial_edges.back()) {
      initial_output = TRAVELDIRECTION_FACE_Z_P;
    }
    if (initial_output != TRAVELDIRECTION_INSIDE) {
      photon.set_position(p);
      return initial_output;
    }

    CoordinateVector< int_fast32_t > cell;
    locate(p, cell);
    double tau_done = 0.;
    const double tau_target = photon.get_target_optical_depth();
    int_fast32_t output = TRAVELDIRECTION_INSIDE;

    while (is_inside(cell) && tau_done < tau_target) {
      const double du = (_u1 - _u0) / _number_of_cells[0];
      const double dv = (_v1 - _v0) / _number_of_cells[1];
      const double ub[2] = {_u0 + cell[0] * du,
                            _u0 + (cell[0] + 1) * du};
      const double vb[2] = {_v0 + cell[1] * dv,
                            _v0 + (cell[1] + 1) * dv};
      const CoordinateVector<> x = p - _centre;
      double best = std::numeric_limits< double >::max();
      int best_axis = -1;
      int best_sign = 0;

      for (int side = 0; side < 2; ++side) {
        const CoordinateVector<> plane =
            axis_u(_face) - ub[side] * normal(_face);
        const double denominator = dot(plane, d);
        if (denominator != 0.) {
          const double distance = -dot(plane, x) / denominator;
          if (distance > epsilon && distance < best) {
            best = distance;
            best_axis = 0;
            best_sign = side ? 1 : -1;
          }
        }
      }
      for (int side = 0; side < 2; ++side) {
        const CoordinateVector<> plane =
            axis_v(_face) - vb[side] * normal(_face);
        const double denominator = dot(plane, d);
        if (denominator != 0.) {
          const double distance = -dot(plane, x) / denominator;
          if (distance > epsilon && distance < best) {
            best = distance;
            best_axis = 1;
            best_sign = side ? 1 : -1;
          }
        }
      }

      // Intersections with the two spherical radial faces.
      for (int side = 0; side < 2; ++side) {
        const double radius = _radial_edges[cell[2] + side];
        const double b = dot(x, d);
        const double discriminant = b * b - dot(x, x) + radius * radius;
        if (discriminant >= 0.) {
          const double root = std::sqrt(discriminant);
          const double candidates[2] = {-b - root, -b + root};
          for (int k = 0; k < 2; ++k) {
            if (candidates[k] > epsilon && candidates[k] < best) {
              best = candidates[k];
              best_axis = 2;
              best_sign = side ? 1 : -1;
            }
          }
        }
      }

      if (!std::isfinite(best)) {
        cmac_error("Could not find a cubed-sphere cell exit.");
      }

      const int_fast32_t active_cell = get_one_index(cell);
      const double tau = get_optical_depth(active_cell, best, photon);
      double distance = best;
      tau_done += tau;
      if (tau_done >= tau_target && tau > 0.) {
        distance *= 1. - (tau_done - tau_target) / tau;
      }
      update_intensity_counters(active_cell, distance, photon);
      p += distance * d;

      if (tau_done >= tau_target ||
          (maximum_distance > 0. &&
           photon.get_distance_travelled() >= maximum_distance)) {
        output = TRAVELDIRECTION_INSIDE;
        break;
      }

      cell[best_axis] += best_sign;
      if (!is_inside(cell)) {
        if (best_axis == 0) {
          output = best_sign > 0 ? TRAVELDIRECTION_FACE_X_P
                                 : TRAVELDIRECTION_FACE_X_N;
        } else if (best_axis == 1) {
          output = best_sign > 0 ? TRAVELDIRECTION_FACE_Y_P
                                 : TRAVELDIRECTION_FACE_Y_N;
        } else {
          output = best_sign > 0 ? TRAVELDIRECTION_FACE_Z_P
                                 : TRAVELDIRECTION_FACE_Z_N;
        }
      } else {
        p += epsilon * d;
      }
    }

    photon.set_target_optical_depth(std::max(0., tau_target - tau_done));
    photon.set_position(p);
    return output;
  }
};

#endif // SPHERICALDENSITYSUBGRID_HPP
