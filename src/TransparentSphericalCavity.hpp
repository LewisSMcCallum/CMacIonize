/** @file TransparentSphericalCavity.hpp
 *  @brief Free flight across the excluded cubed-sphere centre (no opacity).
 */
#ifndef TRANSPARENTSPHERICALCAVITY_HPP
#define TRANSPARENTSPHERICALCAVITY_HPP
#include "PhotonPacket.hpp"
#include "Error.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace TransparentSphericalCavity {
/** Entry is on, or a roundoff displacement inside, the inner boundary.
 *  Preserve optical depth, weight and spectrum; count the vacuum path length.
 *  Return false when the existing photon-distance limit terminates the flight.
 */
inline bool cross(PhotonPacket &photon, const CoordinateVector<> &centre,
                  const double radius, const double maximum_distance) {
  const CoordinateVector<> x = photon.get_position() - centre;
  const CoordinateVector<> d = photon.get_direction();
  const double a = d.norm2();
  const double b = x[0]*d[0] + x[1]*d[1] + x[2]*d[2];
  const double discriminant = b*b - a*(x.norm2()-radius*radius);
  if (!(a > 0.) || !std::isfinite(discriminant)) {
    cmac_error("Invalid photon direction/position at spherical cavity.");
  }
  // Roundoff can put a grazing entry infinitesimally outside the sphere.
  const double tolerance = 1024.*std::numeric_limits<double>::epsilon()*
                           std::max(radius*radius*a, b*b);
  if (discriminant < -tolerance) {
    cmac_error("Photon entering spherical cavity misses its inner sphere.");
  }
  const double distance = std::max(0., (-b+std::sqrt(std::max(0.,discriminant)))/a);
  const double path = distance*std::sqrt(a);
  if (maximum_distance > 0. &&
      photon.get_distance_travelled()+path >= maximum_distance) {
    const double remaining = std::max(0., maximum_distance-photon.get_distance_travelled());
    photon.set_position(photon.get_position()+remaining/std::sqrt(a)*d);
    photon.set_distance_travelled(maximum_distance);
    return false;
  }
  // Only nudge by local floating-point precision, not the outer box size.
  const double epsilon = 256.*std::numeric_limits<double>::epsilon()*radius;
  photon.set_position(photon.get_position()+(distance+epsilon/std::sqrt(a))*d);
  photon.add_distance_travelled(path);
  return true;
}
}
#endif
