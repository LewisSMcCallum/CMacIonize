/*******************************************************************************
 * This file is part of CMacIonize
 ******************************************************************************/

/**
 * @file testSphericalDensitySubGridEntry.cpp
 *
 * @brief Regression test for near-tangent cubed-sphere subgrid entry.
 */

#include "Assert.hpp"
#include "SphericalDensitySubGrid.hpp"

#include <cmath>
#include <limits>
#include <vector>

int main(int argc, char **argv) {

  // Reproduce the radial block and angular boundary involved in a photon that
  // previously bounced forever between the two neighbouring v subgrids.
  const double first_outer_edge = 2.11735e18;
  const double outer_edge = 3.851651256129982e19;
  const double ratio = std::pow(outer_edge / first_outer_edge, 1. / 516.);
  std::vector< double > all_edges(518);
  all_edges[0] = 3.085677581491367e15;
  all_edges[1] = first_outer_edge;
  for (size_t i = 2; i < all_edges.size(); ++i) {
    all_edges[i] = all_edges[i - 1] * ratio;
  }
  std::vector< double > radial_edges(all_edges.begin() + 43 * 11,
                                     all_edges.begin() + 44 * 11 + 1);

  const int_fast32_t angular_subgrids = 18;
  const int_fast32_t iu = 6;
  const int_fast32_t iv = 3;
  const double u0 = -1. + 2. * iu / angular_subgrids;
  const double u1 = -1. + 2. * (iu + 1) / angular_subgrids;
  const double v0 = -1. + 2. * iv / angular_subgrids;
  const double v1 = -1. + 2. * (iv + 1) / angular_subgrids;
  SphericalDensitySubGrid grid(
      CoordinateVector<>(0.), 5, u0, u1, v0, v1, radial_edges,
      CoordinateVector< int_fast32_t >(20, 20, 11));

  PhotonPacket photon;
  photon.set_position(CoordinateVector<>(
      -1.6457336352540963e19, -6.9692830372866918e18,
      -2.4686004528811438e19));
  photon.set_direction(
      CoordinateVector<>(-0.554621, -0.0057935, -0.832083));
  photon.set_target_optical_depth(
      0.25 * std::numeric_limits< double >::max());
  photon.set_weight(0.);
  photon.set_energy(0.);
  photon.set_dust_opacity(0.);
  for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
    photon.set_photoionization_cross_section(ion, 0.);
  }

  const int_fast32_t output =
      grid.interact(photon, TRAVELDIRECTION_FACE_Y_N, -1.);

  // It entered through Y_N, so immediately returning through Y_N is the
  // floating-point bounce that this regression test guards against.
  assert_condition(output != TRAVELDIRECTION_FACE_Y_N);
  assert_condition(photon.get_distance_travelled() > 0.);

  return 0;
}
