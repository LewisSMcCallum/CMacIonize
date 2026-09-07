#include "Assert.hpp"
#include "TransparentSphericalCavity.hpp"
#include <cmath>

int main() {
  const CoordinateVector<> centre(3., -2., 4.);
  for (double offset : {0., 0.6, 0.999999}) {
    PhotonPacket photon;
    const double half_chord = std::sqrt(1.-offset*offset);
    photon.set_position(centre+CoordinateVector<>(-half_chord, offset, 0.));
    photon.set_direction(CoordinateVector<>(1., 0., 0.));
    photon.set_target_optical_depth(0.7);
    photon.set_weight(0.3);
    photon.set_distance_travelled(2.);
    assert_condition(TransparentSphericalCavity::cross(photon, centre, 1., -1.));
    assert_condition(std::abs((photon.get_position()-centre)[0]-half_chord)<1.e-10);
    assert_condition(std::abs(photon.get_distance_travelled()-2.-2.*half_chord)<1.e-10);
    assert_condition(photon.get_target_optical_depth()==0.7);
    assert_condition(photon.get_weight()==0.3);
    assert_condition((photon.get_position()-centre).norm()>=1.);
  }
  PhotonPacket limited;
  limited.set_position(centre+CoordinateVector<>(-1., 0., 0.));
  limited.set_direction(CoordinateVector<>(1., 0., 0.));
  limited.set_target_optical_depth(0.7);
  limited.set_distance_travelled(2.);
  assert_condition(!TransparentSphericalCavity::cross(limited, centre, 1., 2.5));
  assert_condition(limited.get_distance_travelled()==2.5);
  assert_condition(limited.get_target_optical_depth()==0.7);
  return 0;
}
