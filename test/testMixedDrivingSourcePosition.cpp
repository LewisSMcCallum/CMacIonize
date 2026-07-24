/*******************************************************************************
 * This file is part of CMacIonize
 ******************************************************************************/

#include "Assert.hpp"
#include "MixedDrivingPhotonSourceDistribution.hpp"

int main(int argc, char **argv) {
  const Box<> box(CoordinateVector<>(-5., -4., -3.),
                  CoordinateVector<>(10., 8., 6.));
  const CoordinateVector< bool > periodicity(true, true, false);

  CoordinateVector<> position(6., -5., 2.);
  wrap_mixed_driving_source_position(position, box, periodicity);
  assert_condition(position == CoordinateVector<>(-4., 3., 2.));

  CoordinateVector<> unchanged(1., 2., -2.);
  wrap_mixed_driving_source_position(unchanged, box, periodicity);
  assert_condition(unchanged == CoordinateVector<>(1., 2., -2.));

  return 0;
}
