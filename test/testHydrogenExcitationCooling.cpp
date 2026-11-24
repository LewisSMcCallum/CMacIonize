#include "HydrogenExcitationCooling.hpp"
#include "Abundances.hpp"
#include "ChargeTransferRates.hpp"
#include "CollisionalRates.hpp"
#include "DeRijckeRadiativeCooling.hpp"
#include "IonizationVariables.hpp"
#include "LineCoolingData.hpp"
#include "TemperatureCalculator.hpp"
#include "VernerRecombinationRates.hpp"
#include <cmath>

int main() {
  // Published cgs fit at 10^4 K, for 1 electron and 1 H I per cm^3.
  const double loss = hydrogen_excitation_cooling(10000., 1.e6, 1.e6);
  const double expected = 4.12992229945e-24 * 0.1; // erg cm^-3 s^-1 to SI
  if (std::abs(loss / expected - 1.) > 1.e-10) return 1;
  if (hydrogen_excitation_cooling(10000., 0., 1.e6) != 0. ||
      hydrogen_excitation_cooling(10000., 1.e6, 0.) != 0. ||
      hydrogen_excitation_cooling(0., 1.e6, 1.e6) != 0.) return 2;
  if (std::abs(hydrogen_excitation_cooling(10000., 2.e6, 3.e6) /
               (6. * loss) - 1.) > 1.e-12) return 3;
  // Exercise the production balance twice on the same, partially ionized
  // pure-H cell. The tabulated total must not get a second H I contribution.
  const double n = 1.e6, T = 1.e4;
  Abundances abundances;
  IonizationVariables cell;
  cell.set_number_density(n);
#ifdef VARIABLE_ABUNDANCES
  cell.get_abundances().set_abundances(abundances);
#endif
  LineCoolingData lines;
  VernerRecombinationRates recombination;
  ChargeTransferRates charge;
  CollisionalRates collision;
  DeRijckeRadiativeCooling table;
  double j[NUMBER_OF_IONNAMES] = {}, h[NUMBER_OF_HEATINGTERMS] = {};
  j[ION_H_n] = 1.e-13;
  double h0 = 0., he0 = 0., gain = 0., total = 0.;
  auto balance = [&](const DeRijckeRadiativeCooling *curve) {
    TemperatureCalculator::compute_cooling_and_heating_balance(
        h0, he0, gain, total, T, cell, CoordinateVector<>(), j[ION_H_n],
        0., j, abundances, h, 0., 0., 0., lines, recombination, charge,
        collision, curve);
  };
  balance(nullptr);
  if (!(h0 > 0. && h0 < 1.)) return 4;
  const double ne = n * (1. - h0), c = 5.5 - std::log(T);
  const double freefree = 1.42e-40 * (1.1 + .34*std::exp(-c*c/3.)) * std::sqrt(T);
  const double recomb = 2.85e-40 * std::sqrt(T) *
      (5.914 - .5*std::log(T) + .01184*std::cbrt(T));
  const double sum = (freefree + recomb)*ne*ne + expected*(ne/1.e6)*(n*h0/1.e6);
  if (!(std::abs(total/sum - 1.) < 1.e-10)) return 5;
  balance(&table);
  if (!(std::abs(total/(table.get_cooling_rate(T)*n*n) - 1.) < 1.e-12)) return 6;
  return 0;
}
