#ifndef HYDROGENEXCITATIONCOOLING_HPP
#define HYDROGENEXCITATIONCOOLING_HPP

#include <cmath>

/// Electron-impact excitation of neutral H (predominantly Ly-alpha).
/// Katz, Weinberg & Hernquist (1996), Table 1, following Cen (1992).
/// T in K, both densities in m^-3, result in J m^-3 s^-1.
/// For explicit process sums only: do not add to a total cooling table.
inline double hydrogen_excitation_cooling(const double T, const double ne,
                                          const double nHI) {
  if (T <= 0. || ne <= 0. || nHI <= 0.) return 0.;
  // 7.5e-19 erg cm^3 s^-1 = 7.5e-32 J m^3 s^-1.
  return 7.5e-32 * std::exp(-118348. / T) / (1. + std::sqrt(T / 1.e5)) *
         ne * nHI;
}

#endif
