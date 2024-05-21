/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file BenjaminRecombinationRates.cpp
 *
 * @brief RecombinationRates implementation with Benjamin's recombination rates:
 * implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "BenjaminRecombinationRates.hpp"
#include "Error.hpp"
#include "BenjaminRecombinationRatesDataLocation.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include "VernerRecombinationRates.hpp"

/**
 * @brief Constructor.
 *
 * Reads in the data file.
 */
BenjaminRecombinationRates::BenjaminRecombinationRates() {
  std::cout << "HAS " << NUMBER_OF_IONNAMES << " different ions." << std::endl;

    _recomb_rates.resize(NUMBER_OF_IONNAMES);
    _temperatures.resize(250);


    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {

      std::cout << "Setting up data from file - " << get_ion_recombination_filename(i) << std::endl;

      _recomb_rates[i].resize(250);


      std::stringstream filenamestream;
      filenamestream << BENJAMINRECOMBDATALOCATION
                     << get_ion_recombination_filename(i);




      std::ifstream drfile(filenamestream.str());

      // skip the first two lines
      std::string line;
      std::getline(drfile, line);
      std::getline(drfile, line);
      // now parse the remaining lines
      for (uint_fast32_t j = 0; j < 250; ++j) {
        std::getline(drfile, line);

        std::istringstream linestream(line);

        //overwriting temperatures 181 times.... not great but fine

        linestream >> _temperatures[j] >> _recomb_rates[i][j];
        // temperature to K
        _temperatures[j] = std::pow(10.,_temperatures[j]);
        // cm^3/s, leave like this and multiplier comes later

        //std::cout << _temperatures[j] << "kelvin with " << _recomb_rates[i][j] << std::endl;
      }

    }

    _min_logT = std::log(_temperatures[0]);
    _inverse_avg_dlogT =
        250. /
        (std::log(_temperatures[250 - 1]) -
         _min_logT);


  _verner = new VernerRecombinationRates();


}

/**
 * @brief Get the Benjamin recombination rate.
 *
 * This code is identical to the code in Benjamin's rrfit, except for the indices,
 * as C++ starts counting from zero instead of one.
 *
 * @param iz Atomic number.
 * @param in Number of electrons.
 * @param T Temperature (in K).
 * @return Recombination rate (in cm^3s^-1).
 */
double BenjaminRecombinationRates::get_recombination_rate_benjamin(const int ion, const double temperature) const {




      // we need to find the index of the lower limit and the upper limit of
      // the temperature interval that contains the given temperature
      uint_fast32_t ilow, ihigh;

      // first handle the special cases
      if (temperature < _temperatures[0]) {
        // we will linearly extrapolate
        ilow = 0;
        ihigh = 1;
      } else if (temperature >=
                 _temperatures[250 - 1]) {
        ilow = 250 - 2;
        ihigh = 250 - 1;
      } else {

        // normal case
        // first, get a reasonable first guess for ilow and ihigh
        ilow = static_cast< uint_fast32_t >((std::log(temperature) - _min_logT) *
                                            _inverse_avg_dlogT);
        if (temperature < _temperatures[ilow]) {
          ihigh = ilow;
          ilow = 0;
        } else {
          ihigh = 250 - 1;
        }


        cmac_assert(temperature < _temperatures[ihigh]);
        cmac_assert(temperature >= _temperatures[ilow]);

        // now search for the actual indices using bisection
        while ((ihigh - ilow) != 1) {

          uint_fast32_t imid = (ilow + ihigh) >> 1;
          if (imid > 249 || imid<0) {
            ilow = 0;
            ihigh = 249;
            imid = 125;
          }
          if (temperature >= _temperatures[imid]) {
            ilow = imid;
          } else {
            ihigh = imid;
          }
        }
        cmac_assert(ilow < ihigh);
        cmac_assert((ihigh - ilow) == 1);
        cmac_assert(temperature < _temperatures[ihigh]);
        cmac_assert(temperature >= _temperatures[ilow]);
      }


      // we now have the appropriate interval for linear inter/extrapolation
      const double fac = (temperature - _temperatures[ilow]) /
                         (_temperatures[ihigh] - _temperatures[ilow]);

          const double rate =
          (1. - fac) * _recomb_rates[ion][ilow] + fac * _recomb_rates[ion][ihigh];






      return std::max(rate, 0.);


}

/**
 * @brief Get the recombination rate of the given ion at the given
 * temperature.
 *
 * @param ion IonName for a valid ion.
 * @param temperature Temperature (in K).
 * @return Recombination rate (in m^3s^-1).
 */
double BenjaminRecombinationRates::get_recombination_rate(
    const int_fast32_t ion, const double temperature) const {

  double rate = 0.;

  switch (ion) {

  case ION_H_n: {
    // Benjamin & Ferland (1996) formula (4) with values from Table 1 (HI).
    const double T1 = temperature / 3.148;
    const double T2 = temperature / 7.036e5;
    rate = 7.982e-11 / (std::sqrt(T1) * std::pow(1. + std::sqrt(T1), 0.252) *
                        std::pow(1. + std::sqrt(T2), 1.748));
    break;
  }

#ifdef HAS_HELIUM
  case ION_He_n: {
    // Benjamin & Ferland (1996) formula (4) with values from Table 1 (HeIa).
    // Note that we use the second version, which is valid in the range
    // [3 K, 10^10 K].
    const double T1 = temperature / 4.266e-2;
    const double T2 = temperature / 4.677e6;
    rate = 9.356e-10 / (std::sqrt(T1) * std::pow(1. + std::sqrt(T1), 0.2108) *
                        std::pow(1. + std::sqrt(T2), 1.7892));
    break;
  }

  case ION_He_p1: {
    // Benjamin & Ferland (1996) formula (4) with values from Table 1 (HeII).
    // Note that we use the first version, which is only valid in the range
    // [3 K, 10^9 K].
    const double T1 = temperature / 0.937;
    const double T2 = temperature / 2.774e6;
    rate = 1.891e-10 / (std::sqrt(T1) * std::pow(1. + std::sqrt(T1), 0.2476) *
                        std::pow(1. + std::sqrt(T2), 1.7524));
    break;
  }

#endif

#ifdef HAS_CARBON
  case ION_C_p1: {
    // Nussbaumer & Storey (1983) formula (19) with values from Table 1 (C2+).
    // valid in the range [1,000 K; 60,000 K]
    const double T4 = temperature * 1.e-4;
    const double T4_inv = 1. / T4;
    if (temperature > 60000) {
      rate = get_recombination_rate_benjamin(ion,temperature);
    } else {
    rate = get_recombination_rate_benjamin(ion,temperature) +
           1.e-12 *
               (1.8267 * T4_inv + 4.1012 + 4.8443 * T4 + 0.2261 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.5960 * T4_inv);
    }
    break;
  }
  case ION_C_p2: {
    // Nussbaumer & Storey (1983) formula (19) with values from Table 1 (C3+).
    // valid in the range [1,000 K; 60,000 K]
    const double T4 = temperature * 1.e-4;
    const double T4_inv = 1. / T4;
    if (temperature > 60000) {
      rate = get_recombination_rate_benjamin(ion,temperature);
    } else {
    rate = get_recombination_rate_benjamin(ion,temperature) +
           1.e-12 *
               (2.3196 * T4_inv + 10.7328 + 6.8830 * T4 - 0.1824 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.4101 * T4_inv);
      }
    break;
  }
#endif

#ifdef HAS_NITROGEN
  case ION_N_n: {
    // Nussbaumer & Storey (1983) formula (19) with values from Table 1 (N+).
    // valid in the range [1,000 K; 60,000 K]
    const double T4 = temperature * 1.e-4;
    if (temperature > 60000) {
      rate = get_recombination_rate_benjamin(ion,temperature);
    } else {
    rate = get_recombination_rate_benjamin(ion,temperature) +
           1.e-12 * (0.6310 + 0.1990 * T4 - 0.0197 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.4398 / T4);
      }
    break;
  }
  case ION_N_p1: {
    // Nussbaumer & Storey (1983) formula (19) with values from Table 1 (N2+).
    // valid in the range [1,000 K; 60,000 K]
    const double T4 = temperature * 1.e-4;
    const double T4_inv = 1. / T4;
    if (temperature > 60000) {
      rate = get_recombination_rate_benjamin(ion,temperature);
    } else {
    rate = get_recombination_rate_benjamin(ion,temperature) +
           1.e-12 *
               (0.0320 * T4_inv - 0.6624 + 4.3191 * T4 + 0.0003 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.5946 * T4_inv);
    }
    break;
  }
  case ION_N_p2: {
    // Nussbaumer & Storey (1983) formula (19) with values from Table 1 (N3+).
    // valid in the range [1,000 K; 60,000 K]
    const double T4 = temperature * 1.e-4;
    const double T4_inv = 1. / T4;
    if (temperature > 60000) {
      rate = get_recombination_rate_benjamin(ion,temperature);
    } else {
    rate = get_recombination_rate_benjamin(ion,temperature) +
           1.e-12 *
               (-0.8806 * T4_inv + 11.2406 + 30.7066 * T4 - 1.1721 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.6127 * T4_inv);
    }
    break;
  }
#endif

#ifdef HAS_OXYGEN
  case ION_O_n: {
    // Nussbaumer & Storey (1983) formula (19) with values from Table 1 (O+).
    // version 2, valid in the range [1,000 K; 20,000 K]
    const double T4 = temperature * 1.e-4;
    const double T4_inv = 1. / T4;
    if (temperature > 20000){
      rate = get_recombination_rate_benjamin(ion,temperature);
    } else {
    rate = get_recombination_rate_benjamin(ion,temperature) +
           1.e-12 *
               (-0.0001 * T4_inv + 0.0001 + 0.0956 * T4 + 0.0193 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.4106 * T4_inv);
      }
    break;
  }
  case ION_O_p1: {
    // Nussbaumer & Storey (1983) formula (19) with values from Table 1 (O2+).
    // valid in the range [1,000 K; 60,000 K]
    const double T4 = temperature * 1.e-4;
    const double T4_inv = 1. / T4;
    if (temperature > 60000) {
      rate = get_recombination_rate_benjamin(ion,temperature);
    } else {
    rate = get_recombination_rate_benjamin(ion,temperature) +
           1.e-12 *
               (-0.0036 * T4_inv + 0.7519 + 1.5252 * T4 - 0.0838 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.2769 * T4_inv);
      }
    break;
  }
  case ION_O_p2: {
    // Nussbaumer & Storey (1983) formula (19) with values from Table 1 (O3+).
    // valid in the range [1,000 K; 60,000 K]
    const double T4 = temperature * 1.e-4;
    const double T4_inv = 1. / T4;
    if (temperature > 60000) {
      rate = get_recombination_rate_benjamin(ion,temperature);
    } else {
    rate = get_recombination_rate_benjamin(ion,temperature) +
           1.e-12 *
               (0.0 * T4_inv + 21.879 + 16.273 * T4 - 0.702 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-1.1899 * T4_inv);
        }
    break;
  }
  case ION_O_p3: {
    const double T4 = temperature * 1.e-4;
    const double T4_inv = 1. / T4;
    if (temperature < 20000) {
      rate = get_recombination_rate_benjamin(ion,temperature) +
             1.e-12 *
                 (-0.3648 * T4_inv + 7.2698 + 17.2187 * T4 + 9.8335 * T4 * T4) *
                 std::pow(T4, -1.5) * std::exp(0.0166 * T4_inv);
      break;

    } else if (temperature < 60000) {
      rate = get_recombination_rate_benjamin(ion,temperature) +
             1.e-12 *
                 (0-2.5053 * T4_inv + 3.4903 + 67.4128 * T4 - 3.445 * T4 * T4) *
                 std::pow(T4, -1.5) * std::exp(-0.8501 * T4_inv);
      break;
    } else {
      rate = get_recombination_rate_benjamin(ion,temperature);
      break;
    }
  }
#endif

#ifdef HAS_NEON
  case ION_Ne_n:
    // According to Nussbaumer & Storey (1987), dielectronic recombination is
    // negligible for this ion
    rate = get_recombination_rate_benjamin(ion,temperature);
    break;
  case ION_Ne_p1: {
    // Nussbaumer & Storey (1987) formula (7) with values from Table II(b)
    // (Total)
    // valid in the range [1,000 K; 60,000 K]
    // NOTE the sign difference in the first term w.r.t. Kenny's code.
    const double T4 = temperature * 1.e-4;
    const double T4_inv = 1. / T4;
    if (temperature > 60000) {
      rate = get_recombination_rate_benjamin(ion,temperature);
    } else {
    rate = get_recombination_rate_benjamin(ion,temperature) +
           1.e-12 *
               (0.0129 * T4_inv - 0.1779 + 0.9353 * T4 - 0.0682 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.4516 * T4_inv);
        }
    break;
  }
  case ION_Ne_p2: {
    // Nussbaumer & Storey (1987) formula (7) with values from Table III(b)
    // (Total)
    // valid in the range [1,000 K; 60,000 K]
    const double T4 = temperature * 1.e-4;
    const double T4_inv = 1. / T4;
    if (temperature > 60000) {
      rate = get_recombination_rate_benjamin(ion,temperature);
    } else {
    rate = get_recombination_rate_benjamin(ion,temperature) +
           1.e-12 *
               (3.6781 * T4_inv + 14.1481 + 17.1175 * T4 - 0.5017 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.2313 * T4_inv);
        }
    break;
  }
  case ION_Ne_p3: {
    rate = get_recombination_rate_benjamin(ion,temperature);
    break;
  }
#endif

#ifdef HAS_SULPHUR
  case ION_S_p1: {
    // Mazzotta et al. (1998), equation 7, validity range not specified
    // const double T_in_eV = temperature / 1.16045221e4;
    // if (temperature > 60000) {
    //   rate = get_recombination_rate_benjamin(ion,temperature);

    // } else {
    // rate = get_recombination_rate_benjamin(ion,temperature) +
    //        1.37e-9 * std::exp(-14.95 / T_in_eV) * std::pow(T_in_eV, -1.5);
    // }
    // Kaur et al. 2018, need to double check ranges
    const double T_inv = 1. / temperature;
    if (temperature > 9e7) {
      rate = get_recombination_rate_benjamin(ion,temperature);
    } else {
    rate = get_recombination_rate_benjamin(ion,temperature) +
           (3.04e-7 * std::exp(-5.016e1 * T_inv) +
            4.393e-7 * std::exp(-3.266e2 * T_inv) +
            1.609e-6 * std::exp(-3.102e3 * T_inv) +
            4.98e-6 * std::exp(-1.21e4 * T_inv) +
            3.457e-5 * std::exp(-4.969e4 * T_inv) +
            8.617e-3 * std::exp(-2.01e5 * T_inv) +
            9.284e-4 * std::exp(-2.575e5 * T_inv)) *
               std::pow(temperature, -1.5);
    }
    break;
  }
  case ION_S_p2: {
    // Mazzotta et al. (1998), equation 7, validity range not specified
    // const double T_in_eV = temperature / 1.16045221e4;
    // const double T_in_eV_inv = 1. / T_in_eV;
    // if (temperature > 60000) {
    //   rate = get_recombination_rate_benjamin(ion,temperature);
    // } else {
    // rate = get_recombination_rate_benjamin(ion,temperature) +
    //        (8.0729e-9 * std::exp(-17.56 * T_in_eV_inv) +
    //         1.1012e-10 * std::exp(-7.07 * T_in_eV_inv)) *
    //            std::pow(T_in_eV, -1.5);
    // }
    // Abdel-Naby et al. 2012, S3+
    const double T_inv = 1. / temperature;
    if (temperature > 9e7) {
      rate = get_recombination_rate_benjamin(ion,temperature);
    } else {
    rate = get_recombination_rate_benjamin(ion,temperature) +
           (5.817e-7 * std::exp(-362.8 * T_inv) +
            1.391e-6 * std::exp(-1058. * T_inv) +
            1.123e-5 * std::exp(-7160. * T_inv) +
            1.521e-4 * std::exp(-3.26e4 * T_inv) +
            1.875e-3 * std::exp(-1.235e5 * T_inv) +
            2.097e-2 * std::exp(-2.07e5 * T_inv)) *
               std::pow(temperature, -1.5);
    }
    break;
  }
  case ION_S_p3: {
    // Abdel-Naby et al. (2012), equation 3, validity range [90 K; 9x10^7 K].
    // const double T_inv = 1. / temperature;
    // if (temperature > 9e7) {
    //   rate = get_recombination_rate_benjamin(ion,temperature);
    // } else {
    // rate = get_recombination_rate_benjamin(ion,temperature) +
    //        (5.817e-7 * std::exp(-362.8 * T_inv) +
    //         1.391e-6 * std::exp(-1058. * T_inv) +
    //         1.123e-5 * std::exp(-7160. * T_inv) +
    //         1.521e-4 * std::exp(-3.26e4 * T_inv) +
    //         1.875e-3 * std::exp(-1.235e5 * T_inv) +
    //         2.097e-2 * std::exp(-2.07e5 * T_inv)) *
    //            std::pow(temperature, -1.5);
    // }
    // Altun et al. 2007
    const double T_inv = 1. / temperature;
    if (temperature > 9e7) {
      rate = get_recombination_rate_benjamin(ion,temperature);
    } else {
    rate = get_recombination_rate_benjamin(ion,temperature) +
           (9.571e-6 * std::exp(-1.18e3 * T_inv) +
            6.268e-5 * std::exp(-6.443e3 * T_inv) +
            3.807e-4 * std::exp(-2.264e4 * T_inv) +
            1.874e-2 * std::exp(-1.53e5 * T_inv) +
            5.526e-3 * std::exp(-3.564e5 * T_inv)) *
               std::pow(temperature, -1.5);
    }
    break;
  }
#endif

  default:
    cmac_error("Unknown ion: %" PRIiFAST32, ion);
  }
  // convert cm^3s^-1 to m^3s^-1
  rate *= 1.e-6;

  if (temperature < 20000.) {
    rate = _verner->get_recombination_rate(ion,temperature);
  }

  // some rates become negative for large T: make sure we don't use these
  // values
  return std::max(0., rate);
}
