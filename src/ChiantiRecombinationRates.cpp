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
 * @file ChiantiRecombinationRates.cpp
 *
 * @brief RecombinationRates implementation with Chianti's recombination rates:
 * implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "ChiantiRecombinationRates.hpp"
#include "Error.hpp"
#include "ChiantiRecombinationRatesDataLocation.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>


/**
 * @brief Constructor.
 *
 * Reads in the data file.
 */
ChiantiRecombinationRates::ChiantiRecombinationRates() {
  std::cout << "Setting up recombination rates from CHIANTI" << std::endl;
  std::cout << "HAS " << NUMBER_OF_IONNAMES << " different ions." << std::endl;

    _recomb_rates.resize(NUMBER_OF_IONNAMES);
    _temperatures.resize(350);


    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {

      std::cout << "Setting up data from file - " << get_ion_recombination_filename(i) << std::endl;

      _recomb_rates[i].resize(350);


      std::stringstream filenamestream;
      filenamestream << CHIANTIRECOMBDATALOCATION
                     << get_ion_recombination_filename(i);




      std::ifstream drfile(filenamestream.str());

      // skip the first two lines
      std::string line;
      std::getline(drfile, line);
      std::getline(drfile, line);
      // now parse the remaining lines
      for (uint_fast32_t j = 0; j < 350; ++j) {
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
        350. /
        (std::log(_temperatures[350 - 1]) -
         _min_logT);



}

/**
 * @brief Get the Chianti recombination rate.
 *
 * This code is identical to the code in Chianti's rrfit, except for the indices,
 * as C++ starts counting from zero instead of one.
 *
 * @param iz Atomic number.
 * @param in Number of electrons.
 * @param T Temperature (in K).
 * @return Recombination rate (in cm^3s^-1).
 */
double ChiantiRecombinationRates::get_recombination_rate_chianti(const int ion, const double temperature) const {




      // we need to find the index of the lower limit and the upper limit of
      // the temperature interval that contains the given temperature
      uint_fast32_t ilow, ihigh;

      // first handle the special cases
      if (temperature < _temperatures[0]) {
        // we will linearly extrapolate
        ilow = 0;
        ihigh = 1;
      } else if (temperature >=
                 _temperatures[350 - 1]) {
        ilow = 350 - 2;
        ihigh = 350 - 1;
      } else {

        // normal case
        // first, get a reasonable first guess for ilow and ihigh
        ilow = static_cast< uint_fast32_t >((std::log(temperature) - _min_logT) *
                                            _inverse_avg_dlogT);
        if (temperature < _temperatures[ilow]) {
          ihigh = ilow;
          ilow = 0;
        } else {
          ihigh = 350 - 1;
        }


        cmac_assert(temperature < _temperatures[ihigh]);
        cmac_assert(temperature >= _temperatures[ilow]);

        // now search for the actual indices using bisection
        while ((ihigh - ilow) != 1) {

          uint_fast32_t imid = (ilow + ihigh) >> 1;
          if (imid > 349 || imid<0) {
            ilow = 0;
            ihigh = 349;
            imid = 175;
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

          double rate =
          (1. - fac) * _recomb_rates[ion][ilow] + fac * _recomb_rates[ion][ihigh];


      if (temperature >=
                 _temperatures[350 - 1]) {

      rate = _recomb_rates[ion][350 - 1];

      }






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
double ChiantiRecombinationRates::get_recombination_rate(
    const int_fast32_t ion, const double temperature) const {

  double rate = 0.;

  switch (ion) {

  case ION_H_n: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }

#ifdef HAS_HELIUM
  case ION_He_n: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }

  case ION_He_p1: {

    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }

#endif

#ifdef HAS_CARBON
  case ION_C_p1: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }
  case ION_C_p2: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }
#endif

#ifdef HAS_NITROGEN
  case ION_N_n: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }
  case ION_N_p1: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }
  case ION_N_p2: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }
#endif

#ifdef HAS_OXYGEN
  case ION_O_n: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }
  case ION_O_p1: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }
  case ION_O_p2: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }
  case ION_O_p3: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }
#endif

#ifdef HAS_NEON
  case ION_Ne_n:
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  case ION_Ne_p1: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }
  case ION_Ne_p2: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }
  case ION_Ne_p3: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }
#endif

#ifdef HAS_SULPHUR
  case ION_S_p1: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }
  case ION_S_p2: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }
  case ION_S_p3: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }
#endif

#ifdef HAS_MAGNESIUM
  case ION_Mg_p1: {
    rate = get_recombination_rate_chianti(ion,temperature);
    break;
  }
#endif

  default:
    cmac_error("Unknown ion: %" PRIiFAST32, ion);
  }
  // convert cm^3s^-1 to m^3s^-1
  rate *= 1.e-6;


  // some rates become negative for large T: make sure we don't use these
  // values
  return std::max(0., rate);
}
