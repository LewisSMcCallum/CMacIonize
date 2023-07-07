#ifndef COLLISIONALRATES_HPP
#define COLLISIONALRATES_HPP

#include "Error.hpp"
#include "ElementNames.hpp"

#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <iostream>

/**
 * @brief CollisionalRates from Bob Benjamin
 *
 */
class CollisionalRates {
private:
  /*! @brief Temperature values (in K). */
  double _temperatures[250];

  /*! @brief Collisional rate values (in m^-3 s^-1). */
  double _collisional_rates[NUMBER_OF_IONNAMES][250];

  /*! @brief Logarithm of the minimum tabulated temperature
   *  (in log(T / K)). */
  double _min_logT;

  /*! @brief Inverse of the average logarithmic distance between tabulated
   *  temperature values (in 1 / log(T / K)). */
  double _inverse_avg_dlogT;

public:
  CollisionalRates();

  /**
   * @brief Get the cooling rate for the given temperature.
   *
   * @param temperature Temperature (in K).
   * @return Cooling rate (in J m^3 s^-1).
   */
  inline double get_collisional_rate(const int ion, const double temperature) const {


    double single_ion_rates[250];

    for (int_fast32_t i = 0; i < 250; ++i) {

      single_ion_rates[i] = _collisional_rates[ion][i];

    }


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
    const double col_rate =
        (1. - fac) * single_ion_rates[ilow] + fac * single_ion_rates[ihigh];






    return std::max(col_rate, 0.);



  }

  /**
   * @brief Get the lowest tabulated temperature value.
   *
   * @return Lowest tabulated temperature value (in K).
   */
  double get_minimum_temperature() const { return _temperatures[0]; }

  double get_maximum_temperature() const { return _temperatures[250 - 1]; }
};

#endif // DERIJCKERADIATIVECOOLING_HPP
