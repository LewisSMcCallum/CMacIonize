#ifndef COLLISIONALRATES_HPP
#define COLLISIONALRATES_HPP

#include "Error.hpp"
#include "ElementNames.hpp"

#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <iostream>
#include <vector>

/**
 * @brief CollisionalRates from Bob Benjamin
 *
 */
class CollisionalRates {
private:
  /*! @brief Temperature values (in K). */
  std::vector<double> _temperatures;

  /*! @brief Collisional rate values (in m^-3 s^-1). */
  //double _collisional_rates[NUMBER_OF_IONNAMES][350];


  std::vector<std::vector<double>> _collisional_rates;
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

    // we need to find the index of the lower limit and the upper limit of
    // the temperature interval that contains the given temperature
    uint_fast32_t ilow, ihigh;

    // first handle the special cases
    if (temperature < _temperatures[0]) {
      // we will linearly extrapolate
      return 0.0;
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

        double col_rate =
        (1. - fac) * _collisional_rates[ion][ilow] + fac * _collisional_rates[ion][ihigh];


    if (temperature >=
                 _temperatures[350 - 1]) {

      col_rate = _collisional_rates[ion][350 - 1];

    }






    return std::max(col_rate, 0.);



  }

  /**
   * @brief Get the lowest tabulated temperature value.
   *
   * @return Lowest tabulated temperature value (in K).
   */
  double get_minimum_temperature() const { return _temperatures[0]; }

  double get_maximum_temperature() const { return _temperatures[350 - 1]; }
};

#endif // DERIJCKERADIATIVECOOLING_HPP
