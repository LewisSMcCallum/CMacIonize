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
 * @file ChiantiRecombinationRates.hpp
 *
 * @brief RecombinationRates implementation with Chianti's recombination rates:
 * header.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CHIANTIRECOMBINATIONRATES_HPP
#define CHIANTIRECOMBINATIONRATES_HPP

#include "RecombinationRates.hpp"
#include <vector>

/**
 * @brief RecombinationRates implementation with Chianti's recombination rates.
 *

 */
class ChiantiRecombinationRates : public RecombinationRates {
private:
  /*! @brief Temperature values (in K). */
  std::vector<double> _temperatures;

  /*! @brief Collisional rate values (in m^-3 s^-1). */
  //double _collisional_rates[NUMBER_OF_IONNAMES][250];


  std::vector<std::vector<double>> _recomb_rates;
  /*! @brief Logarithm of the minimum tabulated temperature
   *  (in log(T / K)). */
  double _min_logT;

  /*! @brief Inverse of the average logarithmic distance between tabulated
   *  temperature values (in 1 / log(T / K)). */
  double _inverse_avg_dlogT;


public:
  ChiantiRecombinationRates();

  double get_recombination_rate_chianti(const int ion, const double temperature) const;

  virtual double get_recombination_rate(const int_fast32_t ion,
                                        const double temperature) const;
};

#endif // CHIANTIRECOMBINATIONRATES_HPP
