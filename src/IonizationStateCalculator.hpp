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
 * @file IonizationStateCalculator.hpp
 *
 * @brief Class that calculates the ionization state on a grid after the photon
 * shoot loop.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef IONIZATIONSTATECALCULATOR_HPP
#define IONIZATIONSTATECALCULATOR_HPP

#include "DensityGrid.hpp"

class Abundances;
class ChargeTransferRates;
class DensitySubGrid;
class RecombinationRates;
class CollisionalRates;

/**
 * @brief Class that calculates the ionization state on a grid after the photon
 * shoot loop.
 *
 * The procedure in this class is based on Wood, K., Mathis, J. S. & Ercolano,
 * B. 2004, MNRAS, 348, 1337
 * (http://adsabs.harvard.edu/abs/2004MNRAS.348.1337W), section 4 (more
 * specifically equations (13), (14), (18), and (19)).
 */
class IonizationStateCalculator {
private:
  /*! @brief Total ionizing luminosity of all photon sources (in s^-1). */
  double _luminosity;

#ifndef HAVE_HYDROGEN_ONLY
  /*! @brief Abundances. */
  const Abundances &_abundances;
#endif

  /*! @brief Recombination rates used in ionization balance calculation. */
  const RecombinationRates &_recombination_rates;

  /*! @brief Charge transfer recombination rates used in ionization balance
   *  calculation for coolants. */
  const ChargeTransferRates &_charge_transfer_rates;

  const CollisionalRates &_collisional_rates;

public:
  IonizationStateCalculator(const double luminosity,
                            const Abundances &abundances,
                            const RecombinationRates &recombination_rates,
                            const ChargeTransferRates &charge_transfer_rates,
                            const CollisionalRates &collisional_rates);

  void
  calculate_ionization_state(const double jfac, const double hfac,
                             IonizationVariables &ionization_variables, double timestep, bool time_dependent,
                            bool do_metals) const;

  /**
   * @brief Update the total luminosity of the sources.
   *
   * @param luminosity New total luminosity for the sources (in s^-1).
   */
  inline void update_luminosity(const double luminosity) {
    _luminosity = luminosity;
  }

  static void compute_ionization_states_metals(
      const double *j_metals, const double ne, const double T, const double T4,
      const double nh0, const double nhe0, const double nhp,
      const RecombinationRates &recombination_rates,
      const ChargeTransferRates &charge_transfer_rates,
      const CollisionalRates &collisional_rates,
      IonizationVariables &ionization_variables);

  static void compute_time_dependent_metals(
      const double *j_metals, const double ne, const double T, const double T4,
      const double nh0, const double nhe0, const double nhp,
      const RecombinationRates &recombination_rates,
      const ChargeTransferRates &charge_transfer_rates,
      const CollisionalRates &collisional_rates,
      IonizationVariables &ionization_variables, double timestep);

  static void compute_ionization_states_hydrogen_helium(
      const double alphaH, const double alphaHe, const double alphaHe2, const double jH,
      const double jHe, const double nH, const double AHe, const double T,
      double &h0, double &he0, double &hep, double gammaH, double gammaHe1, double gammaHe2);

  static void compute_time_dependent_hydrogen_helium(
      const double alphaH, const double alphaHe, const double alphaHe2, const double jH,
      const double jHe, const double nH, const double AHe, const double T,
      double &h0, double &he0, double &hep, double gammaH, double gammaHe1, double gammaHe2, double ts);

  static double compute_ionization_state_hydrogen(const double alphaH,
                                                  const double jH,
                                                  const double nH,
                                                  const double gammaH,
                                                  const double old_xn,
                                                  double ts);


  static double compute_time_dependent_hydrogen(const double alphaH,
                                                  const double jH,
                                                  const double nH,
                                                  const double gammaH,
                                                  const double old_xn,
                                                  double ts);

  /**
   * @brief Functor used to calculate the ionization state of a single cell.
   */
  class IonizationStateCalculatorFunction {
  private:
    /*! @brief IonizationStateCalculator used to perform the calculation. */
    const IonizationStateCalculator &_calculator;

    /*! @brief Normalization factor used in the IonizationStateCalculator call.
     */
    const double _jfac;

    /*! @brief Normalization factor used in the IonizationStateCalculator call.
     */
    const double _hfac;


    double _timestep;

  public:
    /**
     * @brief Constructor.
     *
     * @param calculator IonizationStateCalculator used to perform the
     * calculation.
     * @param jfac Normalization factor used in the IonizationStateCalculator
     * call.
     * @param hfac Normalization factor used in the IonizationStateCalculator
     * call.
     */
    IonizationStateCalculatorFunction(
        const IonizationStateCalculator &calculator, const double jfac,
        const double hfac, double timestep)
        : _calculator(calculator), _jfac(jfac), _hfac(hfac),_timestep(timestep) {}

    /**
     * @brief Do the ionization state calculation for a single cell.
     *
     * @param cell DensityGrid::iterator pointing to a single cell in the grid.
     */
    inline void operator()(DensityGrid::iterator &cell) {
      _calculator.calculate_ionization_state(_jfac / cell.get_volume(),
                                             _hfac / cell.get_volume(),
                                             cell.get_ionization_variables(),_timestep, 
                                             false, true);
    }
  };

  void
  calculate_ionization_state(const double totweight, DensityGrid &grid,
                             std::pair< cellsize_t, cellsize_t > &block, double timestep) const;

  void calculate_ionization_state(const double totweight,
                                  DensitySubGrid &subgrid, double timestep, bool time_dependent, bool do_metals) const;
};

#endif // IONIZATIONSTATECALCULATOR_HPP
