/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016, 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file IonizationStateCalculator.cpp
 *
 * @brief IonizationStateCalculator implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "IonizationStateCalculator.hpp"
#include "Abundances.hpp"
#include "ChargeTransferRates.hpp"
#include "CollisionalRates.hpp"
#include "DensityGrid.hpp"
#include "DensityGridTraversalJobMarket.hpp"
#include "DensitySubGrid.hpp"
#include "DensityValues.hpp"
#include "Error.hpp"
#include "RecombinationRates.hpp"
#include "WorkDistributor.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

/**
 * @brief Constructor.
 *
 * @param luminosity Total ionizing luminosity of all photon sources (in s^-1).
 * @param abundances Abundances.
 * @param recombination_rates RecombinationRates used in ionization balance
 * calculation.
 * @param charge_transfer_rates ChargeTransferRate used in ionization balance
 * calculation for coolants.
 */
IonizationStateCalculator::IonizationStateCalculator(
    double luminosity, const Abundances &abundances,
    const RecombinationRates &recombination_rates,
    const ChargeTransferRates &charge_transfer_rates,
    const CollisionalRates &collisional_rates)
    : _luminosity(luminosity),
#ifndef HAVE_HYDROGEN_ONLY
      _abundances(abundances),
#endif
      _recombination_rates(recombination_rates),
      _charge_transfer_rates(charge_transfer_rates),
      _collisional_rates(collisional_rates) {
}

/**
 * @brief Does the ionization state calculation for a single cell.
 *
 * @param jfac Normalization factor for the mean intensity integrals in this
 * cell.
 * @param hfac Normalization factor for the heating integrals in this cell.
 * @param ionization_variables Ionization variables for the cell we operate on.
 */
void IonizationStateCalculator::calculate_ionization_state(
    const double jfac, const double hfac,
    IonizationVariables &ionization_variables, double timestep) const {

  // normalize the mean intensity integrals
  const double jH = jfac * ionization_variables.get_mean_intensity(ION_H_n);
  //cmac_assert_message(jH >= 0., "jH: %g, jfac: %g, mean_intensity: %g", jH,
    //                  jfac, ionization_variables.get_mean_intensity(ION_H_n));

#ifdef HAS_HELIUM
  const double jHe = jfac * ionization_variables.get_mean_intensity(ION_He_n);
  //cmac_assert(jHe >= 0.);
#endif

  // normalize the heating integrals (for explicit heating in RHD)
  const double hH = hfac * ionization_variables.get_heating(HEATINGTERM_H);
  ionization_variables.set_heating(HEATINGTERM_H, hH);
#ifdef HAS_HELIUM
  const double hHe = hfac * ionization_variables.get_heating(HEATINGTERM_He);
  ionization_variables.set_heating(HEATINGTERM_He, hHe);
#endif

  // get the number density
  const double ntot = ionization_variables.get_number_density();
  cmac_assert(ntot >= 0.);

  // find the ionization equilibrium for hydrogen and helium
  if (ntot > 0.) {
    const double T = ionization_variables.get_temperature();
    const double alphaH =
        _recombination_rates.get_recombination_rate(ION_H_n, T);
    const double gammaH =
       _collisional_rates.get_collisional_rate(ION_H_n, T);

    cmac_assert(alphaH >= 0.);

#ifdef HAS_HELIUM
#ifdef VARIABLE_ABUNDANCES
    const double AHe =
        ionization_variables.get_abundances().get_abundance(ELEMENT_He);
#else
    const double AHe = _abundances.get_abundance(ELEMENT_He);
#endif
    // h0find
    double h0 = 0.;
    double he0 = 0.;
    double hep = 0.;
    if (AHe != 0.) {
      const double alphaHe =
          _recombination_rates.get_recombination_rate(ION_He_n, T);
      const double gammaHe1 = _collisional_rates.get_collisional_rate(ION_He_n, T);

      const double gammaHe2 = _collisional_rates.get_collisional_rate(ION_He_p1, T);
      const double alphaHe2 = _recombination_rates.get_recombination_rate(ION_He_p1, T);
      compute_ionization_states_hydrogen_helium(alphaH, alphaHe,alphaHe2, jH, jHe, ntot,
                                                AHe, T, h0, he0, hep, gammaH, gammaHe1,
                                                 gammaHe2);


    } else {
      h0 = compute_ionization_state_hydrogen(alphaH, jH, ntot, gammaH, ionization_variables.get_prev_ionic_fraction(ION_H_n), timestep);
    }
#else
    const double h0 = compute_ionization_state_hydrogen(alphaH, jH, ntot, gammaH, ionization_variables.get_prev_ionic_fraction(ION_H_n), timestep);
#endif

    ionization_variables.set_ionic_fraction(ION_H_n, h0);

#ifdef HAS_HELIUM
    ionization_variables.set_ionic_fraction(ION_He_n, he0);
    ionization_variables.set_ionic_fraction(ION_He_p1,hep);
#endif

    // do the coolants
    const double nhp = ntot * (1. - h0);
#ifdef HAS_HELIUM
    const double ne = ntot*(1-h0) + 2.0*AHe*ntot*(1-he0-hep) + ntot*hep*AHe;
#else
    const double ne = nhp;
#endif
    const double T4 = T * 1.e-4;

    const double j_metals[12] = {
#ifdef HAS_CARBON
        jfac*ionization_variables.get_mean_intensity(ION_C_p1),
        jfac*ionization_variables.get_mean_intensity(ION_C_p2),
#else
        0., 0.,
#endif
#ifdef HAS_NITROGEN
        jfac*ionization_variables.get_mean_intensity(ION_N_n),
        jfac*ionization_variables.get_mean_intensity(ION_N_p1),
        jfac*ionization_variables.get_mean_intensity(ION_N_p2),
#else
        0., 0.,
        0.,
#endif
#ifdef HAS_OXYGEN
        jfac*ionization_variables.get_mean_intensity(ION_O_n),
        jfac*ionization_variables.get_mean_intensity(ION_O_p1),
#else
        0., 0.,
#endif
#ifdef HAS_NEON
        jfac*ionization_variables.get_mean_intensity(ION_Ne_n),
        jfac*ionization_variables.get_mean_intensity(ION_Ne_p1),
#else
        0., 0.,
#endif
#ifdef HAS_SULPHUR
        jfac*ionization_variables.get_mean_intensity(ION_S_p1),
        jfac*ionization_variables.get_mean_intensity(ION_S_p2),
        jfac*ionization_variables.get_mean_intensity(ION_S_p3)
#else
        0., 0.,
        0.
#endif
    };

    const double nh0 = ntot * h0;
#ifdef HAS_HELIUM
    const double nhe0 = ntot * he0 * AHe;
#else
    const double nhe0 = 0.;
#endif

    compute_ionization_states_metals(
        j_metals, ne, T, T4, nh0, nhe0, nhp, _recombination_rates,
        _charge_transfer_rates, _collisional_rates, ionization_variables);



  } else {

    // vacuum cell: set all values to 0
      ionization_variables.set_ionic_fraction(ION_H_n, 0.);

#ifdef HAS_HELIUM
      ionization_variables.set_ionic_fraction(ION_He_n, 0.);
#endif

#ifdef HAS_CARBON
      ionization_variables.set_ionic_fraction(ION_C_p1, 0.);
      ionization_variables.set_ionic_fraction(ION_C_p2, 0.);
#endif

#ifdef HAS_NITROGEN
      ionization_variables.set_ionic_fraction(ION_N_n, 0.);
      ionization_variables.set_ionic_fraction(ION_N_p1, 0.);
      ionization_variables.set_ionic_fraction(ION_N_p2, 0.);
#endif

#ifdef HAS_OXYGEN
      ionization_variables.set_ionic_fraction(ION_O_n, 0.);
      ionization_variables.set_ionic_fraction(ION_O_p1, 0.);
#endif

#ifdef HAS_NEON
      ionization_variables.set_ionic_fraction(ION_Ne_n, 0.);
      ionization_variables.set_ionic_fraction(ION_Ne_p1, 0.);
#endif

#ifdef HAS_SULPHUR
      ionization_variables.set_ionic_fraction(ION_S_p1, 0.);
      ionization_variables.set_ionic_fraction(ION_S_p2, 0.);
      ionization_variables.set_ionic_fraction(ION_S_p3, 0.);
#endif

  }

  cmac_assert(ionization_variables.get_ionic_fraction(ION_H_n) >= 0.);

#ifdef DO_OUTPUT_PHOTOIONIZATION_RATES
  // set the mean intensity values to the values in correct physical units
  ionization_variables.set_mean_intensity(ION_H_n, jH);
#ifdef HAS_HELIUM
  ionization_variables.set_mean_intensity(ION_He_n, jHe);
#endif
#endif
}

/**
 * @brief Compute the ionization balance for the metals at the given temperature
 * (and using the given ionizing luminosity integrals).
 *
 *
 * the procedure is always the same: the total density for an element X with
 * ionization states \f$X^0, X^+, X^{2+},...\f$ is
 * \f[
 *   n(X) = n(X^0) + n(X^+) + n(X^{2+}) + ...,
 * \f]
 * while the ionization balance for each ion is given by
 * \f[
 *   n(X^+)R(X^+) = n(X^0)I(X^+),
 * \f]
 * where \f$R(X^+)\f$ is the recombination rate from level \f$X^+\f$ to level
 * \f$X^0\f$, and \f$I(X^+)\f$ is the ionization rate from level \f$X^0\f$ to
 * level \f$X^+\f$.
 *
 * This can be rewritten as
 * \f[
 *   n(X^+) = \frac{n(X^0)I(X^+)}{R(X^+)} = n(X^0)C(X^+).
 * \f]
 * Recombination from \f$X^{2+}\f$ to \f$X^0\f$ happens in two stages, so the
 * recombination rate from \f$X^{2+}\f$ to \f$X^0\f$ is the product of the
 * recombination rates from \f$X^{2+}\f$ to \f$X^+\f$ and from \f$X^+\f$ to
 * \f$X^0\f$.
 *
 * We want the ionic fractions \f$\frac{n(X^+)}{n(X)}\f$, so
 * \f{eqnarray*}{
 *  \frac{n(X^+)}{n(X)} &=& \frac{n(X^0)C(X^+)}{n(X^0) + n(X^+) + n(X^{2+}) +
 *                          ...} \\
 *                      &=& \frac{n(X^0)C(X^+)}{n(X^0) + n(X^0)C(X^+) +
 *                          n(X^+)C(X^{2+}) + ...} \\
 *                      &=& \frac{n(X^0)C(X^+)}{n(X^0) + n(X^0)C(X^+) +
 *                          n(X^0)C(X^+)C(X^{2+}) + ...} \\
 *                      &=& \frac{C(X^+)}{1 + C(X^+) + C(X^+)C(X^{2+}) + ...}.
 * \f}
 *
 * @param j_metals Ionizing luminosity integrals for the metal ions (in s^-1).
 * @param ne Number density of electrons (in m^-3).
 * @param T Temperature (in K).
 * @param T4 Temperature (in 10^4 K).
 * @param nh0 Number density of neutral hydrogen (in m^-3).
 * @param nhe0 Number density of neutral helium (in m^-3).
 * @param nhp Number density of ionized hydrogen (in m^-3).
 * @param recombination_rates RecombinationRates.
 * @param charge_transfer_rates ChargeTransferRates.
 * @param ionization_variables IonizationStateVariables to operate on.
 */
void IonizationStateCalculator::compute_ionization_states_metals(
    const double *j_metals, const double ne, const double T, const double T4,
    const double nh0, const double nhe0, const double nhp,
    const RecombinationRates &recombination_rates,
    const ChargeTransferRates &charge_transfer_rates,
    const CollisionalRates &collisional_rates,
    IonizationVariables &ionization_variables) {

#ifdef HAS_CARBON
  const double jCp1 = j_metals[0];
  const double jCp2 = j_metals[1];
#endif

#ifdef HAS_NITROGEN
  const double jNn = j_metals[2];
  const double jNp1 = j_metals[3];
  const double jNp2 = j_metals[4];
#endif

#ifdef HAS_OXYGEN
  const double jOn = j_metals[5];
  const double jOp1 = j_metals[6];
#endif

#ifdef HAS_NEON
  const double jNen = j_metals[7];
  const double jNep1 = j_metals[8];
#endif

#ifdef HAS_SULPHUR
  const double jSp1 = j_metals[9];
  const double jSp2 = j_metals[10];
  const double jSp3 = j_metals[11];
#endif

#ifdef HAS_CARBON
  const double alphaC[2] = {
      recombination_rates.get_recombination_rate(ION_C_p1, T),
      recombination_rates.get_recombination_rate(ION_C_p2, T)};
#endif

#ifdef HAS_NITROGEN
  const double alphaN[3] = {
      recombination_rates.get_recombination_rate(ION_N_n, T),
      recombination_rates.get_recombination_rate(ION_N_p1, T),
      recombination_rates.get_recombination_rate(ION_N_p2, T)};
#endif

#ifdef HAS_OXYGEN
  const double alphaO[4] = {
      recombination_rates.get_recombination_rate(ION_O_n, T),
      recombination_rates.get_recombination_rate(ION_O_p1, T),
      recombination_rates.get_recombination_rate(ION_O_p2, T),
      recombination_rates.get_recombination_rate(ION_O_p3, T)};
#endif

#ifdef HAS_NEON
  const double alphaNe[4] = {
      recombination_rates.get_recombination_rate(ION_Ne_n, T),
      recombination_rates.get_recombination_rate(ION_Ne_p1, T),
      recombination_rates.get_recombination_rate(ION_Ne_p2, T),
      recombination_rates.get_recombination_rate(ION_Ne_p3, T)};
#endif

#ifdef HAS_SULPHUR
  const double alphaS[3] = {
      recombination_rates.get_recombination_rate(ION_S_p1, T),
      recombination_rates.get_recombination_rate(ION_S_p2, T),
      recombination_rates.get_recombination_rate(ION_S_p3, T)};
#endif

#ifdef HAS_CARBON
  // carbon
  // the charge transfer recombination rates for C+ are negligble

  const double C21 = (jCp1 + ne*collisional_rates.get_collisional_rate(ION_C_p1, T)) / (ne * alphaC[0]);
  const double C32 =
      (jCp2 + ne*collisional_rates.get_collisional_rate(ION_C_p2, T)) /
      (ne * alphaC[1] +
       nh0 * charge_transfer_rates.get_charge_transfer_recombination_rate_H(
                 ION_C_p2, T4) +
       nhe0 * charge_transfer_rates.get_charge_transfer_recombination_rate_He(
                  ION_C_p2, T4));
  const double C31 = C32 * C21;
  const double sumC_inv = 1. / (1. + C21 + C31);
  ionization_variables.set_ionic_fraction(ION_C_p1, C21 * sumC_inv);
  ionization_variables.set_ionic_fraction(ION_C_p2, C31 * sumC_inv);
#endif

#ifdef HAS_NITROGEN
  // nitrogen
  const double N21 =
      (jNn + nhp * charge_transfer_rates.get_charge_transfer_ionization_rate_H(
                       ION_N_n, T4) + ne*collisional_rates.get_collisional_rate(ION_N_n, T)) /
      (ne * alphaN[0] +
       nh0 * charge_transfer_rates.get_charge_transfer_recombination_rate_H(
                 ION_N_n, T4));
  const double N32 =
      (jNp1 + ne*collisional_rates.get_collisional_rate(ION_N_p1, T)) /
      (ne * alphaN[1] +
       nh0 * charge_transfer_rates.get_charge_transfer_recombination_rate_H(
                 ION_N_p1, T4) +
       nhe0 * charge_transfer_rates.get_charge_transfer_recombination_rate_He(
                  ION_N_p1, T4));
  const double N43 =
      (jNp2 + ne*collisional_rates.get_collisional_rate(ION_N_p2, T)) /
      (ne * alphaN[2] +
       nh0 * charge_transfer_rates.get_charge_transfer_recombination_rate_H(
                 ION_N_p2, T4) +
       nhe0 * charge_transfer_rates.get_charge_transfer_recombination_rate_He(
                  ION_N_p2, T4));
  const double N31 = N32 * N21;
  const double N41 = N43 * N31;
  const double sumN_inv = 1. / (1. + N21 + N31 + N41);
  ionization_variables.set_ionic_fraction(ION_N_n, N21 * sumN_inv);
  ionization_variables.set_ionic_fraction(ION_N_p1, N31 * sumN_inv);
  ionization_variables.set_ionic_fraction(ION_N_p2, N41 * sumN_inv);
#endif

#ifdef HAS_OXYGEN
  // Oxygen

  const double O21 =
      (jOn + nhp * charge_transfer_rates.get_charge_transfer_ionization_rate_H(
                       ION_O_n, T4) +
                     ne*collisional_rates.get_collisional_rate(ION_O_n, T)) /
      (ne * alphaO[0] +
       nh0 * charge_transfer_rates.get_charge_transfer_recombination_rate_H(
                 ION_O_n, T4));
  const double O32 =
      (jOp1 + ne*collisional_rates.get_collisional_rate(ION_O_p1, T)) /
      (ne * alphaO[1] +
       nh0 * charge_transfer_rates.get_charge_transfer_recombination_rate_H(
                 ION_O_p1, T4) +
       nhe0 * charge_transfer_rates.get_charge_transfer_recombination_rate_He(
                  ION_O_p1, T4));


  const double O43 =
     (ne*collisional_rates.get_collisional_rate(ION_O_p2, T))/
     (ne*alphaO[2] +
       nh0 * charge_transfer_rates.get_charge_transfer_recombination_rate_H(
                 ION_O_p2, T4) +
       nhe0 * charge_transfer_rates.get_charge_transfer_recombination_rate_He(
                  ION_O_p2, T4));


  const double O54 =
      (ne*collisional_rates.get_collisional_rate(ION_O_p3, T))/
            (ne*alphaO[3] +
          nh0 * charge_transfer_rates.get_charge_transfer_recombination_rate_H(
                                 ION_O_p3, T4) +
            nhe0 * charge_transfer_rates.get_charge_transfer_recombination_rate_He(
                        ION_O_p3, T4));


  const double O31 = O32 * O21;
  const double O41 = O43 * O31;
  const double O51 = O54 * O41;
  const double sumO_inv = 1. / (1. + O21 + O31 + O41 + O51);



  ionization_variables.set_ionic_fraction(ION_O_n, O21 * sumO_inv);
  ionization_variables.set_ionic_fraction(ION_O_p1, O31 * sumO_inv);
  ionization_variables.set_ionic_fraction(ION_O_p2, O41 * sumO_inv);
  ionization_variables.set_ionic_fraction(ION_O_p3, O51 * sumO_inv);
#endif

#ifdef HAS_NEON
  // Neon
  const double Ne21 = (jNen + ne*collisional_rates.get_collisional_rate(ION_Ne_n, T)) / (ne * alphaNe[0]);
  const double Ne32 =
      (jNep1 + ne*collisional_rates.get_collisional_rate(ION_Ne_p1, T)) /
      (ne * alphaNe[1] +
       nh0 * charge_transfer_rates.get_charge_transfer_recombination_rate_H(
                 ION_Ne_p1, T4) +
       nhe0 * charge_transfer_rates.get_charge_transfer_recombination_rate_He(
                  ION_Ne_p1, T4));
  const double Ne43 = (ne*collisional_rates.get_collisional_rate(ION_Ne_p2, T)) /
       (ne*alphaNe[2] +
         nh0 * charge_transfer_rates.get_charge_transfer_recombination_rate_H(
                   ION_Ne_p2, T4) +
         nhe0 * charge_transfer_rates.get_charge_transfer_recombination_rate_He(
                    ION_Ne_p2, T4));

  const double Ne54 = (ne*collisional_rates.get_collisional_rate(ION_Ne_p3, T)) /
        (ne*alphaNe[3]);
  const double Ne31 = Ne32 * Ne21;
  const double Ne41 = Ne43 * Ne31;
  const double Ne51 = Ne54 * Ne41;
  const double sumNe_inv = 1. / (1. + Ne21 + Ne31 + Ne41 + Ne51);

  //if (std::isnan(sumNe_inv) || std::isinf(sumNe_inv)) {
  //  std::cout << "FOR T OF" << T << " " << Ne21 << " " << jNen << " " << alphaNe[0] << " " << ne << std::endl;
  //  std::cout << "COL RATE " << collisional_rates.get_collisional_rate(ION_Ne_n, T) << std::endl;
  //}

  ionization_variables.set_ionic_fraction(ION_Ne_n, Ne21 * sumNe_inv);
  ionization_variables.set_ionic_fraction(ION_Ne_p1, Ne31 * sumNe_inv);
  ionization_variables.set_ionic_fraction(ION_Ne_p2, Ne41 * sumNe_inv);
  ionization_variables.set_ionic_fraction(ION_Ne_p3, Ne51 *sumNe_inv);
#endif

#ifdef HAS_SULPHUR
  // Sulphur
  const double S21 =
      (jSp1 + ne*collisional_rates.get_collisional_rate(ION_S_p1, T)) /
      (ne * alphaS[0] +
       nh0 * charge_transfer_rates.get_charge_transfer_recombination_rate_H(
                 ION_S_p1, T4));
  const double S32 =
      (jSp2 + ne*collisional_rates.get_collisional_rate(ION_S_p2, T)) /
      (ne * alphaS[1] +
       nh0 * charge_transfer_rates.get_charge_transfer_recombination_rate_H(
                 ION_S_p2, T4) +
       nhe0 * charge_transfer_rates.get_charge_transfer_recombination_rate_He(
                  ION_S_p2, T4));
  const double S43 =
      (jSp3 + ne*collisional_rates.get_collisional_rate(ION_S_p3, T)) /
      (ne * alphaS[2] +
       nh0 * charge_transfer_rates.get_charge_transfer_recombination_rate_H(
                 ION_S_p3, T4) +
       nhe0 * charge_transfer_rates.get_charge_transfer_recombination_rate_He(
                  ION_S_p3, T4));
  const double S31 = S32 * S21;
  const double S41 = S43 * S31;
  const double sumS_inv = 1. / (1. + S21 + S31 + S41);
  ionization_variables.set_ionic_fraction(ION_S_p1, S21 * sumS_inv);
  ionization_variables.set_ionic_fraction(ION_S_p2, S31 * sumS_inv);
  ionization_variables.set_ionic_fraction(ION_S_p3, S41 * sumS_inv);
#endif
}

/**
 * @brief Solves the ionization and temperature equations based on the values of
 * the mean intensity integrals in each cell.
 *
 * @param totweight Total weight off all photons used.
 * @param grid DensityGrid for which the calculation is done.
 * @param block Block that should be traversed by the local MPI process.
 */
void IonizationStateCalculator::calculate_ionization_state(
    const double totweight, DensityGrid &grid,
    std::pair< cellsize_t, cellsize_t > &block, double timestep) const {

  // compute the normalization factor for the mean intensity integrals, which
  // depends on the total weight of all photons, and on the volume of each cell
  // the volume of the cell is taken into account on a cell level, since cells
  // don't necessarily have the same volume
  const double jfac = _luminosity / totweight;
  const double hfac =
      jfac * PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK);
  WorkDistributor<
      DensityGridTraversalJobMarket< IonizationStateCalculatorFunction >,
      DensityGridTraversalJob< IonizationStateCalculatorFunction > >
      workers;
  IonizationStateCalculatorFunction do_calculation(*this, jfac, hfac, timestep);
  DensityGridTraversalJobMarket< IonizationStateCalculatorFunction > jobs(
      grid, do_calculation, block);
  workers.do_in_parallel(jobs);
}

/**
 * @brief Calculate the ionization state for all cells in the given subgrid.
 *
 * @param totweight Total weight of all photon packets.
 * @param subgrid DensitySubGrid to work on.
 */
void IonizationStateCalculator::calculate_ionization_state(
    const double totweight, DensitySubGrid &subgrid, double timestep) const {


  double jfac;
  double hfac;

  if (totweight == 0 || _luminosity == 0){
    jfac = 0.0;
    hfac = 0.0;

  } else {

    jfac = _luminosity / totweight;
    hfac =
        jfac * PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK);

  }


  for (auto cellit = subgrid.begin(); cellit != subgrid.end(); ++cellit) {
    calculate_ionization_state(jfac / cellit.get_volume(),
                               hfac / cellit.get_volume(),
                               cellit.get_ionization_variables(),timestep);
  }
}

/**
 * @brief Iteratively find the neutral fractions of hydrogen and helium based on
 * the given value current values, the values of the intensity integrals and the
 * recombination rate.
 *
 * The equation for the ionization balance of hydrogen is
 * @f[
 *   n({\rm{}H}^0)\int_{\nu{}_i}^\infty{} \frac{4\pi{}J_\nu{}}{h\nu{}}
 *   a_\nu{}({\rm{}H}^0) {\rm{}d}\nu{} = n_e n({\rm{}H}^+)
 *   \alpha{}({\rm{}H}^0, T_e) - n_e P({\rm{}H}_{\rm{}OTS})n({\rm{}He}^+)
 *   \alpha{}_{2^1{\rm{}P}}^{\rm{}eff},
 * @f]
 * and that of helium
 * @f[
 *   n({\rm{}He}^0) \int_{\nu{}_i}^\infty{} \frac{4\pi{}J_\nu{}}{h\nu{}}
 *   a_\nu{}({\rm{}He}^0) {\rm{}d}\nu{} = n({\rm{}He}^0) n_e
 *   \alpha{}({\rm{}He}^0, T_e),
 * @f]
 * where the value of the integral on the left hand side of the equations, the
 * temperature @f$T_e@f$, the recombination rates @f$\alpha{}@f$, and the
 * current values of @f$n({\rm{}H}^0)@f$ and @f$n({\rm{}He}^0)@f$ (and hence all
 * other densities) are given. We want to determine the new values for the
 * densities.
 *
 * We start with helium. First, we derive the following expression for the
 * electron density @f$n_e@f$:
 * @f[
 *   n_e = (1-n({\rm{}H}^0)) n_{\rm{}H} + (1-n({\rm{}He}^0)) n_{\rm{}He},
 * @f]
 * which follows from charge conservation (assuming other elements do not
 * contribute free electrons due to their low abundances). Since the total
 * density we store in every cell is the hydrogen density, and abundances are
 * expressed relative to the hydrogen density, we can rewrite this as
 * @f[
 *   n_e = [(1-n({\rm{}H}^0)) + (1-n({\rm{}He}^0)) A_{\rm{}He}] n_{\rm{}tot}.
 * @f]
 * Using this, we can rewrite the helium ionization balance as
 * @f[
 *   n({\rm{}He}^0) = C_{\rm{}He} (1-n({\rm{}He}^0)) \frac{n_e}{n_{\rm{}tot}},
 * @f]
 * with @f$C_{\rm{}He} = \frac{\alpha{}({\rm{}He}^0, T_e) n_{\rm{}tot} }
 * {J_{\rm{}He}}@f$, a constant for a given temperature. Due to the @f$n_e@f$
 * factor, this equation is coupled to the ionization balance of hydrogen.
 * However, if we assume the hydrogen neutral fraction to be known, we can find
 * a closed expression for the helium neutral fraction:
 * @f[
 *   n({\rm{}He}^0) = \frac{-D - \sqrt{D^2 - 4A_{\rm{}He}C_{\rm{}He}B}}
 *   {2A_{\rm{}He}C_{\rm{}He}},
 * @f]
 * where
 * @f[
 *   D = -1 - 2A_{\rm{}He}C_{\rm{}He} - C_{\rm{}He} + C_{\rm{}He}n({\rm{}H}^0),
 * @f]
 * and
 * @f[
 *   B = C_{\rm{}He} + C_{\rm{}He}n({\rm{}H}^0) - A_{\rm{}He} C_{\rm{}He}.
 * @f]
 * This expression can be found by solving the quadratic equation in
 * @f$n({\rm{}He}^0)@f$. We choose the minus sign based on the fact that
 * @f$D < 0@f$ and the requirement that the helium neutral fraction be positive.
 *
 * To find the hydrogen neutral fraction, we have to address the fact that the
 * probability of a  Ly@f$\alpha{}@f$
 * photon being absorbed on the spot (@f$P({\rm{}H}_{\rm{}OTS})@f$) also depends
 * on the neutral fractions. We therefore use an iterative scheme to find the
 * hydrogen neutral fraction. The equation we want to solve is
 * @f[
 *   n({\rm{}H}^0) = C_{\rm{}H} (1-n({\rm{}H}^0)) \frac{n_e}{n_{\rm{}tot}},
 * @f]
 * which (conveniently) has the same form as the equation for helium. But know
 * we have
 * @f[
 *   C_{\rm{}H} = C_{{\rm{}H},1} - C_{{\rm{}H},2} P({\rm{}H}_{\rm{}OTS})
 *   \frac{1-n({\rm{}He}^0)}{1-n({\rm{}H}^0)},
 * @f]
 * with @f$C_{{\rm{}H},1} = \frac{\alpha{}({\rm{}H}^0, T_e) n_{\rm{}tot} }
 * {J_{\rm{}H}}@f$ and @f$C_{{\rm{}H},2} = \frac{A_{\rm{}He}
 * \alpha{}_{2^1{\rm{}P}}^{\rm{}eff} n_{\rm{}tot} }
 * {J_{\rm{}H}}@f$ two constants.
 *
 * For every iteration of the scheme, we calculate an approximate value for
 * @f$C_{\rm{}H}@f$, based on the neutral fractions obtained during the previous
 * iteration. We also calculate the helium neutral fraction, based on the
 * hydrogen neutral fraction from the previous iteration. Using these values,
 * we can then solve for the hydrogen neutral fraction. We repeat the process
 * until the relative difference between the obtained neutral fractions is
 * below some tolerance value.
 *
 * @param alphaH Hydrogen recombination rate (in m^3s^-1).
 * @param alphaHe Helium recombination rate (in m^3s^-1).
 * @param jH Hydrogen intensity integral (in s^-1).
 * @param jHe Helium intensity integral (in s^-1).
 * @param nH Hydrogen number density (in m^-3).
 * @param AHe Helium abundance @f$A_{\rm{}He}@f$ (relative w.r.t. hydrogen).
 * @param T Temperature (in K).
 * @param h0 Variable to store resulting hydrogen neutral fraction in.
 * @param he0 Variable to store resulting helium neutral fraction in.
 */
void IonizationStateCalculator::compute_ionization_states_hydrogen_helium(
    const double alphaH, const double alphaHe, const double alphaHe2, const double jH,
    const double jHe, const double nH, const double AHe, const double T,
    double &h0, double &he0, double &hep, const double gammaH, const double gammaHe1,
    const double gammaHe2) {


  if ((jH == 0) && (jHe == 0) && (gammaH < 1e-25) && (gammaHe1 < 1e-25) && (gammaHe2 < 1e-25)) {
    h0 = 0.999999;
    he0 = 0.999999;
    hep = 0.0;
    return;
  }



  // make sure the input to this function is physical
  cmac_assert(alphaH >= 0.);
  cmac_assert(alphaHe >= 0.);
  //cmac_assert(jH >= 0.);
  //cmac_assert(jHe >= 0.);
  cmac_assert(nH >= 0.);
  cmac_assert(AHe >= 0.);
  cmac_assert(T >= 0.);

  // shortcut: if jH is very small, then the gas is neutral

  //got rid of this because not the case any more
  //if (jH < 1.e-20) {
  //  h0 = 1.;
  //  he0 = 1.;
  //  return;
  //}

  // we multiplied Kenny's value with 1.e-6 to convert from cm^3s^-1 to m^3s^-1
  // NOTE that this is a different expression from the one in Kenny's code!
  const double alpha_e_2sP = 4.17e-20 * std::pow(T * 1.e-4, -0.861);




  // initial guesses for the neutral fractions
  double h0old = 0.99*(1. - std::exp(-1.*(jH+nH*gammaH)/(2.*alphaH*nH)));

  if (h0old != h0old) {
    h0old = 0.9;
  }

  cmac_assert(h0old >= 0. && h0old <= 1.);

  // by enforcing a relative difference of 10%, we make sure we have at least
  // one iteration
  h0 = 0.9 * h0old;

  double he0old = 0.5/(alphaHe*nH/(gammaHe1*nH + jHe));
  he0old = std::min(he0old, 1.);

  if (he0old != he0old) {
    he0old = 0.9;
  }

  // again, by using this value we make sure we have at least one iteration
  he0 = 0.9*he0old;
  uint_fast8_t niter = 0;

  while (std::abs(h0 - h0old) > 1.e-4 * h0old &&
         std::abs(he0 - he0old) > 1.e-4 * he0old) {
    ++niter;
    h0old = h0;
    he0old = he0;




//calculate Helium neutral fraction and helium+ fraction

    const double Bhe = nH*(1-h0 + AHe*alphaHe2/(alphaHe2+gammaHe2) + 2*gammaHe2*AHe/(alphaHe2+gammaHe2));
    const double Che = alphaHe2*AHe/(gammaHe2 + alphaHe2) + 2*gammaHe2*AHe/(alphaHe2+gammaHe2);
    const double Dhe = Che*nH*alphaHe*alphaHe2/(alphaHe2+gammaHe2) + nH*Che*gammaHe1;
    const double Ehe = -Che*nH*alphaHe*alphaHe2/(alphaHe2+gammaHe2) - Bhe*alphaHe*alphaHe2/(alphaHe2+gammaHe2) - Bhe*gammaHe1 - jHe;
    const double Fhe = Bhe*alphaHe*alphaHe2/(alphaHe2+gammaHe2);

    he0 = (-1.0*Ehe - std::sqrt(Ehe * Ehe - 4. * Dhe*Fhe)) /
          (2. * Dhe);

    const double hepp = (1.0-he0)*gammaHe2/(alphaHe2+gammaHe2);

    hep = (1.0 - he0 - hepp);

    double ne = nH*(1 - h0) + 2*hepp*AHe*nH + hep*AHe*nH;



    // calculate a new guess for C_H
    const double pHots = 1. / (1. + 77. * he0 / std::sqrt(T) / h0old);
    // make sure pHots is not NaN
    cmac_assert(pHots == pHots);



    // find the hydrogen neutral fraction -

    double ots = AHe*ne*hep*pHots*alpha_e_2sP;

    h0 = (ne*alphaH - ots)/(jH + ne*alphaH + ne*gammaH);

    if (niter > 20) {
      // if we have a lot of iterations: use the mean value to speed up
      // convergence
      h0 = 0.5 * (h0 + h0old);
      he0 = 0.5 * (he0 + he0old);
    }
    if (niter > 100) {
      if (h0 < 1e-6) {
        break;
      } else {
        std::cout << "UH OH " << std::endl;
        std::cout << alphaH << " " << alphaHe << " " << alphaHe2 << " " << gammaH << " " << gammaHe1 << " " << alphaHe2 << std::endl;
        std::cout << jH << " " << jHe << std::endl;
        std::cout << T << " " << nH << std::endl;
        std::cout << h0 << " " << he0 << std::endl;
        cmac_error("Too many iterations in ionization loop!");

      }

    }
  }
  if (h0 != h0){
    std::cout << alphaH << " " << alphaHe << " " << alphaHe2 << " " << gammaH << " " << gammaHe1 << " " << gammaHe2 << std::endl;
    std::cout << jH << " " << jHe << std::endl;
    std::cout << T << " " << nH << std::endl;
    std::cout << h0 << " " << he0 << std::endl;
    cmac_error("SHTAP");
  }

}

/**
 * @brief find_H0() for a system without helium.
 *
 * We do not need to iterate in this case: the solution is simply given by the
 * solution of a quadratic equation. This can be derived as follows.
 *
 * The ionization balance equation for hydrogen is given by
 * @f[
 *   n_{\rm{}H}^2 \left(1 - x_{\rm{}H}\right)^2 \alpha{}_{\rm{}H} =
 *     n_{\rm{}H} x_{\rm{}H} J_{\rm{}H}.
 * @f]
 *
 * This can be rewritten as
 * @f[
 *   x_{\rm{}H}^2 - \left(2 + C_{\rm{}H}\right) x_{\rm{}H} + 1 = 0,
 * @f]
 * with @f$C_{\rm{}H} = \frac{J_{\rm{}H}}{n_{\rm{}H}\alpha{}_{\rm{}H}}@f$.
 *
 * The solutions of this equation are
 * @f[
 *   x_{\rm{}H} = 1 + \frac{1}{2}C_{\rm{}H} \pm{} \sqrt{\left(1 + \frac{1}{2}
 *     C_{\rm{}H}\right)^2 - 1}.
 * @f]
 *
 * Since @f$C_{\rm{}H} > 0@f$ and @f$0 \leq{} x_{\rm{}H} \leq{} 1@f$, we choose
 * the solution with the negative square root.
 *
 * For normal values of @f$C_{\rm{}H}@f$ (meaning: not too large), the
 * expression reduces to
 * @f[
 *   x_{\rm{}H} = 1 + \frac{1}{2}C_{\rm{}H} \left(1 - \sqrt{\frac{4}{C_{\rm{}H}}
 *     + 1}\right).
 * @f]
 *
 * For very large values of @f$C_{\rm{}H}@f$, the term in between the square
 * root becomes very small and round off error can lead the expression to
 * wrongly evaluate to @f$x_{\rm{}H} = 1@f$. To overcome this, we do a second
 * order Taylor expansion of the square root, which leads to
 * @f[
 *   x_{\rm{}H} = \frac{1}{C_{\rm{}H}}.
 * @f]
 *
 * @param alphaH Hydrogen recombination rate (in m^3s^-1).
 * @param jH Hydrogen intensity integral (in s^-1).
 * @param nH Hydrogen number density (in m^-3).
 * @return Neutral fraction of hydrogen.
 */
double IonizationStateCalculator::compute_ionization_state_hydrogen(
    const double alphaH, const double jH, const double nH, const double gammaH, const double old_xn, double ts) {

  double xn;

  



  if (jH ==0 && nH > 0) {
    xn = alphaH/(alphaH+gammaH);
  } else if (jH > 0 && nH > 0){
    //equation (4) from K+O (2020), in turn from somewhere else...
    double denom = jH + (2*alphaH + gammaH)*nH;
    const double in_root = std::pow(jH + gammaH*nH,2.0) + 4.0*jH*alphaH*nH;
    denom = denom + std::pow(in_root,0.5);
    xn = 2.0*alphaH*nH/denom;

  } else {
    xn = 0;
  }


  if (ts > 0.0) {
    double largest_change = ts*(alphaH*nH*(std::pow(1.- old_xn,2.0)) - old_xn*jH - gammaH*nH*old_xn*(1.-old_xn));
    if (largest_change > 0 && xn < old_xn) {
      // this is a problem...
      cmac_error("Numerical is net recombining, but xn is lower than last step.")
    }
    if (largest_change < 0 && xn > old_xn) {
      // similarly a problem, lets hope this doesnt get called?
      cmac_error("Numerical is net ionizing but xn is increasing from last step.")
    }
    if (largest_change > 0) {
      // we are recombining
      if ((xn-old_xn) > largest_change) {
        //we have over-recombined, add limiter by implementing numerical time dependence
        std::cout << "Applied limiter, old x =" << old_xn << " was gonna go " << xn << " but will go to " << old_xn + largest_change << " instead " << std::endl;
        xn = old_xn + largest_change;
        
      }
    } else if (largest_change < 0) {
      // we are ionizing 
      if ((xn - old_xn) < largest_change) {
        std::cout << "Applied limiter, old x =" << old_xn << " was gonna go " << xn << " but will go to " << old_xn + largest_change << " instead " << std::endl;
        // we have over-ionized for this time, implement limiter
        xn = old_xn + largest_change;
      }
    }
  }



  return std::max(1.e-14,xn);
}
