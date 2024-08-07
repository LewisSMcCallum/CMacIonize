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
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>


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
    IonizationVariables &ionization_variables, double timestep, bool time_dependent, bool do_metals) const {



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
      if (time_dependent) {
      h0 = ionization_variables.get_prev_ionic_fraction(ION_H_n);
      he0 = ionization_variables.get_prev_ionic_fraction(ION_He_n);
      hep = ionization_variables.get_prev_ionic_fraction(ION_He_p1);
      compute_time_dependent_hydrogen_helium(alphaH, alphaHe,alphaHe2, jH, jHe, ntot,
                                                AHe, T, h0, he0, hep, gammaH, gammaHe1,
                                                 gammaHe2,timestep);

      } else {
      compute_ionization_states_hydrogen_helium(alphaH, alphaHe,alphaHe2, jH, jHe, ntot,
                                                AHe, T, h0, he0, hep, gammaH, gammaHe1,
                                                 gammaHe2);
      }


    } else {
      if (time_dependent) {
        h0 = compute_time_dependent_hydrogen(alphaH, jH, ntot, gammaH, ionization_variables.get_prev_ionic_fraction(ION_H_n), timestep);
      } else {
         h0 = compute_ionization_state_hydrogen(alphaH, jH, ntot, gammaH, ionization_variables.get_prev_ionic_fraction(ION_H_n), timestep);
      }
    }
#else
  double h0;
    if (time_dependent) {
      h0 = compute_time_dependent_hydrogen(alphaH, jH, ntot, gammaH, ionization_variables.get_prev_ionic_fraction(ION_H_n), timestep);
    } else {
      h0 = compute_ionization_state_hydrogen(alphaH, jH, ntot, gammaH, ionization_variables.get_prev_ionic_fraction(ION_H_n), timestep);
    }
#endif

    ionization_variables.set_ionic_fraction(ION_H_n, std::min(1.0, std::max(h0, 1e-14)));

#ifdef HAS_HELIUM
    ionization_variables.set_ionic_fraction(ION_He_n, std::min(1.0, std::max(he0, 1e-14)));
    ionization_variables.set_ionic_fraction(ION_He_p1,std::min(1.0, std::max(hep, 1e-14)));
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
  if (do_metals) {
      if (time_dependent) {
      compute_time_dependent_metals(
          j_metals, ne, T, T4, nh0, nhe0, nhp, _recombination_rates,
          _charge_transfer_rates, _collisional_rates, ionization_variables, timestep);

      } else {
      compute_ionization_states_metals(
          j_metals, ne, T, T4, nh0, nhe0, nhp, _recombination_rates,
          _charge_transfer_rates, _collisional_rates, ionization_variables);

      }
  }


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
  ionization_variables.set_ionic_fraction(ION_C_p1, 1.0 - (C21 * sumC_inv) - (C31 * sumC_inv));
  ionization_variables.set_ionic_fraction(ION_C_p2, C21 * sumC_inv);
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
  ionization_variables.set_ionic_fraction(ION_N_n, 1 - (N21 * sumN_inv) - (N31 * sumN_inv) - (N41 * sumN_inv));
  ionization_variables.set_ionic_fraction(ION_N_p1, N21 * sumN_inv);
  ionization_variables.set_ionic_fraction(ION_N_p2, N31 * sumN_inv);
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



  ionization_variables.set_ionic_fraction(ION_O_n, 1.0 - (O21 * sumO_inv) - (O31 * sumO_inv) - (O41 * sumO_inv) - (O51 * sumO_inv));
  ionization_variables.set_ionic_fraction(ION_O_p1, O21 * sumO_inv);
  ionization_variables.set_ionic_fraction(ION_O_p2, O31 * sumO_inv);
  ionization_variables.set_ionic_fraction(ION_O_p3, O41 * sumO_inv);


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

  ionization_variables.set_ionic_fraction(ION_Ne_n, 1.0 - (Ne21 * sumNe_inv) - (Ne31 * sumNe_inv) - (Ne41 * sumNe_inv) - (Ne51 * sumNe_inv));
  ionization_variables.set_ionic_fraction(ION_Ne_p1, Ne21 * sumNe_inv);
  ionization_variables.set_ionic_fraction(ION_Ne_p2, Ne31 * sumNe_inv);
  ionization_variables.set_ionic_fraction(ION_Ne_p3, Ne41 *sumNe_inv);
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
  ionization_variables.set_ionic_fraction(ION_S_p1, 1.0 - (S21 * sumS_inv) - (S31 * sumS_inv) - (S41 * sumS_inv));
  ionization_variables.set_ionic_fraction(ION_S_p2, S21 * sumS_inv);
  ionization_variables.set_ionic_fraction(ION_S_p3, S31 * sumS_inv);
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
    const double totweight, DensitySubGrid &subgrid, double timestep, bool time_dependent, bool do_metals) const {

  
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
                               cellit.get_ionization_variables(),timestep, time_dependent, do_metals);
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


  if ((jH > 0 ) && (jHe ==0) && (gammaHe1 < 1e-25)) {
    he0 = 0.9999;
    hep = 0.0001;
    double hepp = 0.0;

    double h0old = 0.8;
    h0 = 0.9*h0old;
  while (std::abs(h0 - h0old) > 1.e-4 * h0old) {
    h0old = h0;

    double ne = nH*(1 - h0) + 2*hepp*AHe*nH + hep*AHe*nH;
    // calculate a new guess for C_H
    double pHots = 1. / (1. + 77. * he0 / std::sqrt(T) / h0old);

    double alpha_e_2sP = 4.17e-20 * std::pow(T * 1.e-4, -0.861);

    double ots = AHe*ne*hep*pHots*alpha_e_2sP;

    h0 = (ne*alphaH - ots)/(jH + ne*alphaH + ne*gammaH);
  }
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

    double hepp = (1.0-he0)*gammaHe2/(alphaHe2+gammaHe2);

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
    if (niter > 200) {
      if ((h0 < 1e-6)) {
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
  double xrun = old_xn;


  
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

  uint_fast32_t divs = 10;

  

  
  if ((ts > 0.0) & (old_xn > -0.5)) {
    for (uint_fast32_t i=0;i<divs;i++) {
      double largest_change = ts*(alphaH*nH*(std::pow(1.- xrun,2.0)) - xrun*jH - gammaH*nH*xrun*(1.-xrun))/divs;
      // if (largest_change > 0 && xn < xrun) {
      //   // this is a problem...
      //   cmac_error("Numerical is net recombining, but xn is lower than last step.")
      // }
      // if (largest_change < 0 && xn > xrun) {
      //   // similarly a problem, lets hope this doesnt get called?
      //   cmac_error("Numerical is net ionizing but xn is increasing from last step.")
      // }
      if (largest_change > 0) {
        // we are recombining
        if ((xn-xrun) > largest_change) {
          //we have over-recombined, add limiter by implementing numerical time dependence
          //std::cout << "Applied limiter, old x =" << old_xn << " was gonna go " << xn << " but will go to " << old_xn + largest_change << " instead " << std::endl;
          xrun = xrun + largest_change;
          
        } else {
          xrun = xn;
          break;
        }
      } else if (largest_change < 0) {
        // we are ionizing 
        if ((xn - xrun) < largest_change) {
          //std::cout << "Applied limiter, old x =" << old_xn << " was gonna go " << xn << " but will go to " << old_xn + largest_change << " instead " << std::endl;
          // we have over-ionized for this time, implement limiter
          xrun = xrun + largest_change;
        } else {
          xrun = xn;
          break;
        }
      }
  }
  xn = xrun;
  }



  return std::max(1.e-14,xn);
}
int hydrogen_ode_system(double t, const double y[], double f[], void *params) {
    (void)(t); // Avoid unused parameter warning

    // Extract rate coefficients and total density from params
    double *coefficients = static_cast<double*>(params);
    double k_coll = coefficients[0];
    double k_photo = coefficients[1];
    double k_rec = coefficients[2];
    double n_total = coefficients[3];

    // Neutral fraction
    double x = y[0];

    // ODE for the neutral fraction x
    f[0] = -k_coll * x * (1 - x) * n_total - k_photo * x + k_rec * (1 - x) * (1 - x) * n_total;

    return GSL_SUCCESS;
}

double IonizationStateCalculator::compute_time_dependent_hydrogen(
    const double alphaH, const double jH, const double nH, const double gammaH, const double old_xn, double ts) {

  double coefficients[4] = {gammaH, jH, alphaH, nH};
  // Initial conditions: n_H, n_H_plus
  double y[1] = {old_xn}; // Example initial population densities

    // Time domain
  double t = 0.0;

    // Set up the solver
  gsl_odeiv2_system sys = {hydrogen_ode_system, nullptr, 1, coefficients};

  gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rkf45, ts/100., 1e-4, 0.0);


  int status = gsl_odeiv2_driver_apply(driver, &t, ts, y);


  if (status != GSL_SUCCESS) {
      cmac_warning("Error in solver! xn = %g",y[0]);
     // xn = 1e-14;
  }
    

  gsl_odeiv2_driver_free(driver);

  double xn = y[0];



  return std::max(1.e-14,xn);
}

int hydrogen_helium_ode_system(double t, const double y[], double f[], void *params) {
    (void)(t); // Avoid unused parameter warning

    // Extract rate coefficients and total density from params
    double *coefficients = static_cast<double*>(params);
    double alphaH = coefficients[0];
    double alphaHe = coefficients[1];
    double alphaHe2 = coefficients[2];
    double jH = coefficients[3];
    double jHe = coefficients[4];
    double nH = coefficients[5];
    double AHe = coefficients[6];
    double gammaH = coefficients[7];
    double gammaHe1 = coefficients[8];
    double gammaHe2 = coefficients[9];
    double sqrtT = coefficients[10];
    double alpha_e_2sP = coefficients[11];
    

    // Neutral fraction
    double xh = y[0];
    double he0 = y[1];
    double hep = y[2];

    double hepp = 1-he0-hep;

    double ne = (1-xh)*nH + hep*nH*AHe + 2.0*hepp*nH*AHe;

    double pHots = 1. / (1. + 77. * he0 / sqrtT / xh);


    // ODE for the neutral fraction x
    f[0] = -gammaH*ne*xh - jH*xh + alphaH*(1-xh)*ne - ne*pHots*hep*alpha_e_2sP*AHe;
    f[1] = -gammaHe1*ne*he0 - jHe*he0 + alphaHe*ne*hep;
    f[2] = gammaHe1*ne*he0 + jHe*he0 + alphaHe2*ne*hepp - alphaHe*ne*hep - gammaHe2*ne*hep;

    for (int i = 0; i < 3; ++i) {
    // Preventing negative growth for negative values
      if (y[i] < 1e-14) {
          f[i] = std::max(f[i], 0.0);
      } 
      // Preventing positive growth for values greater than 1
      else if (y[i] > 1.0) {
              f[i] = std::min(f[i], 0.0);
      }
    }

    return GSL_SUCCESS;
}


void IonizationStateCalculator::compute_time_dependent_hydrogen_helium(
    const double alphaH, const double alphaHe, const double alphaHe2, const double jH,
    const double jHe, const double nH, const double AHe, const double T,
    double &h0, double &he0, double &hep, const double gammaH, const double gammaHe1,
    const double gammaHe2, double ts) {


  double alpha_e_2sP = 4.17e-20 * std::pow(T * 1.e-4, -0.861);
  double sqrtT = std::sqrt(T);


  double coefficients[12] = {alphaH, alphaHe, alphaHe2, jH, jHe, nH, AHe, gammaH, gammaHe1, gammaHe2, 
              sqrtT,alpha_e_2sP};

  double y[3] = {h0,he0,hep};

  double t = 0.0;

  gsl_odeiv2_system sys = {hydrogen_helium_ode_system, nullptr, 3, coefficients};

  gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rkf45, ts/100., 1e-4, 0.0);

  int status = gsl_odeiv2_driver_apply(driver, &t, ts, y);



  if (status != GSL_SUCCESS) {
      std::cout << h0 << " " << he0 << " " << hep << std::endl;
      std::cout <<  alphaH << " " << alphaHe <<  " " << alphaHe2 << " " << jH << " "  << jHe << " "  << gammaH << " "  << gammaHe1 << " "  << gammaHe2 << std::endl;
      std::cout << nH << " "  << AHe << " "  << T << " "  << std::endl;
      cmac_error("Error in solver!");
  }
    

  gsl_odeiv2_driver_free(driver);

  h0 = std::max(y[0],1e-14);
  he0 = std::max(y[1],1e-14);
  hep = std::max(y[2],1e-14);


  h0 = std::min(h0,1.);
  he0 = std::min(he0,1.);
  hep = std::min(hep,1.);

  if (h0 != h0) {
    cmac_warning("Nan H0 in solver. Setting 0.999");
    h0 = 0.999;
  }
  if (he0 != he0) {
    cmac_warning("Nan He0 in solver. Setting 0.999");
    he0 = 0.999;
  }
  if (hep != hep) {
    cmac_warning("Nan Hep in solver. Setting 0.001");
    hep = 0.001;
  }


}

struct ODEParams {
    std::vector<std::vector<double>> coefficients;
    double ne;
};


inline int metals_ode_system(double t, const double y[], double f[], void *params) {
    (void)(t); // Avoid unused parameter warning
    ODEParams* p = static_cast<ODEParams*>(params);
    std::vector<std::vector<double>>& coefficients = p->coefficients;
    double& ne = p->ne;
    //levels here is one less than total number of states
    size_t levels = coefficients.size();
    // Boundary conditions

    double frac_last = 1.0;
    for (size_t i = 0; i < coefficients.size(); ++i) {
        frac_last -= y[i];
    }

    
//change for first level, and last level, note last level is actually second last level, with the highest level not being tracked explicity
    f[0] = -coefficients[0][0]*y[0]*ne - coefficients[0][1]*y[0] + coefficients[0][2]*y[1]*ne
         -coefficients[0][3]*y[0] + coefficients[0][4]*y[1];
    f[levels-1] = -coefficients[levels-1][0]*y[levels-1]*ne - coefficients[levels-1][1]*y[levels-1] + coefficients[levels-1][2]*frac_last*ne
        +coefficients[levels-2][0]*y[levels-2]*ne + coefficients[levels-2][1]*y[levels-2] - coefficients[levels-2][2]*y[levels-1]*ne
        -coefficients[levels-1][3]*y[levels-1] + coefficients[levels-1][4]*frac_last + coefficients[levels-2][3]*y[levels-2] - coefficients[levels-2][4]*y[levels-1];


    // Middle levels
    for (size_t i = 1; i < levels-1; ++i) {
//ok this is terrible but lets try and comment this, in order the terms are
// loss due to collisions up, loss due to photoionize up, gain due to recombine from above
//  gain due to collisions from below, gain from photoionize from below, loss as recombine down
// loss from CT ionize up, gain from CT recombine from above, gain from CT ionize from below, loss from CT recombine down
      f[i] = -coefficients[i][0]*y[i]*ne - coefficients[i][1]*y[i] + coefficients[i][2]*y[i+1]*ne
        +coefficients[i-1][0]*y[i-1]*ne + coefficients[i-1][1]*y[i-1] - coefficients[i-1][2]*y[i]*ne
        -coefficients[i][3]*y[i] + coefficients[i][4]*y[i+1] + coefficients[i-1][3]*y[i-1] - coefficients[i-1][4]*y[i];
    }

    for (size_t i = 0; i < levels; ++i) {
    // Preventing negative growth for negative values
      if (y[i] < 1e-14) {
          f[i] = std::max(f[i], 0.0);
      } 
      // Preventing positive growth for values greater than 1
      else if (y[i] > 1.0) {
              f[i] = std::min(f[i], 0.0);
      }
    }

    return GSL_SUCCESS;
}



void IonizationStateCalculator::compute_time_dependent_metals(
    const double *j_metals, const double ne, const double T, const double T4,
    const double nh0, const double nhe0, const double nhp,
    const RecombinationRates &recombination_rates,
    const ChargeTransferRates &charge_transfer_rates,
    const CollisionalRates &collisional_rates,
    IonizationVariables &ionization_variables, double ts) {

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


//CARBON
#ifdef HAS_CARBON
{
      const size_t levels_carbon = 3;

      double y[levels_carbon-1] = {ionization_variables.get_ionic_fraction(ION_C_p1), 
                              ionization_variables.get_ionic_fraction(ION_C_p2)};

      ODEParams params;
      params.coefficients = std::vector<std::vector<double>>(levels_carbon-1, std::vector<double>(5, 0.0));
      params.ne = ne;

    //set collisional rates 
      params.coefficients[0][0] = collisional_rates.get_collisional_rate(ION_C_p1, T);
      params.coefficients[1][0] = collisional_rates.get_collisional_rate(ION_C_p2, T);

    //set photoionization rates
      params.coefficients[0][1] = jCp1;
      params.coefficients[1][1] = jCp2;

    //set recombination rates, note rate is rate of tranition which MAKES named ion, i.e. ION_H_n rate is times by H+

      params.coefficients[0][2] = recombination_rates.get_recombination_rate(ION_C_p1, T);
      params.coefficients[1][2] = recombination_rates.get_recombination_rate(ION_C_p2, T);

    // set charge transfer rates, multiply them by appropriate densities here, use as pure rate in solver

    //index [3] total ionization rates, index [4] total recomb rates
     
    // for C, level 1 is neglible, only have recomb for level 2
      
      params.coefficients[1][4] = nh0*charge_transfer_rates.get_charge_transfer_recombination_rate_H(ION_C_p2,T4)
                        + nhe0*charge_transfer_rates.get_charge_transfer_recombination_rate_He(ION_C_p2,T4);


      double t = 0.0;

      gsl_odeiv2_system sys = {metals_ode_system, nullptr, levels_carbon - 1, &params};
      gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
          &sys, gsl_odeiv2_step_rkf45, ts/100., 1.e-4, 0.0);

      int status = gsl_odeiv2_driver_apply(driver, &t, ts, y);

      if (status != GSL_SUCCESS) {
          cmac_error("Error in solver!");
      }
    
      gsl_odeiv2_driver_free(driver);
      //set new 
      ionization_variables.set_ionic_fraction(ION_C_p1, std::min(1.0,std::max(y[0],1e-14)));
      ionization_variables.set_ionic_fraction(ION_C_p2, std::min(1.0,std::max(y[1],1e-14)));
}    
#endif

#ifdef HAS_NITROGEN
{
      const size_t levels_nitrogen = 4;

      double y[levels_nitrogen-1] = {ionization_variables.get_ionic_fraction(ION_N_n), 
                              ionization_variables.get_ionic_fraction(ION_N_p1),
                              ionization_variables.get_ionic_fraction(ION_N_p2)};

      ODEParams params;
      params.coefficients = std::vector<std::vector<double>>(levels_nitrogen-1, std::vector<double>(5, 0.0));
      params.ne = ne;

    //set collisional rates 
      params.coefficients[0][0] = collisional_rates.get_collisional_rate(ION_N_n, T);
      params.coefficients[1][0] = collisional_rates.get_collisional_rate(ION_N_p1, T);
      params.coefficients[2][0] = collisional_rates.get_collisional_rate(ION_N_p2, T);

    //set photoionization rates
      params.coefficients[0][1] = jNn;
      params.coefficients[1][1] = jNp1;
      params.coefficients[2][1] = jNp2;

    //set recombination rates, note rate is rate of tranition which MAKES named ion, i.e. ION_H_n rate is times by H+

      params.coefficients[0][2] = recombination_rates.get_recombination_rate(ION_N_n, T);
      params.coefficients[1][2] = recombination_rates.get_recombination_rate(ION_N_p1, T);
      params.coefficients[2][2] = recombination_rates.get_recombination_rate(ION_N_p2, T);

    // set charge transfer rates, multiply them by appropriate densities here, use as pure rate in solver

    //index [3] total ionization rates, index [4] total recomb rates
     
      params.coefficients[0][3] = nhp*charge_transfer_rates.get_charge_transfer_ionization_rate_H(ION_N_n,T4);
      params.coefficients[0][4] = nh0*charge_transfer_rates.get_charge_transfer_recombination_rate_H(ION_N_n,T4);

      params.coefficients[1][4] = nh0*charge_transfer_rates.get_charge_transfer_recombination_rate_H(ION_N_p1,T4)
                        + nhe0*charge_transfer_rates.get_charge_transfer_recombination_rate_He(ION_N_p1,T4);

      params.coefficients[2][4] = nh0*charge_transfer_rates.get_charge_transfer_recombination_rate_H(ION_N_p2,T4)
                        + nhe0*charge_transfer_rates.get_charge_transfer_recombination_rate_He(ION_N_p2,T4);

      double t = 0.0;

      gsl_odeiv2_system sys = {metals_ode_system, nullptr, levels_nitrogen - 1, &params};
      gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
          &sys, gsl_odeiv2_step_rkf45, ts/100., 1e-4, 0.0);

      int status = gsl_odeiv2_driver_apply(driver, &t, ts, y);

      if (status != GSL_SUCCESS) {
          cmac_error("Error in solver!");
      }
    
      gsl_odeiv2_driver_free(driver);
      //set new 
      ionization_variables.set_ionic_fraction(ION_N_n, std::min(1.0,std::max(y[0],1e-14)));
      ionization_variables.set_ionic_fraction(ION_N_p1, std::min(1.0,std::max(y[1],1e-14)));
      ionization_variables.set_ionic_fraction(ION_N_p2, std::min(1.0,std::max(y[2],1e-14)));
}
#endif

#ifdef HAS_OXYGEN
{
      const size_t levels_oxygen = 5;

      double y[levels_oxygen-1] = {ionization_variables.get_ionic_fraction(ION_O_n), 
                              ionization_variables.get_ionic_fraction(ION_O_p1),
                              ionization_variables.get_ionic_fraction(ION_O_p2),
                              ionization_variables.get_ionic_fraction(ION_O_p3)};

      ODEParams params;
      params.coefficients = std::vector<std::vector<double>>(levels_oxygen-1, std::vector<double>(5, 0.0));
      params.ne = ne;

    //set collisional rates 
      params.coefficients[0][0] = collisional_rates.get_collisional_rate(ION_O_n, T);
      params.coefficients[1][0] = collisional_rates.get_collisional_rate(ION_O_p1, T);
      params.coefficients[2][0] = collisional_rates.get_collisional_rate(ION_O_p2, T);
      params.coefficients[3][0] = collisional_rates.get_collisional_rate(ION_O_p3, T);

    //set photoionization rates
      params.coefficients[0][1] = jOn;
      params.coefficients[1][1] = jOp1;

    //set recombination rates, note rate is rate of tranition which MAKES named ion, i.e. ION_H_n rate is times by H+

      params.coefficients[0][2] = recombination_rates.get_recombination_rate(ION_O_n, T);
      params.coefficients[1][2] = recombination_rates.get_recombination_rate(ION_O_p1, T);
      params.coefficients[2][2] = recombination_rates.get_recombination_rate(ION_O_p2, T);
      params.coefficients[3][2] = recombination_rates.get_recombination_rate(ION_O_p3, T);

    // set charge transfer rates, multiply them by appropriate densities here, use as pure rate in solver

    //index [3] total ionization rates, index [4] total recomb rates
     
      params.coefficients[0][3] = nhp*charge_transfer_rates.get_charge_transfer_ionization_rate_H(ION_O_n,T4);
      params.coefficients[0][4] = nh0*charge_transfer_rates.get_charge_transfer_recombination_rate_H(ION_O_n,T4);

      params.coefficients[1][4] = nh0*charge_transfer_rates.get_charge_transfer_recombination_rate_H(ION_O_p1,T4)
                        + nhe0*charge_transfer_rates.get_charge_transfer_recombination_rate_He(ION_O_p1,T4);

      params.coefficients[2][4] = nh0*charge_transfer_rates.get_charge_transfer_recombination_rate_H(ION_O_p2,T4)
                        + nhe0*charge_transfer_rates.get_charge_transfer_recombination_rate_He(ION_O_p2,T4);

      params.coefficients[3][4] = nh0*charge_transfer_rates.get_charge_transfer_recombination_rate_H(ION_O_p3,T4)
                        + nhe0*charge_transfer_rates.get_charge_transfer_recombination_rate_He(ION_O_p3,T4);

      double t = 0.0;

      gsl_odeiv2_system sys = {metals_ode_system, nullptr, levels_oxygen - 1, &params};
      gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
          &sys, gsl_odeiv2_step_rkf45, ts/100., 1e-4, 0.0);

      int status = gsl_odeiv2_driver_apply(driver, &t, ts, y);

      if (status != GSL_SUCCESS) {
          cmac_error("Error in solver!");
      }
    
      gsl_odeiv2_driver_free(driver);
      //set new 
      ionization_variables.set_ionic_fraction(ION_O_n, std::min(1.0,std::max(y[0],1e-14)));
      ionization_variables.set_ionic_fraction(ION_O_p1, std::min(1.0,std::max(y[1],1e-14)));
      ionization_variables.set_ionic_fraction(ION_O_p2, std::min(1.0,std::max(y[2],1e-14)));
      ionization_variables.set_ionic_fraction(ION_O_p3, std::min(1.0,std::max(y[3],1e-14)));
}
#endif

#ifdef HAS_NEON
{
      const size_t levels_neon = 5;

      double y[levels_neon-1] = {ionization_variables.get_ionic_fraction(ION_Ne_n), 
                              ionization_variables.get_ionic_fraction(ION_Ne_p1),
                              ionization_variables.get_ionic_fraction(ION_Ne_p2),
                              ionization_variables.get_ionic_fraction(ION_Ne_p3)};

      ODEParams params;
      params.coefficients = std::vector<std::vector<double>>(levels_neon-1, std::vector<double>(5, 0.0));
      params.ne = ne;

    //set collisional rates 
      params.coefficients[0][0] = collisional_rates.get_collisional_rate(ION_Ne_n, T);
      params.coefficients[1][0] = collisional_rates.get_collisional_rate(ION_Ne_p1, T);
      params.coefficients[2][0] = collisional_rates.get_collisional_rate(ION_Ne_p2, T);
      params.coefficients[3][0] = collisional_rates.get_collisional_rate(ION_Ne_p3, T);

    //set photoionization rates
      params.coefficients[0][1] = jNen;
      params.coefficients[1][1] = jNep1;

    //set recombination rates, note rate is rate of tranition which MAKES named ion, i.e. ION_H_n rate is times by H+

      params.coefficients[0][2] = recombination_rates.get_recombination_rate(ION_Ne_n, T);
      params.coefficients[1][2] = recombination_rates.get_recombination_rate(ION_Ne_p1, T);
      params.coefficients[2][2] = recombination_rates.get_recombination_rate(ION_Ne_p2, T);
      params.coefficients[3][2] = recombination_rates.get_recombination_rate(ION_Ne_p3, T);

    // set charge transfer rates, multiply them by appropriate densities here, use as pure rate in solver

    //index [3] total ionization rates, index [4] total recomb rates
     
      params.coefficients[1][4] = nh0*charge_transfer_rates.get_charge_transfer_recombination_rate_H(ION_Ne_p1,T4)
                        + nhe0*charge_transfer_rates.get_charge_transfer_recombination_rate_He(ION_Ne_p1,T4);

      params.coefficients[2][4] = nh0*charge_transfer_rates.get_charge_transfer_recombination_rate_H(ION_Ne_p2,T4)
                        + nhe0*charge_transfer_rates.get_charge_transfer_recombination_rate_He(ION_Ne_p2,T4);


      double t = 0.0;

      gsl_odeiv2_system sys = {metals_ode_system, nullptr, levels_neon - 1, &params};
      gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
          &sys, gsl_odeiv2_step_rkf45, ts/100., 1e-4, 0.0);

      int status = gsl_odeiv2_driver_apply(driver, &t, ts, y);

      if (status != GSL_SUCCESS) {
          cmac_error("Error in solver!");
      }
    
      gsl_odeiv2_driver_free(driver);
      //set new 
      ionization_variables.set_ionic_fraction(ION_Ne_n, std::min(1.0,std::max(y[0],1e-14)));
      ionization_variables.set_ionic_fraction(ION_Ne_p1, std::min(1.0,std::max(y[1],1e-14)));
      ionization_variables.set_ionic_fraction(ION_Ne_p2, std::min(1.0,std::max(y[2],1e-14)));
      ionization_variables.set_ionic_fraction(ION_Ne_p3, std::min(1.0,std::max(y[3],1e-14)));
}
#endif

#ifdef HAS_SULPHUR
{
      const size_t levels_sulphur = 4;

      double y[levels_sulphur-1] = {ionization_variables.get_ionic_fraction(ION_S_p1), 
                              ionization_variables.get_ionic_fraction(ION_S_p2),
                              ionization_variables.get_ionic_fraction(ION_S_p3)};

      ODEParams params;
      params.coefficients = std::vector<std::vector<double>>(levels_sulphur-1, std::vector<double>(5, 0.0));
      params.ne = ne;

    //set collisional rates 
      params.coefficients[0][0] = collisional_rates.get_collisional_rate(ION_S_p1, T);
      params.coefficients[1][0] = collisional_rates.get_collisional_rate(ION_S_p2, T);
      params.coefficients[2][0] = collisional_rates.get_collisional_rate(ION_S_p3, T);

    //set photoionization rates
      params.coefficients[0][1] = jSp1;
      params.coefficients[1][1] = jSp2;
      params.coefficients[2][1] = jSp3;

    //set recombination rates, note rate is rate of tranition which MAKES named ion, i.e. ION_H_n rate is times by H+

      params.coefficients[0][2] = recombination_rates.get_recombination_rate(ION_S_p1, T);
      params.coefficients[1][2] = recombination_rates.get_recombination_rate(ION_S_p2, T);
      params.coefficients[2][2] = recombination_rates.get_recombination_rate(ION_S_p3, T);

    // set charge transfer rates, multiply them by appropriate densities here, use as pure rate in solver

    //index [3] total ionization rates, index [4] total recomb rates
     
    // for C, level 1 is neglible, only have recomb for level 2

      params.coefficients[0][4] = nh0*charge_transfer_rates.get_charge_transfer_recombination_rate_H(ION_S_p1,T4);
      
      params.coefficients[1][4] = nh0*charge_transfer_rates.get_charge_transfer_recombination_rate_H(ION_S_p2,T4)
                        + nhe0*charge_transfer_rates.get_charge_transfer_recombination_rate_He(ION_S_p2,T4);

      params.coefficients[2][4] = nh0*charge_transfer_rates.get_charge_transfer_recombination_rate_H(ION_S_p3,T4)
                        + nhe0*charge_transfer_rates.get_charge_transfer_recombination_rate_He(ION_S_p3,T4);


      double t = 0.0;

      gsl_odeiv2_system sys = {metals_ode_system, nullptr, levels_sulphur - 1, &params};
      gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
          &sys, gsl_odeiv2_step_rkf45, ts/100., 1e-4, 0.0);

      int status = gsl_odeiv2_driver_apply(driver, &t, ts, y);

      if (status != GSL_SUCCESS) {
          cmac_error("Error in solver!");
      }
    
      gsl_odeiv2_driver_free(driver);
      //set new 
      ionization_variables.set_ionic_fraction(ION_S_p1, std::min(1.0,std::max(y[0],1e-14)));
      ionization_variables.set_ionic_fraction(ION_S_p2, std::min(1.0,std::max(y[1],1e-14)));
      ionization_variables.set_ionic_fraction(ION_S_p3, std::min(1.0,std::max(y[2],1e-14)));

}
#endif


    
    
  }







