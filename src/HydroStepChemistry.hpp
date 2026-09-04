#ifndef HYDROSTEPCHEMISTRY_HPP
#define HYDROSTEPCHEMISTRY_HPP

#include "IonizationStateCalculator.hpp"
#include "PhysicalConstants.hpp"

namespace HydroStepChemistry {

/** Radiation iterations predict opacity over ONE hydro dt. They do not commit
 * chemistry. On the last iteration restore the physical starting state and
 * retain all photoionization/heating rates in SI, ready for hydro-step solves.
 * Use a trial copy because optional output code overwrites H/He estimators. */
inline void radiation_trial(IonizationVariables &cell,
    const IonizationStateCalculator &calculator, double jfac, double dt,
    bool last_iteration) {
  const double hfac = jfac * PhysicalConstants::get_physical_constant(
                                PHYSICALCONSTANT_PLANCK);
  IonizationVariables trial = cell;
  calculator.calculate_ionization_state(jfac, hfac, trial, dt, true, false);
  for (int ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
    cell.set_ionic_fraction(ion, last_iteration
        ? cell.get_prev_ionic_fraction(ion) : trial.get_ionic_fraction(ion));
    if (last_iteration) {
      cell.set_mean_intensity(ion, jfac * cell.get_mean_intensity(ion));
    }
  }
  if (last_iteration) {
    for (int h = 0; h < NUMBER_OF_HEATINGTERMS; ++h) {
      cell.set_heating(h, hfac * cell.get_heating(h));
    }
  }
}

/** Advance the post-advection state exactly once. Rates remain Eulerian and
 * fixed until the next photon shoot; density and temperature are current. */
inline void advance(IonizationVariables &cell,
    const IonizationStateCalculator &calculator, double dt) {
  cell.copy_previous_fractions();
  calculator.calculate_ionization_state(1., 1., cell, dt, true, true);
}
}
#endif
