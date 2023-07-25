/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file PowerLawPhotonSourceSpectrum.hpp
 *
 * @brief PhotonSourceSpectrum implementation for a power law spectrum.
 *
 * @author Lewis McCallum (lm261@st-andrews.ac.uk)
 */
#ifndef POWERLAWPHOTONSOURCESPECTRUM_HPP
#define POWERLAWPHOTONSOURCESPECTRUM_HPP

#include "Error.hpp"
#include "PhotonSourceSpectrum.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief PhotonSourceSpectrum implementation for a uniform spectrum.
 */
class PowerLawPhotonSourceSpectrum : public PhotonSourceSpectrum {
private:
    double _alpha;
    double _v1;
    double _v2;
public:
    PowerLawPhotonSourceSpectrum(double alpha, Log *log = nullptr):
    _alpha(alpha),_v1(13.6),_v2(54.0) {
     //go from energy power law to photon sampling spectrum
        _alpha += 1.0;
        if (_alpha == 1.0) {
            cmac_error("Choose a slightly different alpha for power law spectrum. 0 breaks the formula.")
        }
    }

    PowerLawPhotonSourceSpectrum(ParameterFile &params, Log *log = nullptr):
        PowerLawPhotonSourceSpectrum(
            params.get_value<double>("PhotonSourceSpectrum:alpha",2.0), log) {} 


   /**
   * @brief Virtual destructor.
   */
  virtual ~PowerLawPhotonSourceSpectrum() {}

  /**
   * @brief Get a random frequency from the spectrum.
   *
   * @param random_generator RandomGenerator to use.
   * @param temperature Temperature of the gas (for reemission spectra) (in K).
   * @return Random frequency, distributed according to the spectrum (in Hz).
   */
  virtual double get_random_frequency(RandomGenerator &random_generator,
                                      double temperature = 0.) const {
       double ep = random_generator.get_uniform_random_double();
       double vr = std::pow(_v1,1.0-_alpha) + ep*(std::pow(_v2,1.0-_alpha) - std::pow(_v1,1.0-_alpha));
       vr = std::pow(vr,1./(1.-_alpha));
    return vr*3.289e15;
  }

  /**
   * @brief Get the total ionizing flux emitted by the spectrum.
   *
   * @return Total ionizing flux (in m^-2 s^-1).
   */
  virtual double get_total_flux() const {
    cmac_error("This function should not be used!");
    return 0.;
  }
};

#endif // POWERLAWPHOTONSOURCESPECTRUM_HPP
