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
 * @file HomogeneousDensityFunction.hpp
 *
 * @brief Homogeneous DensityFunction implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HOMOGENEOUSDENSITYFUNCTION_HPP
#define HOMOGENEOUSDENSITYFUNCTION_HPP

#include "DensityFunction.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

/**
 * @brief DensityFunction that returns a constant value for all coordinates,
 * corresponding to a homogeneous density field.
 */
class HomogeneousDensityFunction : public DensityFunction {
private:
  /*! @brief Single density value for the entire box (in m^-3). */
  const double _density;

  /*! @brief Single temperature value for the entire box (in K). */
  const double _temperature;

  /*! @brief Initial hydrogen neutral fraction for the entire box. */
  const double _neutral_fraction_H;

  const double _dust_to_gas;

  const double _fraction_silicates;

public:
  /**
   * @brief Constructor.
   *
   * @param density Single density value for the entire box (in m^-3).
   * @param temperature Single temperature value for the entire box (in K).
   * @param neutral_fraction_H Single hydrogen neutral fraction value for the
   * entire box.
   * @param log Log to write logging information to.
   */
  HomogeneousDensityFunction(double density = 1., double temperature = 8000.,
                             double neutral_fraction_H = 1.e-6,
                             double dust_to_gas = 0.0, double fraction_silicates = 0.6,
                             Log *log = nullptr)
      : _density(density), _temperature(temperature),
        _neutral_fraction_H(neutral_fraction_H),
        _dust_to_gas(dust_to_gas),_fraction_silicates(fraction_silicates) {

    if (log) {
      log->write_status(
          "Created HomogeneousDensityFunction with constant density ", _density,
          " m^-3 and constant temperature ", _temperature, " K. Dust gas ratio of ",
           _dust_to_gas, " and fraction silicates ", _fraction_silicates);
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - density: Constant number density value (default: 100. cm^-3)
   *  - temperature: Constant initial temperature value (default: 8000. K)
   *  - neutral fraction H: Constant initial neutral fraction value
   *    (default: 1.e-6)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging information to.
   */
  HomogeneousDensityFunction(ParameterFile &params, Log *log = nullptr)
      : HomogeneousDensityFunction(
            params.get_physical_value< QUANTITY_NUMBER_DENSITY >(
                "DensityFunction:density", "100. cm^-3"),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "DensityFunction:temperature", "8000. K"),
            params.get_value< double >("DensityFunction:neutral fraction H",
                                       1.e-6),
            params.get_value<double>("DensityFunction:dust to gas",0.0),
            params.get_value<double>("DensityFunction:fraction silicates",0.6),
            log) {}

  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) {
    DensityValues values;
    values.set_number_density(_density);
    values.set_temperature(_temperature);
    values.set_ionic_fraction(ION_H_n, _neutral_fraction_H);
    values.set_dust_gas_ratio(_dust_to_gas);
    values.set_fraction_silicates(_fraction_silicates);
#ifdef HAS_HELIUM
    values.set_ionic_fraction(ION_He_n, 0.999);
    values.set_ionic_fraction(ION_He_p1, 0.001);

#endif
#ifdef HAS_CARBON
    values.set_ionic_fraction(ION_C_p1, 0.999);
    values.set_ionic_fraction(ION_C_p2, 0.001);
#endif
#ifdef HAS_NITROGEN
    values.set_ionic_fraction(ION_N_n, 0.999);
    values.set_ionic_fraction(ION_N_p1, 0.001);
    values.set_ionic_fraction(ION_N_p2,0.000);
#endif
#ifdef HAS_OXYGEN
    values.set_ionic_fraction(ION_O_n, 0.999);
    values.set_ionic_fraction(ION_O_p1, 0.001);
    values.set_ionic_fraction(ION_O_p2, 0.000);
    values.set_ionic_fraction(ION_O_p3, 0.000);
#endif
#ifdef HAS_NEON
    values.set_ionic_fraction(ION_Ne_n, 0.999);
    values.set_ionic_fraction(ION_Ne_p1, 0.001);
    values.set_ionic_fraction(ION_Ne_p2, 0.000);
    values.set_ionic_fraction(ION_Ne_p3, 0.000);
#endif
#ifdef HAS_SULPHUR
    values.set_ionic_fraction(ION_S_p1, 0.999);
    values.set_ionic_fraction(ION_S_p2, 0.001);
    values.set_ionic_fraction(ION_S_p3, 0.000);
#endif
    return values;
  }
};

#endif // HOMOGENEOUSDENSITYFUNCTION_HPP
