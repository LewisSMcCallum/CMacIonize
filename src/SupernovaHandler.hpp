/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file SupernovaHandler.hpp
 *
 * @brief Supernova handling functions.
 *
 * @author Lewis McCallum (lm261@st-andrews.ac.uk)
 */
#ifndef SUPERNOVAHANDLER_HPP
#define SUPERNOVAHANDLER_HPP

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <tuple>
#include <vector>


/**
 * @brief Handler with useful functions.
 */
class SupernovaHandler {
private:

  const double _sne_energy;
  bool _use_tigress_like_injection;

  static constexpr double parsec() { return 3.085677581491367e16; }
  static constexpr double solar_mass() { return 1.98847e30; }
  static constexpr double proton_mass() { return 1.67262192369e-27; }
  static constexpr double injection_radius_in_cells() { return 3.; }
  static constexpr double resolved_mass_fraction() { return 0.1; }
  static constexpr double thermal_energy_fraction() { return 0.72; }
  static constexpr double kinetic_energy_fraction() { return 0.28; }

  inline double get_feedback_mass(const double r_inj,
                                  const double nbar) const {
    const double volume = 4. * M_PI * r_inj * r_inj * r_inj / 3.;
    return nbar * 1.e6 * proton_mass() * volume;
  }

  /**
   * @brief Shell-formation mass used by the TIGRESS++ feedback criterion.
   *
   * Kim & Ostriker calibrate this for a two-phase medium as
   * M_sf = 1344 Msun E_51^0.87 n_H^-0.26.
   */
  inline double get_shell_formation_mass(const double nbar) const {
    if (!(nbar > 0.)) {
      return std::numeric_limits< double >::infinity();
    }
    const double energy_51 = _sne_energy / 1.e44;
    return 1344. * solar_mass() * std::pow(energy_51, 0.87) *
           std::pow(nbar, -0.26);
  }

  inline bool use_momentum_injection(const double r_inj,
                                     const double r_st, const double dx,
                                     const double nbar) const {
    if (!_use_tigress_like_injection) {
      return r_st < 4. * dx;
    }
    const double feedback_mass = get_feedback_mass(r_inj, nbar);
    const double shell_formation_mass = get_shell_formation_mass(nbar);
    return feedback_mass >=
           resolved_mass_fraction() * shell_formation_mass;
  }

  /**
   * @brief Report the injection branch selected by the resolution criterion.
   *
   * This is called once when the injection radius is measured, rather than
   * from inject_sne(), which is called once for every subgrid.
   */
  inline void report_injection_mode(const double r_inj, const double r_st,
                                    const double dx, const double nbar,
                                    const size_t num_cells) const {
    const double feedback_mass = get_feedback_mass(r_inj, nbar);
    const double shell_formation_mass = get_shell_formation_mass(nbar);
    const bool use_momentum =
        use_momentum_injection(r_inj, r_st, dx, nbar);
    std::cout << "SN_INJECTION_MODE mode="
              << (use_momentum ? "momentum" : "thermal")
              << " scheme="
              << (_use_tigress_like_injection ? "TIGRESS_like"
                                               : "SILCC_like")
              << " r_inj_pc=" << r_inj / parsec()
              << " r_sf_pc=" << r_st / parsec()
              << " dx_pc=" << dx / parsec()
              << " r_st_over_dx=" << r_st / dx
              << " m_fb_msun=" << feedback_mass / solar_mass()
              << " m_sf_msun=" << shell_formation_mass / solar_mass()
              << " m_fb_over_m_sf="
              << feedback_mass / shell_formation_mass
              << " criterion="
              << (_use_tigress_like_injection ? "m_fb_ge_0.1_m_sf"
                                               : "r_sf_lt_4dx")
              << " thermal_fraction="
              << (_use_tigress_like_injection ? thermal_energy_fraction() : 1.)
              << " kinetic_fraction="
              << (_use_tigress_like_injection ? kinetic_energy_fraction() : 0.)
              << " mass_return_msun=0"
              << " nbar_cm^-3=" << nbar << " cells=" << num_cells
              << std::endl;
  }



public:
  /**
   * @brief Constructor.
   *
   * @param sne_energy Energy per supernova (in J).
   * @param use_tigress_like_injection True for the TIGRESS-like prescription,
   * false for the legacy SILCC-like prescription.
   */
  SupernovaHandler(const double sne_energy,
                   const bool use_tigress_like_injection = true)
      : _sne_energy(sne_energy),
        _use_tigress_like_injection(use_tigress_like_injection) {

        //constructor functionality


  }

  inline std::tuple<double,double,double,double> get_r_inj(DensitySubGridCreator< HydroDensitySubGrid > *grid_creator,
                                         CoordinateVector<> sne_loc) {
      HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(sne_loc);
      double cell_vol =  subgrid.get_cell(sne_loc).get_volume();
      double dx = std::pow(cell_vol,1./3.);
      double r_run = _use_tigress_like_injection
                         ? injection_radius_in_cells() * dx
                         : 4. * dx;
      std::vector<std::pair<uint_fast32_t,uint_fast32_t>> vec =
          grid_creator->cells_within_radius(sne_loc, r_run);
      double mtot = 0.0;
      size_t num_cells = 0;
      for (auto & pair : vec) {
        HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(std::get<0>(pair));
        const double cell_mass =
            (subgrid.hydro_begin() + std::get<1>(pair))
                .get_hydro_variables()
                .get_conserved_mass();
        if (cell_mass > 0.) {
          mtot += cell_mass;
          ++num_cells;
        }
      }

      // The legacy SILCC-like prescription expands the feedback region until
      // it encloses at least 1000 Msun. TIGRESS-like injection always keeps
      // the feedback radius fixed at three cells.
      const double legacy_target_mass = 1000. * solar_mass();
      while (!_use_tigress_like_injection && mtot < legacy_target_mass) {
        r_run += 0.25 * dx;
        vec = grid_creator->cells_within_radius(sne_loc, r_run);
        mtot = 0.;
        num_cells = 0;
        for (auto &pair : vec) {
          HydroDensitySubGrid &radius_subgrid =
              *grid_creator->get_subgrid(std::get<0>(pair));
          const double cell_mass =
              (radius_subgrid.hydro_begin() + std::get<1>(pair))
                  .get_hydro_variables()
                  .get_conserved_mass();
          if (cell_mass > 0.) {
            mtot += cell_mass;
            ++num_cells;
          }
        }
      }

      const double injection_volume = 4. * M_PI * std::pow(r_run, 3.) / 3.;
      const double rho = mtot / injection_volume;
      const double nbar = 1.e-6 * rho / proton_mass();
      double r_st = std::numeric_limits< double >::infinity();
      if (nbar > 0.) {
        if (_use_tigress_like_injection) {
          const double shell_formation_mass =
              get_shell_formation_mass(nbar);
          r_st = std::cbrt(3. * shell_formation_mass /
                           (4. * M_PI * rho));
        } else {
          r_st = parsec() * 19.1 *
                 std::pow(_sne_energy * 1.e-44, 5. / 17.) *
                 std::pow(nbar, -7. / 17.);
        }
      }
      report_injection_mode(r_run, r_st, dx, nbar, num_cells);
      return std::make_tuple(r_run, r_st, nbar, num_cells);

   }



   inline void inject_sne(HydroDensitySubGrid &subgrid, Hydro &hydro, CoordinateVector<double> position,
         double r_inj, double r_st, double nbar, int numcells) {

      const double dx = std::cbrt(subgrid.hydro_begin().get_volume());
      const bool use_momentum =
          use_momentum_injection(r_inj, r_st, dx, nbar);
      const double feedback_mass = get_feedback_mass(r_inj, nbar);
      if (!(feedback_mass > 0.)) {
        return;
      }

      const double energy_51 = _sne_energy / 1.e44;
      const double terminal_momentum =
          (_use_tigress_like_injection
               ? 2.8e5 * std::pow(energy_51, 0.93) *
                     std::pow(nbar, -0.17)
               : 2.6e5 * std::pow(nbar, -2. / 17.) *
                     std::pow(energy_51, 16. / 17.)) *
          solar_mass() * 1.e3;

      for (auto cellit = subgrid.hydro_begin();
           cellit != subgrid.hydro_end(); ++cellit) {

           CoordinateVector<> cellpos = cellit.get_cell_midpoint();

          // is cell within injeciton radius of SNe?
          const double distance = (cellpos - position).norm();
          if (distance <= r_inj) {
                const double cell_mass =
                    cellit.get_hydro_variables().get_conserved_mass();
                if (!(cell_mass > 0.)) {
                    //dont add energy to cell without mass...
                    continue;
                }
                const double mass_fraction = cell_mass / feedback_mass;
                const double legacy_cell_fraction =
                    numcells > 0 ? 1. / numcells : 0.;
                CoordinateVector<> velocity =
                    cellit.get_hydro_variables().get_primitives_velocity();

                if (distance == 0.) {
                  // There is no radial direction at the source. For resolved
                  // events, retain this cell's share as thermal energy.
                  if (!use_momentum) {
                    cellit.get_hydro_variables().set_energy_term(
                        cellit.get_hydro_variables().get_energy_term() +
                        _sne_energy *
                            (_use_tigress_like_injection
                                 ? mass_fraction
                                 : legacy_cell_fraction));
                  }
                  continue;
                }

                const CoordinateVector<> direction =
                    (cellpos - position) / distance;
                if (use_momentum) {
                  // A common radial velocity increment distributes the
                  // terminal scalar momentum in proportion to cell mass.
                  velocity += (terminal_momentum / feedback_mass) * direction;
                  cellit.get_hydro_variables().set_primitives_velocity(velocity);
                  hydro.set_conserved_variables(cellit.get_hydro_variables(),
                                                cellit.get_volume());
                } else if (_use_tigress_like_injection) {
                  // Add exactly the requested kinetic-energy increment even
                  // when the ambient gas already has radial motion. Solving
                  // 0.5*m*dv^2 + m*v_r*dv = dE_kin gives the positive kick.
                  const double radial_velocity =
                      CoordinateVector<>::dot_product(velocity, direction);
                  const double specific_kinetic_energy =
                      kinetic_energy_fraction() * _sne_energy / feedback_mass;
                  const double velocity_increment =
                      -radial_velocity +
                      std::sqrt(radial_velocity * radial_velocity +
                                2. * specific_kinetic_energy);
                  velocity += velocity_increment * direction;
                  cellit.get_hydro_variables().set_primitives_velocity(velocity);
                  hydro.set_conserved_variables(cellit.get_hydro_variables(),
                                                cellit.get_volume());

                  // Accumulate the thermal part so overlapping SNe remain
                  // additive. The simulation applies this term immediately
                  // after all source-distribution calls.
                  cellit.get_hydro_variables().set_energy_term(
                      cellit.get_hydro_variables().get_energy_term() +
                      thermal_energy_fraction() * _sne_energy * mass_fraction);
                } else {
                  // Legacy SILCC-like resolved events are purely thermal.
                  cellit.get_hydro_variables().set_energy_term(
                      cellit.get_hydro_variables().get_energy_term() +
                      _sne_energy * legacy_cell_fraction);
                }
           }

        }




   }

  inline bool uses_tigress_like_injection() const {
    return _use_tigress_like_injection;
  }

  inline void set_tigress_like_injection(const bool value) {
    _use_tigress_like_injection = value;
  }

  /**
   * @brief Destructor.
   */
  ~SupernovaHandler() {}
};

#endif // SUPERNOVAHANDLER_HPP
