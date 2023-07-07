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
 * @file PeakDrivingPhotonSourceDistribution.hpp
 *
 * @brief Peak Driving PhotonSourceDistribution.
 *
 * @author Lewis McCallum (lm261@st-andrews.ac.uk)
 */
#ifndef PEAKDRIVINGPHOTONSOURCEDISTRIBUTION_HPP
#define PEAKDRIVINGPHOTONSOURCEDISTRIBUTION_HPP

#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhotonSourceDistribution.hpp"
#include "RandomGenerator.hpp"
#include "DensitySubGridCreator.hpp"

#include <algorithm>
#include <cinttypes>
#include <fstream>
#include <unistd.h>
#include <vector>

/**
 * @brief Disc patch PhotonSourceDistribution.
 */
class PeakDrivingPhotonSourceDistribution : public PhotonSourceDistribution {
private:
  /*! @brief Lifetime of a source (in s). */
  const double _star_formation_rate;

  /*! @brief Update time interval (in s). */
  const double _update_interval;

  /*! @brief Positions of the sources (in m). */
  std::vector< CoordinateVector<> > _source_positions;

  /*! @brief Remaining lifetime of the sources (in s). */
  std::vector< double > _source_lifetimes;

  std::vector< double > _source_luminosities;

  std::vector < double > _cum_imf;
  std::vector < double > _mass_range;

  /*! @brief Output file for the sources (if applicable). */
  std::ofstream *_output_file;

  std::ofstream *_output_file2;

  /*! @brief Number of updates since the start of the simulation. */
  uint_fast32_t _number_of_updates;

  /*! @brief Indices of the sources (if output is enabled). */
  std::vector< uint_fast32_t > _source_indices;

  /*! @brief Index of the next source to add (if output is enabled). */
  uint_fast32_t _next_index;

  uint_fast32_t _num_sne = 0;



  std::vector< CoordinateVector<> > _to_do_feedback;

  std::vector< double > _r_inj;
  std::vector< double > _r_st;
  std::vector< double > _num_cells_injected;
  std::vector< double > _nbar;


  const double _sne_energy = 1.e44;


  const double _lum_adjust;






  double _excess_mass = 0;

  const double _scaleheight;

  /*! @brief Pseudo-random number generator. */
  RandomGenerator _random_generator;




  double get_r_inj(DensitySubGridCreator< HydroDensitySubGrid > *grid_creator,
                                         CoordinateVector<> sne_loc) {



      HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(sne_loc);

      double cell_vol =  subgrid.get_cell(sne_loc).get_volume();


      double dx = std::pow(cell_vol,1./3.);

      double r_run = 4*dx;


      std::vector<std::pair<uint_fast32_t,uint_fast32_t>> vec;


      vec = grid_creator->cells_within_radius(sne_loc,r_run);


      double mtot = 0.0;
      for (auto & pair : vec) {
        HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(std::get<0>(pair));
        mtot = mtot + (subgrid.hydro_begin() + std::get<1>(pair)).get_hydro_variables().get_conserved_mass();

      }
      if (mtot > 1.988e+33) {
        _num_cells_injected.push_back(268);
        return r_run;
      }

      while (mtot < 1.988e+33) {
        r_run = r_run+(0.25*dx);
        vec = grid_creator->cells_within_radius(sne_loc,r_run);
        mtot = 0.0;
        for (auto & pair : vec) {
          HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(std::get<0>(pair));
          double cell_mass = (subgrid.hydro_begin() + std::get<1>(pair)).get_hydro_variables().get_conserved_mass();
          mtot = mtot + cell_mass;
        }
      }
      _num_cells_injected.push_back(vec.size());
      return r_run;

   }

   double get_r_st(DensitySubGridCreator< HydroDensitySubGrid > *grid_creator,
                                 CoordinateVector<> sne_loc, double r_inj) {

       std::vector<std::pair<uint_fast32_t,uint_fast32_t>> vec;

        vec = grid_creator->cells_within_radius(sne_loc,r_inj);
        double mtot = 0.0;
        for (auto & pair : vec) {
          HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(std::get<0>(pair));
          mtot = mtot + (subgrid.hydro_begin() + std::get<1>(pair)).get_hydro_variables().get_conserved_mass();
        }

        double inj_vol = 1.3333*3.14159265*std::pow(r_inj,3.0);
        double rho = mtot/inj_vol;
        double nbar = 1.e-6*rho/1.67262192e-27;
        _nbar.push_back(nbar);

        double r_st = 3.086e+16 * 19.1 * std::pow(_sne_energy*1.e-44,5./17.) * std::pow(nbar,-7./17);

        return r_st;

      }
    static double kroupa_imf(double mass) {
      if (mass > 0.5) {
        return std::pow(mass,-2.3);
      } else if (mass < 0.08){
        return 2*std::pow(mass,-1.3);
      } else {
        return 25*std::pow(mass,-0.3);
      }
    }

    double integral(double (*f)(double), double a, double b, int n) {
    double step = (b - a) / n;  // width of each small rectangle
    double area = 0.0;  // signed area
    for (int i = 0; i < n; i ++) {
        area += f(a + (i + 0.5) * step) * step; // sum up each small rectangle
    }
    return area;
    }

    double get_single_mass(std::vector<double> mass_range,
           std::vector<double> cum_imf, double rand_num) {

       int Nup = mass_range.size()-1;
       int Nlow=0;
       int mid=(Nup + Nlow)/2;
       while(Nup - Nlow > 1){
         mid=(Nup + Nlow)/2;
         if (rand_num > cum_imf[mid]){
           Nlow = mid;
         } else {
           Nup = mid;
         }
       }

     return (mass_range[Nup] + mass_range[Nlow])/2.0;


    }

    double lum_from_mass(double mass) {
      double lum;
      if (mass < 19.3) {
        lum=0.0;
      } else if (mass < 21.2) {
        lum=std::pow(10,47.865);
      } else if (mass < 23.3) {
        lum=std::pow(10,48.14);
      } else if (mass < 25.4) {
        lum=std::pow(10,48.365);
      } else if (mass < 28.0) {
        lum=std::pow(10,48.54);
      } else if (mass < 30.8) {
        lum=std::pow(10,48.68);
      } else if (mass < 34.1) {
        lum=std::pow(10,48.835);
      } else if (mass < 37.7) {
        lum=std::pow(10,48.99);
      } else if (mass < 41.0) {
        lum=std::pow(10,49.12);
      } else if (mass < 45.2) {
        lum=std::pow(10,49.235);
      } else if (mass < 50.4) {
        lum=std::pow(10,49.34);
      } else if (mass < 56.6) {
        lum=std::pow(10,49.44);
      } else if (mass < 62.3) {
        lum=std::pow(10,49.54);
      } else if (mass < 68.9) {
        lum=std::pow(10,49.635);
      } else if (mass < 87.6) {
        lum=std::pow(10,49.775);
      } else {
        lum=std::pow(10,49.87);
      }


    //  if (mass > 40) {
    //    lum = std::pow(10,49);
    //  } else {
    //    lum = 0.0;
    //  }

      lum = lum*_lum_adjust;
      return lum;

    }



public:
  /**
   * @brief Constructor.
   *
   * @param source_lifetime Lifetime of a source (in s).
   * @param source_luminosity Ionising luminosity of a single source (in s^-1).
   * @param average_number Average number of sources at any given time.
   * @param anchor_x x component of the anchor of the rectangular disk (in m).
   * @param sides_x x side length of the rectangular disk (in m).
   * @param anchor_y  y component of the anchor of the rectangular disk (in m).
   * @param sides_y y side length of the rectangular disk (in m).
   * @param origin_z Origin of the Gaussian disk height distribution (in m).
   * @param scaleheight_z Scale height of the Gaussian disk height distribution
   * (in m).
   * @param seed Seed for the pseudo-random number generator.
   * @param update_interval Time interval in between successive source
   * distribution updates (in s).
   * @param starting_time Start time of the simulation. The distribution is
   * evolved forward in time to this point before it is used (in s).
   * @param output_sources Should the source positions be written to a file?
   */
  inline PeakDrivingPhotonSourceDistribution(
      const double star_formation_rate,
      const int_fast32_t seed, const double update_interval,
      const double starting_time, bool output_sources = false,
      const double sne_energy = 1.e44,
      const double lum_adjust=1.0,
      const double scaleheight=0.0)
      : _star_formation_rate(star_formation_rate), _update_interval(update_interval),
        _output_file(nullptr), _number_of_updates(1), _next_index(0),
        _sne_energy(sne_energy), _lum_adjust(lum_adjust), _scaleheight(scaleheight),
        _random_generator(seed) {




    // form cumulative IMF
    double imf_start = 8.0;
    double imf_end = 120;

    double full_area = integral(kroupa_imf, imf_start, imf_end, 10000);


    uint_fast32_t range_length = 10000;
    for (uint_fast32_t i=0; i< range_length; ++i){
      double step = (imf_end-imf_start)/range_length;
      _mass_range.push_back(imf_start + step*i);
      double part_integral = integral(kroupa_imf, imf_start,imf_start + step*i,10000);
      _cum_imf.push_back(part_integral/full_area);
    }



    if (output_sources) {
      _output_file = new std::ofstream("PeakDriving_source_positions.txt");
      *_output_file << "#time (s)\tx (m)\ty (m)\tz (m)\tevent\tindex\tluminosity\n";
      _output_file->flush();

      _output_file2 = new std::ofstream("TotalLuminosity.txt");
      *_output_file2 << "time (s)\tlum (s^-1)\tnumsne\n";
      _output_file2->flush();

    }

  }



  // -----------------------------------------------

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - source lifetime: Lifetime of a source (default: 20. Myr)
   *  - source luminosity: Ionising luminosity of a single source
   *    (default: 1.e48 s^-1)
   *  - average number of sources: Average number of sources (default: 24)
   *  - anchor x: X position of the anchor of the 2D disc (default: -1. kpc)
   *  - sides x: X side length of the 2D disc (default: 2. kpc)
   *  - anchor y: Y position of the anchor of the 2D disc (default: -1. kpc)
   *  - sides y: Y side length of the 2D disc (default: 2. kpc)
   *  - origin z: Origin of the exponential disc profile in the z direction
   *    (default: 0. pc)
   *  - scaleheight z: Vertical scale height of the exponential disc profile
   *    (default: 63. pc)
   *  - random seed: Random seed used to initialize the random generator that
   *    is used to sample the individual positions (default: 42)
   *  - update interval: Time interval in between successive distribution
   *    updates (default: 0.1 Myr)
   *  - starting time: Starting time of the simulation. The distribution is
   *    evolved forward in time to this point before it is used
   *    (default: 0. Myr)
   *  - output sources: Whether or not to write the source positions to a file
   *    (default: false)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  PeakDrivingPhotonSourceDistribution(ParameterFile &params, Log *log = nullptr)
      : PeakDrivingPhotonSourceDistribution(
            params.get_physical_value< QUANTITY_MASS_RATE >(
                "PhotonSourceDistribution:star formation rate", "0.01 Msol yr^-1"),
            params.get_value< int_fast32_t >(
                "PhotonSourceDistribution:random seed", 42),
            params.get_physical_value< QUANTITY_TIME >(
                "PhotonSourceDistribution:update interval", "0.1 Myr"),
            params.get_physical_value< QUANTITY_TIME >(
                "PhotonSourceDistribution:starting time", "0. Myr"),
            params.get_value< bool >("PhotonSourceDistribution:output sources",
                                     false),
            params.get_physical_value< QUANTITY_ENERGY > (
                "PhotonSourceDistribution:supernova energy", "1.e51 erg"),
            params.get_value< double >("PhotonSourceDistribution:luminosity adjust",1.0),
            params.get_physical_value<QUANTITY_LENGTH> (
              "PhotonSourceDistribution:scale height","0.0 m")) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~PeakDrivingPhotonSourceDistribution() {}

  /**
   * @brief Get the number of sources contained within this distribution.
   *
   * The PhotonSourceDistribution will return exactly this number of valid
   * and unique positions by successive application of operator().
   *
   * @return Number of sources.
   */
  virtual photonsourcenumber_t get_number_of_sources() const {
    return _source_positions.size();
  }


  /**
   * @brief Will the distribution do stellar feedback at the given time?
   *
   * @param current_time Current simulation time (in s).
   * @return True if the star has not exploded yet and its lifetime has been
   * exceeded.
   */
  virtual bool do_stellar_feedback(const double current_time) const {
    return (_to_do_feedback.size() > 0);
  }



  virtual void get_sne_radii(DensitySubGridCreator< HydroDensitySubGrid > &grid_creator) {

       for (uint_fast32_t i = 0; i < _to_do_feedback.size(); ++i) {

         double r_inj = get_r_inj(&grid_creator,_to_do_feedback[i]);

         double r_st = get_r_st(&grid_creator,_to_do_feedback[i],r_inj);


         _r_inj.push_back(r_inj);
         _r_st.push_back(r_st);


       }
  }



  virtual void add_stellar_feedback(HydroDensitySubGrid &subgrid, Hydro &hydro) {



    for (uint_fast32_t i = 0; i < _to_do_feedback.size(); ++i) {

      for (auto cellit = subgrid.hydro_begin();
           cellit != subgrid.hydro_end(); ++cellit) {

           CoordinateVector<> cellpos = cellit.get_cell_midpoint();

           if (cellit.get_hydro_variables().get_primitives_density() == 0) {
             //dont add energy to cell without mass...
             continue;
           }

          // is cell within injeciton radius of SNe?
           if ((cellpos - _to_do_feedback[i]).norm() < _r_inj[i]) {


             double dx = std::pow(cellit.get_volume(),1./3.);
             if (_r_st[i] < 4.*dx) {


              CoordinateVector<> vel_prior =
                       cellit.get_hydro_variables().get_primitives_velocity();

               // Blondin et al
               double mom_to_inj = 2.6e5*std::pow(_nbar[i],-2./17) * std::pow(_sne_energy*1.e-44,16./17.);
               // Msol km/s to kg m/s
               mom_to_inj = mom_to_inj * 2.e30 * 1.e3;

               double m_tot = (_nbar[i]*1e6*1.67e-27)*(4.*3.14159265*std::pow(_r_inj[i],3)/3.);

               double vel_to_inj = mom_to_inj/m_tot;

               CoordinateVector<> direction = (cellpos-_to_do_feedback[i])/((cellpos-_to_do_feedback[i]).norm());


               CoordinateVector<> vel_new = vel_prior + vel_to_inj*direction;

               cellit.get_hydro_variables().set_primitives_velocity(vel_new);


               double density = cellit.get_hydro_variables().get_primitives_density();

               double xH = cellit.get_ionization_variables().get_ionic_fraction(ION_H_n);


               double pressure = 8254.397014*1.e4*density*2./(1.+xH);

              cellit.get_ionization_variables().set_temperature(1.e4);

              cellit.get_hydro_variables().set_primitives_pressure(pressure);


            //  if (vel_new.norm() > 1e6) {
            //    double divisor = vel_new.norm()/1e6;
            //    cellit.get_hydro_variables().set_primitives_velocity(vel_new/divisor);
//
            //  }

              hydro.set_conserved_variables(cellit.get_hydro_variables(), cellit.get_volume());






             }
             else {



               cellit.get_hydro_variables().set_energy_term(_sne_energy/_num_cells_injected[i]);
             }
           }

        }
    }
  }

  virtual void done_stellar_feedback() {

    for (uint_fast32_t i=0; i<_to_do_feedback.size();i++) {

    std::cout << "\n SNe INJECTION HERE: R_inj = " << _r_inj[i] << " R_st = " <<  _r_st[i]
       << " num_cells = " <<  _num_cells_injected[i] << " nbar = "  << _nbar[i] << "\n";
    }



    _to_do_feedback.clear();
    _r_inj.clear();
    _r_st.clear();
    _num_cells_injected.clear();
    _nbar.clear();

  }



  /**
   * @brief Get a valid position from the distribution.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return CoordinateVector of a valid and photon source position (in m).
   */
  virtual CoordinateVector<> get_position(photonsourcenumber_t index) {
    return _source_positions[index];
  }

  /**
   * @brief Get the weight of a photon source.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return Weight of the photon source, used to determine how many photons are
   * emitted from this particular source.
   */
  virtual double get_weight(photonsourcenumber_t index) const {
    return _source_luminosities[index] / get_total_luminosity();
  }

  /**
   * @brief Get the total luminosity of all sources together.
   *
   * @return Total luminosity (in s^-1).
   */


  virtual double get_total_luminosity() const {
    double tot_lum = 0.0;
    for (uint_fast32_t i=0;i<_source_luminosities.size();++i) {
      tot_lum += _source_luminosities[i];
    }
    return tot_lum;
  }


  /**
   * @brief Update the distribution after the system moved to the given time.
   *
   * @param simulation_time Current simulation time (in s).
   * @return True if the distribution changed, false otherwise.
   */
   virtual bool update(DensitySubGridCreator< HydroDensitySubGrid > *grid_creator) override {

    double total_time = _number_of_updates*_update_interval;

    if (_output_file2 != nullptr) {
      double totallum = get_total_luminosity();
      *_output_file2 << total_time << "\t" << totallum << "\t" << _num_sne << "\n";
      _output_file2->flush();

    }




    // clear out sources which no longer exist and add them to SNe todo list
    size_t i = 0;
    while (i < _source_lifetimes.size()) {
      _source_lifetimes[i] -= _update_interval;
      if (_source_lifetimes[i] <= 0.) {
        // remove the element
        if (_output_file != nullptr) {
          *_output_file << total_time << "\t0.\t0.\t0.\t2\t"
                        << _source_indices[i] << "\t0\n";
          _source_indices.erase(_source_indices.begin() + i);
        }

        _to_do_feedback.push_back(_source_positions[i]);
        _source_positions.erase(_source_positions.begin() + i);
        _source_lifetimes.erase(_source_lifetimes.begin() + i);
        _source_luminosities.erase(_source_luminosities.begin() + i);
        _num_sne = _num_sne + 1;



      } else {
        // check the next element
        ++i;
      }
    }



      // form cumulative mass structure

      size_t total_cells = grid_creator->number_of_cells();



      std::vector<double> cumulative_mass(total_cells);

      AtomicValue< size_t > igrid(0);
      i = 0;
      double running_mass = 0.0;
      while (igrid.value() < grid_creator->number_of_original_subgrids()) {
        const size_t this_igrid = igrid.post_increment();
        if (this_igrid < grid_creator->number_of_original_subgrids()) {
          HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(this_igrid);
          for (auto it = subgrid.hydro_begin(); it != subgrid.hydro_end();
               ++it) {

            double cell_mass = it.get_hydro_variables().get_conserved_mass();

            if (_scaleheight == 0.0) {
              running_mass+= cell_mass;
            } else {
              double z_abs = std::abs(it.get_cell_midpoint()[2]);
              running_mass += cell_mass*std::exp(-z_abs/_scaleheight);
            }

            cumulative_mass[i] = running_mass;

            i += 1;


          }
        }
      }






      for (size_t i=0;i<total_cells;i++) {

        cumulative_mass[i] = cumulative_mass[i]/running_mass;
      }




      // 0.073 factor is to take into account we only form stars over 8Msol
      // mass_to_generate in units of Msol to match IMF
      double mass_to_generate = _update_interval*_star_formation_rate/1.988e30*0.073;


       std::cout << "SHOULD BE GENERATING " << mass_to_generate - _excess_mass<< std::endl;
      double mass_generated = 0.0;
      while (mass_generated < mass_to_generate - _excess_mass){
        double m_cur = get_single_mass(_mass_range,_cum_imf,
               _random_generator.get_uniform_random_double());
          std::cout << "MAKING STAR OF MASS " << m_cur <<  std::endl;

          double source_pos_val = _random_generator.get_uniform_random_double();
          CoordinateVector<> cell_midpoint;
          double cell_length = 0;

          AtomicValue< size_t > igrid(0);
          uint_fast32_t i = 0;

          while (igrid.value() < grid_creator->number_of_original_subgrids()) {
            const size_t this_igrid = igrid.post_increment();
            if (this_igrid < grid_creator->number_of_original_subgrids()) {
              HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(this_igrid);
              for (auto it = subgrid.hydro_begin(); it != subgrid.hydro_end();
                   ++it) {

                if (cumulative_mass[i] >= source_pos_val){

                  cell_midpoint = it.get_cell_midpoint();
                  cell_length = std::pow(it.get_volume(),1./3.);
                  goto afterloop;


                }

                i += 1;

              }
            }
          }

         afterloop:

        CoordinateVector<> blur;
        blur[0] = _random_generator.get_uniform_random_double()*cell_length - (0.5*cell_length);
        blur[1] = _random_generator.get_uniform_random_double()*cell_length - (0.5*cell_length);
        blur[2] = _random_generator.get_uniform_random_double()*cell_length - (0.5*cell_length);

        _source_positions.push_back(cell_midpoint + blur);
        double lifetime = 1.e10 * std::pow(m_cur,-2.5) * 3.154e+7;
        double offset =
              _random_generator.get_uniform_random_double() * _update_interval;
        _source_lifetimes.push_back(lifetime-offset);
        _source_luminosities.push_back(lum_from_mass(m_cur));
        if (_output_file != nullptr) {
          _source_indices.push_back(_next_index);
          ++_next_index;
          const CoordinateVector<> &pos = _source_positions.back();
          *_output_file << total_time << "\t" << pos.x() << "\t" << pos.y()
                        << "\t" << pos.z() << "\t1\t"
                        << _source_indices.back() << "\t"
                        << _source_luminosities.back() << "\n";
        }

        mass_generated += m_cur;
      }
      if (mass_generated == 0) {
        _excess_mass = _excess_mass - mass_to_generate;
        std::cout << "Still over mass to generate, decreased excess to " << _excess_mass << std::endl;

      } else {
        _excess_mass = mass_generated - mass_to_generate + _excess_mass;
        std::cout << "OVER SHOT, saving excess of" << _excess_mass << std::endl;

      }









      if (_output_file != nullptr) {
        _output_file->flush();
      }

      ++_number_of_updates;



    return true;
  }


// --------------------------------------

  /**
   * @brief Write the distribution to the given restart file.
   *
   * @param restart_writer RestartWriter to use.
   */
  virtual void write_restart_file(RestartWriter &restart_writer) const {

    restart_writer.write(_star_formation_rate);
    restart_writer.write(_update_interval);
    restart_writer.write(_lum_adjust);
    restart_writer.write(_excess_mass);
    restart_writer.write(_scaleheight);
    _random_generator.write_restart_file(restart_writer);
    {
      const auto size = _source_positions.size();
      restart_writer.write(size);
      for (std::vector< CoordinateVector<> >::size_type i = 0; i < size; ++i) {
        _source_positions[i].write_restart_file(restart_writer);
      }
    }
    {
      const auto size = _source_lifetimes.size();
      restart_writer.write(size);
      for (std::vector< double >::size_type i = 0; i < size; ++i) {
        restart_writer.write(_source_lifetimes[i]);
      }
    }
    {
      const auto size = _source_luminosities.size();
      restart_writer.write(size);
      for (std::vector< double >::size_type i = 0; i < size; ++i) {
        restart_writer.write(_source_luminosities[i]);
      }

    }
    restart_writer.write(_number_of_updates);
    const bool has_output = (_output_file != nullptr);
    restart_writer.write(has_output);
    if (has_output) {
      // store current position in the std::ofstream
      // we want to be able to continue writing from that point
      const auto filepos = _output_file->tellp();
      restart_writer.write(filepos);
      const auto filepos2 = _output_file2->tellp();
      restart_writer.write(filepos2);
      {
        const auto size = _source_indices.size();
        restart_writer.write(size);
        for (std::vector< uint_fast32_t >::size_type i = 0; i < size; ++i) {
          restart_writer.write(_source_indices[i]);
        }
      }
      restart_writer.write(_next_index);
    }
  }

  /**
   * @brief Restart constructor.
   *
   * @param restart_reader Restart file to read from.
   */
  inline PeakDrivingPhotonSourceDistribution(RestartReader &restart_reader)
      : _star_formation_rate(restart_reader.read< double >()),
        _update_interval(restart_reader.read< double >()),
        _lum_adjust(restart_reader.read< double >()),
        _excess_mass(restart_reader.read<double>()),
        _scaleheight(restart_reader.read<double>()),
        _random_generator(restart_reader) {

    {
      const std::vector< CoordinateVector<> >::size_type size =
          restart_reader.read< std::vector< CoordinateVector<> >::size_type >();
      _source_positions.resize(size);
      for (std::vector< CoordinateVector<> >::size_type i = 0; i < size; ++i) {
        _source_positions[i] = CoordinateVector<>(restart_reader);
      }
    }
    {
      const std::vector< double >::size_type size =
          restart_reader.read< std::vector< double >::size_type >();
      _source_lifetimes.resize(size);
      for (std::vector< double >::size_type i = 0; i < size; ++i) {
        _source_lifetimes[i] = restart_reader.read< double >();
      }
    }
    {
      const std::vector< double >::size_type size =
          restart_reader.read< std::vector< double >::size_type >();
      _source_luminosities.resize(size);
      for (std::vector< double >::size_type i = 0; i < size; ++i) {
        _source_luminosities[i] = restart_reader.read< double >();
      }
    }
    _number_of_updates = restart_reader.read< uint_fast32_t >();
    const bool has_output = restart_reader.read< bool >();
    if (has_output) {
      const std::streampos filepos = restart_reader.read< std::streampos >();
      // truncate the original file to the size we were at
      if (truncate("PeakDriving_source_positions.txt", filepos) != 0) {
        cmac_error("Error while truncating output file!");
      }
      // now open the file in append mode
      _output_file = new std::ofstream("PeakDriving_source_positions.txt",
                                       std::ios_base::app);

      const std::streampos filepos2 = restart_reader.read< std::streampos >();
                                       // truncate the original file to the size we were at
      if (truncate("TotalLuminosity.txt", filepos2) != 0) {
              cmac_error("Error while truncating output file!");
            }
                                       // now open the file in append mode
      _output_file2 = new std::ofstream("TotalLuminosity.txt",
                                            std::ios_base::app);


      {
        const std::vector< uint_fast32_t >::size_type size =
            restart_reader.read< std::vector< uint_fast32_t >::size_type >();
        _source_indices.resize(size);
        for (std::vector< uint_fast32_t >::size_type i = 0; i < size; ++i) {
          _source_indices[i] = restart_reader.read< uint_fast32_t >();
        }
      }
      _next_index = restart_reader.read< uint_fast32_t >();
    }


        // form cumulative IMF
        double imf_start = 8.0;
        double imf_end = 120;

        double full_area = integral(kroupa_imf, imf_start, imf_end, 10000);


        uint_fast32_t range_length = 10000;
        for (uint_fast32_t i=0; i< range_length; ++i){
          double step = (imf_end-imf_start)/range_length;
          _mass_range.push_back(imf_start + step*i);
          double part_integral = integral(kroupa_imf, imf_start,imf_start + step*i,10000);
          _cum_imf.push_back(part_integral/full_area);
        }
  }
};

#endif // PEAKDRIVINGPHOTONSOURCEDISTRIBUTION_HPP
