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
 * @file MixedDrivingPhotonSourceDistribution.hpp
 *
 * @brief Mixed Driving PhotonSourceDistribution.
 *
 * @author Lewis McCallum (lm261@st-andrews.ac.uk)
 */
#ifndef MIXEDDRIVINGPHOTONSOURCEDISTRIBUTION_HPP
#define MIXEDDRIVINGPHOTONSOURCEDISTRIBUTION_HPP

#include "Log.hpp"
#include "ExternalPotential.hpp"
#include "GalacticShearingBox.hpp"
#ifdef HAVE_HDF5
#include "HDF5Tools.hpp"
#endif
#include "ParameterFile.hpp"
#include "PhotonSourceDistribution.hpp"
#include "RandomGenerator.hpp"
#include "DensitySubGridCreator.hpp"
#include "SupernovaHandler.hpp"
#include "WMBasicPhotonSourceSpectrum.hpp"
#include "PowerLawPhotonSourceSpectrum.hpp"
#include "Pegase3PhotonSourceSpectrum.hpp"
#include "PhotonSourceSpectrum.hpp"

#include <algorithm>
#include <cinttypes>
#include <fstream>
#include <unistd.h>
#include <vector>

/**
 * @brief Wrap a moving source through ordinary periodic box boundaries.
 *
 * @param position Source position to update.
 * @param box Simulation box.
 * @param periodicity Periodicity flags for each axis.
 */
inline void wrap_mixed_driving_source_position(
    CoordinateVector<> &position, const Box<> &box,
    const CoordinateVector< bool > &periodicity) {
  const CoordinateVector<> &anchor = box.get_anchor();
  const CoordinateVector<> &sides = box.get_sides();
  for (uint_fast8_t axis = 0; axis < 3; ++axis) {
    if (periodicity[axis]) {
      double offset = std::fmod(position[axis] - anchor[axis], sides[axis]);
      if (offset < 0.) {
        offset += sides[axis];
      }
      position[axis] = anchor[axis] + offset;
    }
  }
}



/**
 * @brief Disc patch PhotonSourceDistribution.
 */
class MixedDrivingPhotonSourceDistribution : public PhotonSourceDistribution {
private:
  /*! @brief Lifetime of a source (in s). */
  const double _star_formation_rate;

  /*! @brief Update time interval (in s). */
  const double _update_interval;

  /*! @brief Positions of the sources (in m). */
  std::vector< CoordinateVector<> > _source_positions;

  /*! @brief Velocities inherited from the gas at source birth (in m s^-1). */
  std::vector< CoordinateVector<> > _source_velocities;

  /*! @brief Remaining lifetime of the sources (in s). */
  std::vector< double > _source_lifetimes;

  std::vector< double > _source_luminosities;

  /*! @brief Initial stellar masses (in Msol). */
  std::vector< double > _source_masses;

  /*! @brief Supernovae that occurred since the preceding HDF5 snapshot. */
  std::vector< CoordinateVector<> > _snapshot_supernova_positions;
  std::vector< double > _snapshot_supernova_times;

  std::vector<int> _to_delete;



  std::vector < double > _cum_imf;
  std::vector < double > _mass_range;

  /*! @brief Output file for the sources (if applicable). */
  std::ofstream *_output_file;

  std::ofstream *_output_file2;

  /*! @brief Number of updates since the start of the simulation. */
  uint_fast32_t _number_of_updates;

  /*! @brief Indices of the sources (if output is enabled). */
  std::vector< uint_fast32_t > _source_indices;

  std::vector<int> _spectrum_index;
  std::vector<PhotonSourceSpectrum*> _all_spectra;

  std::vector<double> stellarMasses = {57.95, 46.94, 38.08, 34.39, 30.98, 28.0, 25.29, 22.90, 20.76, 18.80, 17.08, 15.55};
  std::vector<double> temperatures = {44852, 42857, 40862, 39865, 38867, 37870, 36872, 35874, 34877, 33879, 32882, 31884};

  std::vector<double> avail_temps = {32000, 34000, 34000, 35000, 36000, 37000, 39000, 39000, 40000, 41000,42000,43000,44000,45000};

  /*! @brief Index of the next source to add (if output is enabled). */
  uint_fast32_t _next_index;


  std::vector< CoordinateVector<> > _to_do_feedback;

  std::vector< double > _r_inj;
  std::vector< double > _r_st;
  std::vector< double > _num_cells_injected;
  std::vector< double > _nbar;

  const double _sne_energy = 1.e44;

  const double _lum_adjust;

  double _excess_mass = 0;

  const double _scaleheight;

  const double _peak_fraction;

  /*! @brief Density exponent used for the density-selected birth CDF. */
  const double _clustering_factor;

  /*! @brief Whether sources move after formation. */
  const bool _float_sources;

  /*! @brief Pseudo-random number generator. */

  double init_running_mass;
  uint_fast32_t _num_sne = 0;

  double _holmes_time;
  double _holmes_sh;
  double _holmes_lum;
  uint_fast32_t _number_of_holmes;


  bool _read_file;
  std::string _filename;
  double _time;

  int type1done = 0;


  double _total_time = 0.;

  bool _holmes_added = false;

  double _last_sf = 0.;

  /*! @brief Whether moving sources changed the source-copy hierarchy. */
  bool _sources_changed = false;

  /*! @brief Whether file-loaded sources still need their local gas velocity. */
  bool _initialize_imported_source_velocities = false;

  RandomGenerator _random_generator;

  SupernovaHandler *novahandler;

  Log *_log;


    static double kroupa_imf(double mass) {
      if (mass > 0.5) {
        return std::pow(mass,-2.3);
      } else if (mass > 0.08){
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

// Function to perform linear interpolation for descending xVals
double interpolate(double x, const std::vector<double>& xVals, const std::vector<double>& yVals) {
    // Ensure inputs are valid
    if (xVals.size() != yVals.size() || xVals.empty()) {
        throw std::invalid_argument("Invalid input: xVals and yVals must have the same size and cannot be empty.");
    }


    if (x > xVals.front()){
      return yVals.front();
    } else if (x < xVals.back()) {
      return yVals.back();
    }

    // Find the interval containing x
    for (size_t i = 0; i < xVals.size() - 1; ++i) {
        if (xVals[i] >= x && x >= xVals[i + 1]) {
            // Perform linear interpolation
            double x1 = xVals[i];
            double x2 = xVals[i + 1];
            double y1 = yVals[i];
            double y2 = yVals[i + 1];
            return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
        }
    }

    // If we reach here, x was not in any valid interval
    throw std::logic_error("Could not interpolate: x is not within any interval.");
}


size_t findClosestIndex(double value, const std::vector<double>& values) {
    if (values.empty()) {
        throw std::invalid_argument("The values vector cannot be empty.");
    }

    size_t closestIndex = 0;
    double minDifference = std::abs(value - values[0]);

    for (size_t i = 1; i < values.size(); ++i) {
        double difference = std::abs(value - values[i]);
        if (difference < minDifference) {
            minDifference = difference;
            closestIndex = i;
        }
    }

    return closestIndex;
}

    double lum_from_mass(double mass) {

      std::vector<double> tablemasses= {57.95, 46.94, 38.08, 34.39, 30.98, 28.0, 25.29, 22.90, 20.76, 18.80, 17.08, 15.55};
      std::vector<double> tablelums = {49.64,49.44,49.22,49.10,48.99,48.88,48.75,48.61,48.44,48.27,48.06,47.88};
      double lum = 0.0;

      if (mass > tablemasses.front()){
        lum = tablelums.front();
        lum = std::pow(10,lum);
        lum = lum*_lum_adjust;
        return lum;
      } else if (mass < tablemasses.back()) {
        return 0.0;
      }


      // Find the interval containing x
      for (size_t i = 0; i < tablemasses.size() - 1; ++i) {
          if (tablemasses[i] >= mass && mass >= tablemasses[i + 1]) {
              // Perform linear interpolation
              double x1 = tablemasses[i];
              double x2 = tablemasses[i + 1];
              double y1 = tablelums[i];
              double y2 = tablelums[i + 1];
              lum = y1 + (mass - x1) * (y2 - y1) / (x2 - x1);
          }
      }

      lum = std::pow(10,lum);
      lum = lum*_lum_adjust;
      return lum;

    }

    /** Recover a best available mass for sources in legacy restart files. */
    double mass_from_luminosity(const double luminosity) const {
      if (luminosity <= 0. || _lum_adjust <= 0. ||
          (_holmes_lum > 0. &&
           std::abs(luminosity - _holmes_lum) < 1.e-12 * _holmes_lum)) {
        return 0.;
      }
      const std::vector<double> masses =
          {57.95, 46.94, 38.08, 34.39, 30.98, 28.0, 25.29, 22.90,
           20.76, 18.80, 17.08, 15.55};
      const std::vector<double> log_luminosities =
          {49.64, 49.44, 49.22, 49.10, 48.99, 48.88, 48.75, 48.61,
           48.44, 48.27, 48.06, 47.88};
      const double log_luminosity =
          std::log10(std::max(luminosity / _lum_adjust, 1.e-99));
      if (log_luminosity >= log_luminosities.front()) {
        return masses.front();
      }
      for (size_t i = 0; i + 1 < masses.size(); ++i) {
        if (log_luminosities[i] >= log_luminosity &&
            log_luminosity >= log_luminosities[i + 1]) {
          const double fraction =
              (log_luminosity - log_luminosities[i]) /
              (log_luminosities[i + 1] - log_luminosities[i]);
          return masses[i] + fraction * (masses[i + 1] - masses[i]);
        }
      }
      return 0.;
    }

    void initialize_spectra(Log *log) {
      _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(32000,25,log));
      _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(34000,25,log));
      _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(34000,25,log));
      _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(35000,40,log));
      _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(36000,25,log));
      _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(37000,25,log));
      _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(39000,25,log));
      _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(39000,25,log));
      _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(40000,25,log));
      _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(41000,40,log));
      _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(42000,40,log));
      _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(43000,40,log));
      _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(44000,40,log));
      _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(45000,40,log));
      _all_spectra.push_back(new Pegase3PhotonSourceSpectrum(1e10,0.02,log));
    }

    /** Recover the spectrum bin from a luminosity in legacy restart files. */
    size_t spectrum_index_from_luminosity(const double luminosity) {
      if (_holmes_lum > 0. &&
          std::abs(luminosity - _holmes_lum) < 1.e-12 * _holmes_lum) {
        return 14;
      }
      const double mass = mass_from_luminosity(luminosity);
      return findClosestIndex(interpolate(mass, stellarMasses, temperatures),
                              avail_temps);
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
  inline MixedDrivingPhotonSourceDistribution(
      const double star_formation_rate,
      const int_fast32_t seed, const double update_interval,
      const double starting_time, bool output_sources = false,
      const double sne_energy = 1.e44,
      const double lum_adjust=1.0,
      const double scaleheight=0.0,
      const double peak_fraction=0.5,
      const double clustering_factor=2.0,
      const bool float_sources=false,
      const double holmes_time=0.0,
      const double holmes_sh=3e18,
      const double holmes_lum=5e46,
      const uint_fast32_t number_of_holmes=200,
      const bool read_file=false,
      const std::string filename="sources.txt",
      const double time=100.,
      Log *log=nullptr)
      : _star_formation_rate(star_formation_rate), _update_interval(update_interval),
        _output_file(nullptr), _number_of_updates(1), _next_index(0),
        _sne_energy(sne_energy), _lum_adjust(lum_adjust), _scaleheight(scaleheight),
        _peak_fraction(peak_fraction), _clustering_factor(clustering_factor),
        _float_sources(float_sources), _holmes_time(holmes_time),
        _holmes_sh(holmes_sh),_holmes_lum(holmes_lum),_number_of_holmes(number_of_holmes),
        _read_file(read_file), _filename(filename), _time(time),
        _random_generator(seed), _log(log){

    if (_clustering_factor <= 0.) {
      cmac_error("The star-formation clustering factor must be positive.");
    }

    novahandler = new SupernovaHandler(_sne_energy);

    

    initialize_spectra(log);






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
      _output_file = new std::ofstream("MixedDriving_source_positions.txt");
      *_output_file << "#time (s)\tx (m)\ty (m)\tz (m)\tevent\tindex\tluminosity\tMass\ttype\n";
      _output_file->flush();

      _output_file2 = new std::ofstream("TotalLuminosity.txt");
      *_output_file2 << "time (s)\tlum (s^-1)\tnumsne\tSFR\n";
      _output_file2->flush();

    }

    if (_read_file){
      std::ifstream file;
      file.open(_filename);
      if (!file.is_open()) {
        cmac_error("Could not open file \"%s\"!", _filename.c_str());
      } else {
        std::cout << "Opened file - " << _filename << std::endl;
      }


      double time_val,posx,posy,posz,luminosity,mass;
      int event,index;


      std::string dummyLine,star_type;

      std::getline(file, dummyLine);

      time_val = 0.0;

      while (!file.eof() && time_val <= time) {
        file >> time_val >> posx >> posy >> posz >> event >> index >> luminosity >> mass >> star_type;

        if (event == 2) {
          _to_delete.push_back(index);
        }
      }



    file.close();

    file.open(_filename);


    std::getline(file, dummyLine);

    double a0z = 9.955209529401348;
    double a1z = -3.3370109454102326;
    double a2z = 0.8116654874025604;
    time_val = 0.0;
    while (!file.eof() && time_val <= time) {
      file >> time_val >> posx >> posy >> posz >> event >> index >> luminosity >> mass >> star_type;
      if (event == 1) {
        if (std::find(_to_delete.begin(), _to_delete.end(), index) == _to_delete.end()) {
            
            _source_positions.push_back(CoordinateVector<double>(posx,posy,posz));
            // Historical source files do not contain velocities. These are
            // initialized from the local gas on the first grid update.
            _source_velocities.push_back(CoordinateVector<double>());

            _source_luminosities.push_back(luminosity);
            _source_masses.push_back(mass);
            
            double lifetime = a0z + a1z*std::log10(mass) + a2z*(std::log10(mass)*std::log10(mass));
            lifetime = std::pow(10.0,lifetime);
            lifetime = lifetime*3.154e+7;
            lifetime -= (_time - time_val);
            _source_lifetimes.push_back(lifetime);
            _source_indices.push_back(_next_index);
            if (star_type == "HOLMES") {
              _spectrum_index.push_back(14);
              if (_output_file != nullptr) {
                *_output_file << _total_time << "\t" << posx << "\t" << posy
                          << "\t" << posz << "\t1\t"
                          << _source_indices.back() << "\t"
                          << _source_luminosities.back() << "\t"
                          << mass << "\t"
                          << "HOLMES\n";
              }
            } else {
              double interpolatedTemp = interpolate(mass, stellarMasses, temperatures);
              size_t closestIndex = findClosestIndex(interpolatedTemp, avail_temps);
              std::cout << "Adding star of mass " << mass << " temp of " << interpolatedTemp << " for spec index " << closestIndex << std::endl;
              _spectrum_index.push_back(closestIndex);
              if (_output_file != nullptr) {
                *_output_file << _total_time << "\t" << posx << "\t" << posy
                          << "\t" << posz << "\t1\t"
                          << _source_indices.back() << "\t"
                          << _source_luminosities.back() << "\t"
                          << mass << "\t"
                          << "OSTAR\n";
              }
            }
            ++_next_index;

      }
    }
  }



  file.close();
  _initialize_imported_source_velocities = true;
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
   *  - clustering factor: Density exponent for the CDF used to select
   *    density-born stars (default: 2; one restores mass weighting)
   *  - float sources: Move stars with the gas gravitational acceleration and
   *    Galactic-frame source terms (default: false)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  MixedDrivingPhotonSourceDistribution(ParameterFile &params, Log *log = nullptr)
      : MixedDrivingPhotonSourceDistribution(
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
              "PhotonSourceDistribution:scale height","0.0 m"),
            params.get_value< double >("PhotonSourceDistribution:peak fraction",0.5),
            params.get_value< double >("PhotonSourceDistribution:clustering factor",2.0),
            params.get_value< bool >("PhotonSourceDistribution:float sources",false),
            params.get_physical_value<QUANTITY_TIME> (
                "PhotonSourceDistribution:holmes time","50 Myr"),
            params.get_physical_value<QUANTITY_LENGTH>(
                "PhotonSourceDistribution:holmes height","700 pc"),
            params.get_physical_value<QUANTITY_FREQUENCY>(
                "PhotonSourceDistribution:holmes luminosity","5e46 s^-1"),
            params.get_value<double>("PhotonSourceDistribution:number of holmes",200),
            params.get_value<bool>("PhotonSourceDistribution:read file",false),
            params.get_value<std::string>("PhotonSourceDistribution:filename","SourceFile.txt"),
            params.get_physical_value<QUANTITY_TIME>("PhotonSourceDistribution:time","0.0 Myr"),
            log) {
    novahandler->set_tigress_like_injection(params.get_value< bool >(
        "SupernovaHandler:TIGRESS like injection", true));
  }

  /**
   * @brief Virtual destructor.
   */
  virtual ~MixedDrivingPhotonSourceDistribution() {}

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



  virtual void set_tigress_like_supernova_injection(const bool value) {
    novahandler->set_tigress_like_injection(value);
  }

  virtual void get_sne_radii(DensitySubGridCreator< HydroDensitySubGrid > &grid_creator) {

       for (uint_fast32_t i = 0; i < _to_do_feedback.size(); ++i) {

        double r_inj,r_st,nbar,num_inj;


        std::tie(r_inj,r_st,nbar,num_inj) = novahandler->get_r_inj(&grid_creator,_to_do_feedback[i]);

         _r_inj.push_back(r_inj);
         _r_st.push_back(r_st);
         _nbar.push_back(nbar);
         _num_cells_injected.push_back(num_inj);


       }
  }



  virtual void add_stellar_feedback(HydroDensitySubGrid &subgrid, Hydro &hydro) {



    for (uint_fast32_t i = 0; i < _to_do_feedback.size(); ++i) {




      novahandler->inject_sne(subgrid, hydro, _to_do_feedback[i], _r_inj[i],_r_st[i],_nbar[i],_num_cells_injected[i]);


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
   * @brief Append live stars and supernovae since the previous snapshot.
   */
  virtual void write_snapshot_metadata(const std::string &filename,
                                       const double simulation_time) override {
#ifdef HAVE_HDF5
    if (filename.empty()) {
      return;
    }
    cmac_assert(_source_positions.size() == _source_luminosities.size());
    cmac_assert(_source_positions.size() == _source_masses.size());
    HDF5Tools::HDF5File file =
        HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_APPEND);

    HDF5Tools::HDF5Group sources =
        HDF5Tools::create_group(file, "PhotonSources");
    uint32_t number_of_sources = _source_positions.size();
    std::string coordinate_units = "m";
    std::string luminosity_units = "s^-1";
    std::string mass_units = "Msol";
    HDF5Tools::write_attribute< uint32_t >(
        sources, "NumberOfSources", number_of_sources);
    HDF5Tools::write_attribute< std::string >(
        sources, "CoordinateUnits", coordinate_units);
    HDF5Tools::write_attribute< std::string >(
        sources, "IonizingLuminosityUnits", luminosity_units);
    HDF5Tools::write_attribute< std::string >(
        sources, "MassUnits", mass_units);
    if (!_source_positions.empty()) {
      HDF5Tools::write_dataset< CoordinateVector<> >(
          sources, "Coordinates", _source_positions);
      HDF5Tools::write_dataset< double >(
          sources, "IonizingLuminosity", _source_luminosities);
      HDF5Tools::write_dataset< double >(sources, "Mass", _source_masses);
    }
    HDF5Tools::close_group(sources);

    HDF5Tools::HDF5Group supernovae =
        HDF5Tools::create_group(file, "SupernovaEvents");
    uint32_t number_of_supernovae =
        _snapshot_supernova_positions.size();
    std::string time_units = "s";
    double snapshot_time = simulation_time;
    HDF5Tools::write_attribute< uint32_t >(
        supernovae, "NumberOfEvents", number_of_supernovae);
    HDF5Tools::write_attribute< std::string >(
        supernovae, "CoordinateUnits", coordinate_units);
    HDF5Tools::write_attribute< std::string >(
        supernovae, "TimeUnits", time_units);
    HDF5Tools::write_attribute< double >(
        supernovae, "SnapshotTime", snapshot_time);
    if (!_snapshot_supernova_positions.empty()) {
      HDF5Tools::write_dataset< CoordinateVector<> >(
          supernovae, "Coordinates", _snapshot_supernova_positions);
      HDF5Tools::write_dataset< double >(
          supernovae, "Time", _snapshot_supernova_times);
    }
    HDF5Tools::close_group(supernovae);
    HDF5Tools::close_file(file);

    // Only clear after the HDF5 file was closed successfully.
    _snapshot_supernova_positions.clear();
    _snapshot_supernova_times.clear();
#else
    (void)filename;
    (void)simulation_time;
#endif
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


  double get_photon_frequency(RandomGenerator &random_generator,
  photonsourcenumber_t index) {

    return _all_spectra[_spectrum_index[index]]->get_random_frequency(random_generator,0.0);

}

  double get_photon_frequency_weighted(RandomGenerator &random,
      photonsourcenumber_t index, double uniform_fraction, double &weight) {
    return _all_spectra[_spectrum_index[index]]->get_random_frequency_weighted(
        random, uniform_fraction, weight);
  }


  /**
   * @brief Update the distribution after the system moved to the given time.
   *
   * @param simulation_time Current simulation time (in s).
   * @return True if the distribution changed, false otherwise.
   */
   virtual bool update(DensitySubGridCreator< HydroDensitySubGrid > *grid_creator, double actual_timestep) override {

    _total_time += actual_timestep;

    // Historical source catalogues contain positions but no velocities.  A
    // zero velocity makes every imported star fall towards the centre as soon
    // as source floating is enabled.  Give these stars the same local gas
    // velocity that newly formed stars inherit at birth.
    if (_initialize_imported_source_velocities) {
      for (size_t isource = 0; isource < _source_positions.size(); ++isource) {
        if (grid_creator->get_box().inside(_source_positions[isource])) {
          HydroDensitySubGrid &subgrid =
              *grid_creator->get_subgrid(_source_positions[isource]);
          const auto cell = subgrid.get_hydro_cell(_source_positions[isource]);
          _source_velocities[isource] =
              cell.get_hydro_variables().get_primitives_velocity();
        }
      }
      _initialize_imported_source_velocities = false;
    }


    bool updated = _sources_changed;
    _sources_changed = false;





    // clear out sources which no longer exist and add them to SNe todo list
    size_t i = 0;
    while (i < _source_lifetimes.size()) {
      _source_lifetimes[i] -= actual_timestep;
      if (_source_lifetimes[i] <= 0.) {
        // remove the element
        if (_output_file != nullptr) {
          *_output_file << _total_time << "\t" << _source_positions[i].x()
                        << "\t" << _source_positions[i].y() << "\t"
                        << _source_positions[i].z() << "\t2\t"
                        << _source_indices[i] << "\t0\t0\tSNe\n";
        }
        _to_do_feedback.push_back(_source_positions[i]);
        _snapshot_supernova_positions.push_back(_source_positions[i]);
        _snapshot_supernova_times.push_back(_total_time);
        _source_positions.erase(_source_positions.begin() + i);
        _source_velocities.erase(_source_velocities.begin() + i);
        _source_lifetimes.erase(_source_lifetimes.begin() + i);
        _source_luminosities.erase(_source_luminosities.begin() + i);
        _source_masses.erase(_source_masses.begin() + i);
        _spectrum_index.erase(_spectrum_index.begin() + i);
        _source_indices.erase(_source_indices.begin() + i);
        _num_sne = _num_sne + 1;
        updated = true;



      } else {
        // check the next element
        ++i;
      }
    }

    //get simulation box limits

    double _anchor_x = grid_creator->get_box().get_anchor()[0];
    double _anchor_y = grid_creator->get_box().get_anchor()[1];

    double _sides_x = grid_creator->get_box().get_sides()[0];
    double _sides_y = grid_creator->get_box().get_sides()[1];


    double area_kpc = _sides_x*_sides_y/(3.086e+19)/(3.086e+19);


    int should_have_done = int(4*area_kpc*_total_time/3.15576e13)*0;

    int do_type1 = should_have_done-type1done;

    if (do_type1 > 0) {

    for (int i=0;i<do_type1;i++) {
      double x =
       _anchor_x + _random_generator.get_uniform_random_double() * _sides_x;
      double y =
       _anchor_y + _random_generator.get_uniform_random_double() * _sides_y;

      double z =
       325 * 3.086e16 *
           std::sqrt(-2. *
                     std::log(_random_generator.get_uniform_random_double())) *
           std::cos(2. * M_PI *
                    _random_generator.get_uniform_random_double());
            if (grid_creator->get_box().inside(CoordinateVector<double>(x,y,z))) {
              const CoordinateVector<> position(x, y, z);
              _to_do_feedback.push_back(position);
              _snapshot_supernova_positions.push_back(position);
              _snapshot_supernova_times.push_back(_total_time);
      }

      //dotype1
      type1done += 1;
    }

  }



    if ((_total_time > _holmes_time) && (!_holmes_added)) {

      for (uint_fast32_t i=0; i<_number_of_holmes; ++i) {

        double x =
         _anchor_x + _random_generator.get_uniform_random_double() * _sides_x;
        double y =
         _anchor_y + _random_generator.get_uniform_random_double() * _sides_y;
     // we use the Box-Muller method to sample the Gaussian
        double z =
         _holmes_sh *
             std::sqrt(-2. *
                       std::log(_random_generator.get_uniform_random_double())) *
             std::cos(2. * M_PI *
                      _random_generator.get_uniform_random_double());
        if (std::abs(z) >= grid_creator->get_box().get_sides()[2]/2.) {
          continue;
        }

        _source_positions.push_back(CoordinateVector<double>(x,y,z));
        auto grid = grid_creator->get_subgrid(_source_positions.back());
        auto cell = (*grid).get_hydro_cell(_source_positions.back());
        _source_velocities.push_back(
            cell.get_hydro_variables().get_primitives_velocity());

        double lifetime = 1e99;

        _source_lifetimes.push_back(lifetime);
        _source_luminosities.push_back(_holmes_lum);
        // HOLMES represents a population rather than one massive star.
        _source_masses.push_back(0.);
        _source_indices.push_back(_next_index);
        ++_next_index;
        if (_output_file != nullptr) {
          const CoordinateVector<> &pos = _source_positions.back();
          *_output_file << _total_time << "\t" << pos.x() << "\t" << pos.y()
                        << "\t" << pos.z() << "\t1\t"
                        << _source_indices.back() << "\t"
                        << _source_luminosities.back() << "\t"
                        << "HOLMES\n";
        }

      }
      _holmes_added = true;
      updated = true;
    }


    if (_total_time - _last_sf >= _update_interval) {

      // Use the actual elapsed time so that a late hydro step does not discard
      // star formation between the scheduled update and the current time.
      const double star_formation_interval = _total_time - _last_sf;




      // form cumulative mass structure

      size_t total_cells = grid_creator->number_of_cells();



      std::vector<double> cumulative_weight(total_cells);

      AtomicValue< size_t > igrid(0);
      i = 0;
      double running_mass = 0.0;
      double running_weight = 0.0;
      while (igrid.value() < grid_creator->number_of_original_subgrids()) {
        const size_t this_igrid = igrid.post_increment();
        if (this_igrid < grid_creator->number_of_original_subgrids()) {
          HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(this_igrid);
          for (auto it = subgrid.hydro_begin(); it != subgrid.hydro_end();
               ++it) {

            double cell_mass = it.get_hydro_variables().get_conserved_mass();
            double cell_z = std::abs(it.get_cell_midpoint()[2]);

            if (cell_z > 6.171e18) {
              cell_mass = 0.0;
            }

            running_mass+= cell_mass;
            // A density exponent of one reproduces the original mass-weighted
            // CDF.  The volume factor makes larger exponents resolution
            // independent: weight is proportional to rho^p dV.
            double cell_weight = 0.;
            if (cell_mass > 0.) {
              cell_weight = it.get_volume() * std::pow(
                  std::max(0., it.get_hydro_variables().get_primitives_density()),
                  _clustering_factor);
            }
            running_weight += cell_weight;
            cumulative_weight[i] = running_weight;

            i += 1;


          }
        }
      }
      if (_number_of_updates == 1) {
          init_running_mass = running_mass;
      }

      if (running_mass <= 0. || init_running_mass <= 0. ||
          running_weight <= 0.) {
        if (_log != nullptr) {
          _log->write_warning(
              "No positive mass or density weight available for star formation.");
        }
        _last_sf = _total_time;
        ++_number_of_updates;
        return updated;
      }


      std::cout << "SFR Ratio from init =  " << running_mass/init_running_mass << std::endl;








      for (size_t i=0;i<total_cells;i++) {

        cumulative_weight[i] = cumulative_weight[i] / running_weight;
      }


      if (_output_file2 != nullptr) {
        double totallum = get_total_luminosity();
        *_output_file2 << _total_time << "\t" << totallum << "\t" << _num_sne << "\t" << _star_formation_rate*std::pow(running_mass/init_running_mass,1.4) << "\n";
        _output_file2->flush();

      }






      // 0.207 is the Kroupa IMF mass fraction in stars above 8 Msol.
      // mass_to_generate in units of Msol to match IMF
      const double mass_to_generate =
          star_formation_interval * _star_formation_rate / 1.988e30 * 0.207 *
          std::pow(running_mass / init_running_mass, 1.4);


       std::cout << "SHOULD BE GENERATING " << mass_to_generate - _excess_mass<< std::endl;
      double mass_generated = 0.0;
      while (mass_generated < mass_to_generate - _excess_mass){
        double m_cur = get_single_mass(_mass_range,_cum_imf,
               _random_generator.get_uniform_random_double());
          double use_density = _random_generator.get_uniform_random_double();
          if(use_density < _peak_fraction) {






          double source_pos_val = _random_generator.get_uniform_random_double();
          CoordinateVector<> cell_midpoint;
          CoordinateVector<> source_velocity;
          double cell_length = 0;

          AtomicValue< size_t > igrid(0);
          uint_fast32_t i = 0;

          while (igrid.value() < grid_creator->number_of_original_subgrids()) {
            const size_t this_igrid = igrid.post_increment();
            if (this_igrid < grid_creator->number_of_original_subgrids()) {
              HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(this_igrid);
              for (auto it = subgrid.hydro_begin(); it != subgrid.hydro_end();
                   ++it) {

                if (cumulative_weight[i] >= source_pos_val){

                  cell_midpoint = it.get_cell_midpoint();
                  source_velocity =
                      it.get_hydro_variables().get_primitives_velocity();
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
        _source_velocities.push_back(source_velocity);

      } else {
        double x =
         _anchor_x + _random_generator.get_uniform_random_double() * _sides_x;
        double y =
         _anchor_y + _random_generator.get_uniform_random_double() * _sides_y;
     // we use the Box-Muller method to sample the Gaussian
        double z =
         _scaleheight *
             std::sqrt(-2. *
                       std::log(_random_generator.get_uniform_random_double())) *
             std::cos(2. * M_PI *
                      _random_generator.get_uniform_random_double());

        _source_positions.push_back(CoordinateVector<double>(x,y,z));
        CoordinateVector<> source_velocity;
        if (grid_creator->get_box().inside(_source_positions.back())) {
          auto grid = grid_creator->get_subgrid(_source_positions.back());
          auto cell = (*grid).get_hydro_cell(_source_positions.back());
          source_velocity =
              cell.get_hydro_variables().get_primitives_velocity();
        }
        _source_velocities.push_back(source_velocity);

      }
        double a0z = 9.955209529401348;
        double a1z = -3.3370109454102326;
        double a2z = 0.8116654874025604;
       // double lifetime = 1.e10 * std::pow(m_cur,-2.5) * 3.154e+7;
       double lifetime = a0z + a1z*std::log10(m_cur) + a2z*(std::log10(m_cur)*std::log10(m_cur));
       lifetime = std::pow(10.0,lifetime);
       lifetime = lifetime*3.154e+7;
        double offset = _random_generator.get_uniform_random_double() *
                        star_formation_interval;
        _source_lifetimes.push_back(lifetime-offset);
        _source_luminosities.push_back(lum_from_mass(m_cur));
        _source_masses.push_back(m_cur);
        _source_indices.push_back(_next_index);
        ++_next_index;
        double interpolatedTemp = interpolate(m_cur, stellarMasses, temperatures);
        size_t closestIndex = findClosestIndex(interpolatedTemp, avail_temps);
        _spectrum_index.push_back(closestIndex);
        std::cout << "MAKING STAR OF MASS " << m_cur << " temp = " << interpolatedTemp << " specindex = " << closestIndex <<  std::endl;
        if (_output_file != nullptr) {
          const CoordinateVector<> &pos = _source_positions.back();
          *_output_file << _total_time << "\t" << pos.x() << "\t" << pos.y()
                        << "\t" << pos.z() << "\t1\t"
                        << _source_indices.back() << "\t"
                        << _source_luminosities.back() << "\t"
                        << m_cur << "\t"
                        << "OSTAR\n";
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

        _last_sf = _total_time;
        // Rebuilding the radiation-copy hierarchy is only necessary when a
        // source was actually created. The bookkeeping update itself does not
        // change the source distribution.
        updated = updated || mass_generated > 0.;
        ++_number_of_updates;
    }









      if (_output_file != nullptr) {
        _output_file->flush();
      }

      



    return updated;
  }


  /**
   * @brief Move stellar sources with the external potential and local shear.
   *
   * Gas supplies the initial velocity at formation.  Thereafter the source is
   * treated as a collisionless particle: first the local gas gravitational
   * acceleration (external potential plus self gravity, when enabled), then
   * the same Coriolis/tidal rotation used for gas, followed by drift.
   */
  virtual void float_sources(
      DensitySubGridCreator< HydroDensitySubGrid > *grid_creator,
      const double timestep, const ExternalPotential *external_potential,
      const GalacticShearingBox *galactic_shearing_box,
      const CoordinateVector< bool > &periodicity) override {
    if (!_float_sources || timestep <= 0.) {
      return;
    }
    size_t i = 0;
    while (i < _source_positions.size()) {
      const bool was_inside =
          grid_creator->get_box().inside(_source_positions[i]);
      size_t old_subgrid = 0;
      if (was_inside) {
        old_subgrid = grid_creator->get_subgrid(_source_positions[i]).get_index();
      }
      if (grid_creator->get_box().inside(_source_positions[i])) {
        HydroDensitySubGrid &subgrid =
            *grid_creator->get_subgrid(_source_positions[i]);
        const auto cell = subgrid.get_hydro_cell(_source_positions[i]);
        _source_velocities[i] +=
            cell.get_hydro_variables().get_gravitational_acceleration() *
            timestep;
      } else if (external_potential != nullptr) {
        // A source can leave through a vertical boundary before disappearing.
        _source_velocities[i] += external_potential->get_acceleration(
                                    _source_positions[i]) *
                                timestep;
      }
      if (galactic_shearing_box != nullptr) {
        galactic_shearing_box->apply_to_source(_source_positions[i],
                                                _source_velocities[i],
                                                timestep);
      }
      _source_positions[i] += _source_velocities[i] * timestep;

      // Source particles obey the same ordinary periodic boundaries as the
      // gas.  (The Galactic shearing box still has no shearing remap.)
      const Box<> box = grid_creator->get_box();
      wrap_mixed_driving_source_position(_source_positions[i], box,
                                         periodicity);
      if (!grid_creator->get_box().inside(_source_positions[i])) {
        if (_log != nullptr) {
          _log->write_warning(
              "Removing MixedDriving source after it left the simulation box.");
        }
        _source_positions.erase(_source_positions.begin() + i);
        _source_velocities.erase(_source_velocities.begin() + i);
        _source_lifetimes.erase(_source_lifetimes.begin() + i);
        _source_luminosities.erase(_source_luminosities.begin() + i);
        _source_masses.erase(_source_masses.begin() + i);
        _spectrum_index.erase(_spectrum_index.begin() + i);
        if (i < _source_indices.size()) {
          _source_indices.erase(_source_indices.begin() + i);
        }
        _sources_changed = true;
      } else {
        if (was_inside &&
            old_subgrid !=
                grid_creator->get_subgrid(_source_positions[i]).get_index()) {
          _sources_changed = true;
        }
        ++i;
      }
    }
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
    restart_writer.write(_peak_fraction);
    restart_writer.write(_clustering_factor);
    restart_writer.write(_float_sources);
    restart_writer.write(init_running_mass);
    restart_writer.write(_num_sne);
    restart_writer.write(_holmes_time);
    restart_writer.write(_holmes_sh);
    restart_writer.write(_holmes_lum);
    restart_writer.write(_number_of_holmes);
    restart_writer.write(type1done);
    restart_writer.write(_total_time);
    restart_writer.write(_holmes_added);
    restart_writer.write(_last_sf);
    _random_generator.write_restart_file(restart_writer);
    {
      const auto size = _source_positions.size();
      restart_writer.write(size);
      for (std::vector< CoordinateVector<> >::size_type i = 0; i < size; ++i) {
        _source_positions[i].write_restart_file(restart_writer);
      }
    }
    {
      const auto size = _source_velocities.size();
      restart_writer.write(size);
      for (std::vector< CoordinateVector<> >::size_type i = 0; i < size; ++i) {
        _source_velocities[i].write_restart_file(restart_writer);
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
    {
      const auto size = _source_masses.size();
      restart_writer.write(size);
      for (const double mass : _source_masses) {
        restart_writer.write(mass);
      }
    }
    {
      const auto size = _snapshot_supernova_positions.size();
      restart_writer.write(size);
      for (const CoordinateVector<> &position :
           _snapshot_supernova_positions) {
        position.write_restart_file(restart_writer);
      }
      for (const double time : _snapshot_supernova_times) {
        restart_writer.write(time);
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
    }
    {
      const auto size = _source_indices.size();
      restart_writer.write(size);
      for (std::vector< uint_fast32_t >::size_type i = 0; i < size; ++i) {
        restart_writer.write(_source_indices[i]);
      }
    }
    restart_writer.write(_next_index);
  }

  /**
   * @brief Restart constructor.
   *
   * @param restart_reader Restart file to read from.
   */
  inline MixedDrivingPhotonSourceDistribution(RestartReader &restart_reader,
                                               const bool extended = false)
      : _star_formation_rate(restart_reader.read< double >()),
        _update_interval(restart_reader.read< double >()),
        _output_file(nullptr), _output_file2(nullptr), _number_of_updates(0),
        _next_index(0),
        _lum_adjust(restart_reader.read< double >()),
        _excess_mass(restart_reader.read<double>()),
        _scaleheight(restart_reader.read<double>()),
        _peak_fraction(restart_reader.read<double>()),
        _clustering_factor(restart_reader.read<double>()),
        _float_sources(restart_reader.read<bool>()),
        init_running_mass(restart_reader.read<double>()),
        _num_sne(restart_reader.read<uint_fast32_t>()),
        _holmes_time(restart_reader.read<double>()),
        _holmes_sh(restart_reader.read<double>()),
        _holmes_lum(restart_reader.read<double>()),
        _number_of_holmes(restart_reader.read<uint_fast32_t>()),
        _read_file(false), _filename(), _time(0.),
        type1done(restart_reader.read<int>()),
        _total_time(restart_reader.read<double>()),
        _holmes_added(restart_reader.read<bool>()),
        _last_sf(restart_reader.read<double>()),
        _random_generator(restart_reader), novahandler(nullptr), _log(nullptr) {

    {
      const std::vector< CoordinateVector<> >::size_type size =
          restart_reader.read< std::vector< CoordinateVector<> >::size_type >();
      _source_positions.resize(size);
      for (std::vector< CoordinateVector<> >::size_type i = 0; i < size; ++i) {
        _source_positions[i] = CoordinateVector<>(restart_reader);
      }
    }
    {
      const std::vector< CoordinateVector<> >::size_type size =
          restart_reader.read< std::vector< CoordinateVector<> >::size_type >();
      _source_velocities.resize(size);
      for (std::vector< CoordinateVector<> >::size_type i = 0; i < size; ++i) {
        _source_velocities[i] = CoordinateVector<>(restart_reader);
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
    if (extended) {
      const std::vector< double >::size_type size =
          restart_reader.read< std::vector< double >::size_type >();
      _source_masses.resize(size);
      for (double &mass : _source_masses) {
        mass = restart_reader.read< double >();
      }
      const std::vector< CoordinateVector<> >::size_type number_of_supernovae =
          restart_reader
              .read< std::vector< CoordinateVector<> >::size_type >();
      _snapshot_supernova_positions.resize(number_of_supernovae);
      _snapshot_supernova_times.resize(number_of_supernovae);
      for (CoordinateVector<> &position : _snapshot_supernova_positions) {
        position = CoordinateVector<>(restart_reader);
      }
      for (double &time : _snapshot_supernova_times) {
        time = restart_reader.read< double >();
      }
    } else {
      _source_masses.reserve(_source_luminosities.size());
      for (const double luminosity : _source_luminosities) {
        _source_masses.push_back(mass_from_luminosity(luminosity));
      }
    }
    _number_of_updates = restart_reader.read< uint_fast32_t >();
    const bool has_output = restart_reader.read< bool >();
    if (has_output) {
      const std::streampos filepos = restart_reader.read< std::streampos >();
      // truncate the original file to the size we were at
      if (truncate("MixedDriving_source_positions.txt", filepos) != 0) {
        cmac_error("Error while truncating output file!");
      }
      // now open the file in append mode
      _output_file = new std::ofstream("MixedDriving_source_positions.txt",
                                       std::ios_base::app);

      const std::streampos filepos2 = restart_reader.read< std::streampos >();
                                       // truncate the original file to the size we were at
      if (truncate("TotalLuminosity.txt", filepos2) != 0) {
              cmac_error("Error while truncating output file!");
            }
                                       // now open the file in append mode
      _output_file2 = new std::ofstream("TotalLuminosity.txt",
                                            std::ios_base::app);

      if (!extended) {
        const std::vector< uint_fast32_t >::size_type size =
            restart_reader.read< std::vector< uint_fast32_t >::size_type >();
        _source_indices.resize(size);
        for (std::vector< uint_fast32_t >::size_type i = 0; i < size; ++i) {
          _source_indices[i] = restart_reader.read< uint_fast32_t >();
        }
        _next_index = restart_reader.read< uint_fast32_t >();
      }
    }
    if (extended) {
      const std::vector< uint_fast32_t >::size_type size =
          restart_reader.read< std::vector< uint_fast32_t >::size_type >();
      _source_indices.resize(size);
      for (uint_fast32_t &index : _source_indices) {
        index = restart_reader.read< uint_fast32_t >();
      }
      _next_index = restart_reader.read< uint_fast32_t >();
    } else if (!has_output) {
      // Legacy restarts only stored source indices when text output was on.
      _source_indices.resize(_source_positions.size());
      for (uint_fast32_t i = 0; i < _source_indices.size(); ++i) {
        _source_indices[i] = i;
      }
      _next_index = _source_indices.size();
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

  initialize_spectra(nullptr);
  for (const double luminosity : _source_luminosities) {
    _spectrum_index.push_back(spectrum_index_from_luminosity(luminosity));
  }
  novahandler = new SupernovaHandler(_sne_energy);
  }
};

#endif // MIXEDDRIVINGPHOTONSOURCEDISTRIBUTION_HPP
