/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2026 Lewis McCallum
 *
 * CMacIonize is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 ******************************************************************************/

/**
 * @file SwiggumFilePhotonSourceDistribution.hpp
 *
 * @brief Time-dependent stellar sources read from a Swiggum history file.
 */
#ifndef SWIGGUMFILEPHOTONSOURCEDISTRIBUTION_HPP
#define SWIGGUMFILEPHOTONSOURCEDISTRIBUTION_HPP

#include "Box.hpp"
#include "PhotonSourceDistribution.hpp"

#include <cstdint>
#include <string>
#include <vector>

class Log;
class ParameterFile;
class PhotonSourceSpectrum;
class RestartReader;
class SupernovaHandler;

/**
 * @brief Stellar photon sources and supernovae following tabulated tracks.
 */
class SwiggumFilePhotonSourceDistribution : public PhotonSourceDistribution {
private:
  struct TrackPoint {
    double time;
    CoordinateVector<> position;
  };

  struct Star {
    uint_fast64_t id;
    double birth_time;
    double death_time;
    double mass;
    std::vector< TrackPoint > track;
  };

  std::string _filename;
  Box<> _box;
  double _history_start_time;
  double _history_time;
  double _luminosity_adjustment;
  double _supernova_energy;
  Log *_log;

  std::vector< Star > _stars;
  std::vector< CoordinateVector<> > _source_positions;
  std::vector< double > _source_luminosities;
  std::vector< uint_fast64_t > _source_ids;
  std::vector< size_t > _spectrum_indices;
  std::vector< PhotonSourceSpectrum * > _spectra;

  std::vector< CoordinateVector<> > _feedback_positions;
  std::vector< double > _injection_radii;
  std::vector< double > _sedov_taylor_radii;
  std::vector< double > _injection_cell_counts;
  std::vector< double > _injection_densities;
  SupernovaHandler *_supernova_handler;

  void load_file();
  void initialize_spectra();
  void rebuild_sources();
  CoordinateVector<> interpolate_position(const Star &star,
                                           const double time) const;
  double luminosity_from_mass(const double mass) const;
  size_t spectrum_index_from_mass(const double mass) const;

public:
  SwiggumFilePhotonSourceDistribution(ParameterFile &params, Log *log = nullptr);
  SwiggumFilePhotonSourceDistribution(RestartReader &restart_reader,
                                      Log *log = nullptr);
  virtual ~SwiggumFilePhotonSourceDistribution();

  virtual photonsourcenumber_t get_number_of_sources() const;
  virtual CoordinateVector<> get_position(photonsourcenumber_t index);
  virtual double get_weight(photonsourcenumber_t index) const;
  virtual double get_total_luminosity() const;
  virtual double get_photon_frequency(RandomGenerator &random_generator,
                                      photonsourcenumber_t index);

  virtual bool update(
      DensitySubGridCreator< HydroDensitySubGrid > *grid_creator,
      double actual_timestep);

  virtual bool do_stellar_feedback(const double current_time) const;
  virtual void get_sne_radii(
      DensitySubGridCreator< HydroDensitySubGrid > &grid_creator);
  virtual void add_stellar_feedback(HydroDensitySubGrid &subgrid, Hydro &hydro);
  virtual void done_stellar_feedback();

  virtual void write_restart_file(RestartWriter &restart_writer) const;
};

#endif // SWIGGUMFILEPHOTONSOURCEDISTRIBUTION_HPP
