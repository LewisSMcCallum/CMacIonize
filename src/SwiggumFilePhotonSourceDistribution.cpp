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
 * @file SwiggumFilePhotonSourceDistribution.cpp
 *
 * @brief Time-dependent stellar sources read from a Swiggum history file.
 */

#include "SwiggumFilePhotonSourceDistribution.hpp"

#include "Error.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "RandomGenerator.hpp"
#include "RestartReader.hpp"
#include "RestartWriter.hpp"
#include "SupernovaHandler.hpp"
#include "UnitConverter.hpp"
#include "WMBasicPhotonSourceSpectrum.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>

namespace {

std::vector< std::string > split_csv_line(const std::string &line) {
  std::vector< std::string > fields;
  std::string field;
  bool quoted = false;
  for (size_t i = 0; i < line.size(); ++i) {
    const char character = line[i];
    if (character == '"') {
      if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
        field.push_back('"');
        ++i;
      } else {
        quoted = !quoted;
      }
    } else if (character == ',' && !quoted) {
      fields.push_back(field);
      field.clear();
    } else if (character != '\r') {
      field.push_back(character);
    }
  }
  fields.push_back(field);
  return fields;
}

std::string quote_for_shell(const std::string &value) {
  std::string quoted("'");
  for (const char character : value) {
    if (character == '\'') {
      quoted += "'\\''";
    } else {
      quoted.push_back(character);
    }
  }
  quoted.push_back('\'');
  return quoted;
}

Box<> read_box(RestartReader &restart_reader) {
  const CoordinateVector<> anchor(restart_reader);
  const CoordinateVector<> sides(restart_reader);
  return Box<>(anchor, sides);
}

} // namespace

void SwiggumFilePhotonSourceDistribution::
    set_tigress_like_supernova_injection(const bool value) {
  _supernova_handler->set_tigress_like_injection(value);
}

SwiggumFilePhotonSourceDistribution::SwiggumFilePhotonSourceDistribution(
    ParameterFile &params, Log *log)
    : _filename(params.get_value< std::string >(
          "PhotonSourceDistribution:filename", "swiggum_history.csv.gz")),
      _box(params.get_physical_vector< QUANTITY_LENGTH >(
               "SimulationBox:anchor", "[-500. pc, -500. pc, -500. pc]"),
           params.get_physical_vector< QUANTITY_LENGTH >(
               "SimulationBox:sides", "[1. kpc, 1. kpc, 1. kpc]")),
      _history_start_time(params.get_physical_value< QUANTITY_TIME >(
          "PhotonSourceDistribution:history start time", "-100. Myr")),
      _history_time(_history_start_time),
      _luminosity_adjustment(params.get_value< double >(
          "PhotonSourceDistribution:luminosity adjust", 1.)),
      _supernova_energy(params.get_physical_value< QUANTITY_ENERGY >(
          "PhotonSourceDistribution:supernova energy", "1.e51 erg")),
      _log(log), _supernova_handler(nullptr) {
  initialize_spectra();
  _supernova_handler = new SupernovaHandler(_supernova_energy);
  _supernova_handler->set_tigress_like_injection(params.get_value< bool >(
      "SupernovaHandler:TIGRESS like injection", true));
  load_file();
  rebuild_sources();
}

SwiggumFilePhotonSourceDistribution::SwiggumFilePhotonSourceDistribution(
    RestartReader &restart_reader, Log *log)
    : _filename(restart_reader.read< std::string >()),
      _box(read_box(restart_reader)),
      _history_start_time(restart_reader.read< double >()),
      _history_time(restart_reader.read< double >()),
      _luminosity_adjustment(restart_reader.read< double >()),
      _supernova_energy(restart_reader.read< double >()), _log(log),
      _supernova_handler(nullptr) {
  initialize_spectra();
  _supernova_handler = new SupernovaHandler(_supernova_energy);
  load_file();
  rebuild_sources();
}

SwiggumFilePhotonSourceDistribution::~SwiggumFilePhotonSourceDistribution() {
  for (PhotonSourceSpectrum *spectrum : _spectra) {
    delete spectrum;
  }
  delete _supernova_handler;
}

void SwiggumFilePhotonSourceDistribution::initialize_spectra() {
  const double temperatures[] = {32000., 34000., 34000., 35000., 36000.,
                                 37000., 39000., 39000., 40000., 41000.,
                                 42000., 43000., 44000., 45000.};
  const double gravities[] = {25., 25., 25., 40., 25., 25., 25.,
                              25., 25., 40., 40., 40., 40., 40.};
  for (size_t i = 0; i < sizeof(temperatures) / sizeof(temperatures[0]); ++i) {
    _spectra.push_back(
        new WMBasicPhotonSourceSpectrum(temperatures[i], gravities[i], _log));
  }
}

void SwiggumFilePhotonSourceDistribution::load_file() {
  std::vector< std::string > lines;
  if (_filename.size() >= 3 &&
      _filename.substr(_filename.size() - 3) == ".gz") {
    const std::string command = "gzip -cd -- " + quote_for_shell(_filename);
    FILE *pipe = popen(command.c_str(), "r");
    if (pipe == nullptr) {
      cmac_error("Unable to open compressed Swiggum history: %s",
                 _filename.c_str());
    }
    char buffer[65536];
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
      lines.push_back(buffer);
    }
    if (pclose(pipe) != 0) {
      cmac_error("Unable to decompress Swiggum history: %s",
                 _filename.c_str());
    }
  } else {
    std::ifstream file(_filename.c_str());
    if (!file) {
      cmac_error("Unable to open Swiggum history: %s", _filename.c_str());
    }
    std::string line;
    while (std::getline(file, line)) {
      lines.push_back(line);
    }
  }

  if (lines.empty()) {
    cmac_error("Swiggum history file is empty: %s", _filename.c_str());
  }

  const std::vector< std::string > header = split_csv_line(lines[0]);
  std::map< std::string, size_t > columns;
  for (size_t i = 0; i < header.size(); ++i) {
    columns[header[i]] = i;
  }
  const char *required[] = {"trajectory_time_myr", "SNTime",
                            "stellar_birth_time_myr", "star_id", "X", "Y",
                            "Z", "Stellar_Mass"};
  for (const char *name : required) {
    if (columns.count(name) == 0) {
      cmac_error("Required column '%s' is absent from %s", name,
                 _filename.c_str());
    }
  }

  const double pc = UnitConverter::to_SI< QUANTITY_LENGTH >(1., "pc");
  const double Myr = UnitConverter::to_SI< QUANTITY_TIME >(1., "Myr");
  std::map< uint_fast64_t, Star > stars;
  for (size_t iline = 1; iline < lines.size(); ++iline) {
    if (lines[iline].empty()) {
      continue;
    }
    const std::vector< std::string > fields = split_csv_line(lines[iline]);
    if (fields.size() != header.size()) {
      cmac_error("Malformed CSV row %zu in %s", iline + 1,
                 _filename.c_str());
    }
    const uint_fast64_t id =
        static_cast< uint_fast64_t >(std::stoull(fields[columns["star_id"]]));
    Star &star = stars[id];
    star.id = id;
    star.birth_time =
        std::stod(fields[columns["stellar_birth_time_myr"]]) * Myr;
    star.death_time = std::stod(fields[columns["SNTime"]]) * Myr;
    star.mass = std::stod(fields[columns["Stellar_Mass"]]);
    TrackPoint point;
    point.time = std::stod(fields[columns["trajectory_time_myr"]]) * Myr;
    point.position =
        CoordinateVector<>(std::stod(fields[columns["X"]]) * pc,
                           std::stod(fields[columns["Y"]]) * pc,
                           std::stod(fields[columns["Z"]]) * pc);
    star.track.push_back(point);
  }

  _stars.reserve(stars.size());
  for (auto &entry : stars) {
    Star &star = entry.second;
    std::sort(star.track.begin(), star.track.end(),
              [](const TrackPoint &left, const TrackPoint &right) {
                return left.time < right.time;
              });
    star.track.erase(
        std::unique(star.track.begin(), star.track.end(),
                    [](const TrackPoint &left, const TrackPoint &right) {
                      return left.time == right.time;
                    }),
        star.track.end());
    _stars.push_back(star);
  }

  if (_log) {
    _log->write_status("Loaded ", _stars.size(),
                       " stellar trajectories from ", _filename, ".");
  }
}

CoordinateVector<> SwiggumFilePhotonSourceDistribution::interpolate_position(
    const Star &star, const double time) const {
  cmac_assert(!star.track.empty());
  if (time <= star.track.front().time) {
    return star.track.front().position;
  }
  if (time >= star.track.back().time) {
    return star.track.back().position;
  }
  const auto upper = std::lower_bound(
      star.track.begin(), star.track.end(), time,
      [](const TrackPoint &point, const double value) {
        return point.time < value;
      });
  const TrackPoint &right = *upper;
  const TrackPoint &left = *(upper - 1);
  const double fraction = (time - left.time) / (right.time - left.time);
  return left.position + fraction * (right.position - left.position);
}

double SwiggumFilePhotonSourceDistribution::luminosity_from_mass(
    const double mass) const {
  const double masses[] = {57.95, 46.94, 38.08, 34.39, 30.98, 28.0,
                           25.29, 22.90, 20.76, 18.80, 17.08, 15.55};
  const double log_luminosities[] = {49.64, 49.44, 49.22, 49.10, 48.99, 48.88,
                                     48.75, 48.61, 48.44, 48.27, 48.06, 47.88};
  const size_t count = sizeof(masses) / sizeof(masses[0]);
  if (mass < masses[count - 1]) {
    return 0.;
  }
  double log_luminosity = log_luminosities[0];
  if (mass < masses[0]) {
    for (size_t i = 0; i + 1 < count; ++i) {
      if (masses[i] >= mass && mass >= masses[i + 1]) {
        const double fraction =
            (mass - masses[i]) / (masses[i + 1] - masses[i]);
        log_luminosity = log_luminosities[i] +
                         fraction *
                             (log_luminosities[i + 1] - log_luminosities[i]);
        break;
      }
    }
  }
  return _luminosity_adjustment * std::pow(10., log_luminosity);
}

size_t SwiggumFilePhotonSourceDistribution::spectrum_index_from_mass(
    const double mass) const {
  const double masses[] = {57.95, 46.94, 38.08, 34.39, 30.98, 28.0,
                           25.29, 22.90, 20.76, 18.80, 17.08, 15.55};
  const double temperatures[] = {44852., 42857., 40862., 39865., 38867., 37870.,
                                 36872., 35874., 34877., 33879., 32882., 31884.};
  const double available[] = {32000., 34000., 34000., 35000., 36000., 37000.,
                              39000., 39000., 40000., 41000., 42000., 43000.,
                              44000., 45000.};
  const size_t count = sizeof(masses) / sizeof(masses[0]);
  double temperature = temperatures[0];
  if (mass <= masses[count - 1]) {
    temperature = temperatures[count - 1];
  } else if (mass < masses[0]) {
    for (size_t i = 0; i + 1 < count; ++i) {
      if (masses[i] >= mass && mass >= masses[i + 1]) {
        const double fraction =
            (mass - masses[i]) / (masses[i + 1] - masses[i]);
        temperature = temperatures[i] +
                      fraction * (temperatures[i + 1] - temperatures[i]);
        break;
      }
    }
  }
  size_t closest = 0;
  for (size_t i = 1; i < sizeof(available) / sizeof(available[0]); ++i) {
    if (std::abs(available[i] - temperature) <
        std::abs(available[closest] - temperature)) {
      closest = i;
    }
  }
  return closest;
}

void SwiggumFilePhotonSourceDistribution::rebuild_sources() {
  _source_positions.clear();
  _source_luminosities.clear();
  _source_ids.clear();
  _spectrum_indices.clear();
  for (const Star &star : _stars) {
    if (_history_time < star.birth_time || _history_time >= star.death_time) {
      continue;
    }
    const double luminosity = luminosity_from_mass(star.mass);
    if (luminosity <= 0.) {
      continue;
    }
    const CoordinateVector<> position = interpolate_position(star, _history_time);
    if (_box.inside(position)) {
      _source_positions.push_back(position);
      _source_luminosities.push_back(luminosity);
      _source_ids.push_back(star.id);
      _spectrum_indices.push_back(spectrum_index_from_mass(star.mass));
    }
  }
}

photonsourcenumber_t
SwiggumFilePhotonSourceDistribution::get_number_of_sources() const {
  return _source_positions.size();
}

CoordinateVector<> SwiggumFilePhotonSourceDistribution::get_position(
    photonsourcenumber_t index) {
  return _source_positions[index];
}

double SwiggumFilePhotonSourceDistribution::get_weight(
    photonsourcenumber_t index) const {
  return _source_luminosities[index] / get_total_luminosity();
}

double SwiggumFilePhotonSourceDistribution::get_total_luminosity() const {
  double total = 0.;
  for (const double luminosity : _source_luminosities) {
    total += luminosity;
  }
  return total;
}

double SwiggumFilePhotonSourceDistribution::get_photon_frequency(
    RandomGenerator &random_generator, photonsourcenumber_t index) {
  return _spectra[_spectrum_indices[index]]->get_random_frequency(
      random_generator, 0.);
}

bool SwiggumFilePhotonSourceDistribution::update(
    DensitySubGridCreator< HydroDensitySubGrid > *grid_creator,
    const double actual_timestep) {
  const double old_time = _history_time;
  const std::vector< uint_fast64_t > old_ids = _source_ids;
  const std::vector< CoordinateVector<> > old_positions = _source_positions;
  _history_time += actual_timestep;

  for (const Star &star : _stars) {
    if (old_time < star.death_time && star.death_time <= _history_time) {
      const CoordinateVector<> position =
          interpolate_position(star, star.death_time);
      if (_box.inside(position)) {
        _feedback_positions.push_back(position);
      }
    }
  }
  rebuild_sources();

  if (old_ids != _source_ids || old_positions.size() != _source_positions.size()) {
    return true;
  }
  for (size_t i = 0; i < old_positions.size(); ++i) {
    if (grid_creator->get_subgrid(old_positions[i]).get_index() !=
        grid_creator->get_subgrid(_source_positions[i]).get_index()) {
      return true;
    }
  }
  return false;
}

bool SwiggumFilePhotonSourceDistribution::do_stellar_feedback(
    const double current_time) const {
  (void)current_time;
  return !_feedback_positions.empty();
}

void SwiggumFilePhotonSourceDistribution::get_sne_radii(
    DensitySubGridCreator< HydroDensitySubGrid > &grid_creator) {
  for (const CoordinateVector<> &position : _feedback_positions) {
    double injection_radius, sedov_taylor_radius, density, cell_count;
    std::tie(injection_radius, sedov_taylor_radius, density, cell_count) =
        _supernova_handler->get_r_inj(&grid_creator, position);
    _injection_radii.push_back(injection_radius);
    _sedov_taylor_radii.push_back(sedov_taylor_radius);
    _injection_densities.push_back(density);
    _injection_cell_counts.push_back(cell_count);
  }
}

void SwiggumFilePhotonSourceDistribution::add_stellar_feedback(
    HydroDensitySubGrid &subgrid, Hydro &hydro) {
  for (size_t i = 0; i < _feedback_positions.size(); ++i) {
    _supernova_handler->inject_sne(
        subgrid, hydro, _feedback_positions[i], _injection_radii[i],
        _sedov_taylor_radii[i], _injection_densities[i],
        _injection_cell_counts[i]);
  }
}

void SwiggumFilePhotonSourceDistribution::done_stellar_feedback() {
  _feedback_positions.clear();
  _injection_radii.clear();
  _sedov_taylor_radii.clear();
  _injection_densities.clear();
  _injection_cell_counts.clear();
}

void SwiggumFilePhotonSourceDistribution::write_restart_file(
    RestartWriter &restart_writer) const {
  restart_writer.write(_filename);
  _box.get_anchor().write_restart_file(restart_writer);
  _box.get_sides().write_restart_file(restart_writer);
  restart_writer.write(_history_start_time);
  restart_writer.write(_history_time);
  restart_writer.write(_luminosity_adjustment);
  restart_writer.write(_supernova_energy);
}
