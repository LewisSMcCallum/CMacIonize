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
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>

namespace {

const std::string restart_prefix = "CMAC_SWIGGUM_RESTART_V2:";

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
    } else if (character != '\r' && character != '\n') {
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

std::vector< std::string > read_text_lines(const std::string &filename,
                                           const std::string &description) {
  std::vector< std::string > lines;
  if (filename.size() >= 3 &&
      filename.substr(filename.size() - 3) == ".gz") {
    const std::string command = "gzip -cd -- " + quote_for_shell(filename);
    FILE *pipe = popen(command.c_str(), "r");
    if (pipe == nullptr) {
      cmac_error("Unable to open compressed %s: %s", description.c_str(),
                 filename.c_str());
    }
    char buffer[65536];
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
      lines.push_back(buffer);
    }
    if (pclose(pipe) != 0) {
      cmac_error("Unable to decompress %s: %s", description.c_str(),
                 filename.c_str());
    }
  } else {
    std::ifstream file(filename.c_str());
    if (!file) {
      cmac_error("Unable to open %s: %s", description.c_str(),
                 filename.c_str());
    }
    std::string line;
    while (std::getline(file, line)) {
      lines.push_back(line);
    }
  }
  if (lines.empty()) {
    cmac_error("%s is empty: %s", description.c_str(), filename.c_str());
  }
  return lines;
}

std::string quote_csv(const std::string &value) {
  std::string result("\"");
  for (const char character : value) {
    result.push_back(character);
    if (character == '"') {
      result.push_back('"');
    }
  }
  result.push_back('"');
  return result;
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
      _log(log),
      _snap_clusters_to_gas(params.get_value< bool >(
          "PhotonSourceDistribution:snap clusters to gas", false)),
      _cluster_birth_filename(params.get_value< std::string >(
          "PhotonSourceDistribution:cluster birth filename",
          "cluster_births_july9_imf_sample0.csv.gz")),
      _snap_search_radius(params.get_physical_value< QUANTITY_LENGTH >(
          "PhotonSourceDistribution:snap search radius", "100. pc")),
      _snap_maximum_displacement(
          params.get_physical_value< QUANTITY_LENGTH >(
              "PhotonSourceDistribution:snap maximum displacement",
              "100. pc")),
      _snap_minimum_number_density(
          params.get_physical_value< QUANTITY_NUMBER_DENSITY >(
              "PhotonSourceDistribution:snap minimum number density",
              "1.0 cm^-3")),
      _snap_minimum_density_contrast(params.get_value< double >(
          "PhotonSourceDistribution:snap minimum density contrast", 2.)),
      _snap_minimum_neutral_fraction(params.get_value< double >(
          "PhotonSourceDistribution:snap minimum neutral fraction", 0.8)),
      _snap_maximum_temperature(
          params.get_physical_value< QUANTITY_TEMPERATURE >(
              "PhotonSourceDistribution:snap maximum temperature", "5000. K")),
      _snap_distance_scale(params.get_physical_value< QUANTITY_LENGTH >(
          "PhotonSourceDistribution:snap distance scale", "50. pc")),
      _snap_centroid_radius(params.get_physical_value< QUANTITY_LENGTH >(
          "PhotonSourceDistribution:snap centroid radius", "15.625 pc")),
      _snap_log_filename(params.get_value< std::string >(
          "PhotonSourceDistribution:snap log filename",
          "cluster_snap_events.csv")),
      _last_snap_subgrids_inspected(0), _supernova_handler(nullptr) {
  initialize_spectra();
  _supernova_handler = new SupernovaHandler(_supernova_energy);
  _supernova_handler->set_tigress_like_injection(params.get_value< bool >(
      "SupernovaHandler:TIGRESS like injection", true));
  load_file();
  if (_snap_clusters_to_gas) {
    validate_snap_configuration();
    load_cluster_birth_file();
    validate_snap_clusters();
    open_snap_log();
  }
  rebuild_sources();
}

SwiggumFilePhotonSourceDistribution::SwiggumFilePhotonSourceDistribution(
    RestartReader &restart_reader, Log *log)
    : _history_start_time(0.), _history_time(0.), _luminosity_adjustment(1.),
      _supernova_energy(0.), _log(log), _snap_clusters_to_gas(false),
      _snap_search_radius(0.), _snap_maximum_displacement(0.),
      _snap_minimum_number_density(0.), _snap_minimum_density_contrast(0.),
      _snap_minimum_neutral_fraction(0.), _snap_maximum_temperature(0.),
      _snap_distance_scale(0.), _snap_centroid_radius(0.),
      _last_snap_subgrids_inspected(0), _supernova_handler(nullptr) {
  const std::string stored_filename = restart_reader.read< std::string >();
  const bool version_two =
      stored_filename.compare(0, restart_prefix.size(), restart_prefix) == 0;
  _filename = version_two ? stored_filename.substr(restart_prefix.size())
                          : stored_filename;
  _box = read_box(restart_reader);
  _history_start_time = restart_reader.read< double >();
  _history_time = restart_reader.read< double >();
  _luminosity_adjustment = restart_reader.read< double >();
  _supernova_energy = restart_reader.read< double >();
  if (version_two) {
    _snap_clusters_to_gas = restart_reader.read< bool >();
    _cluster_birth_filename = restart_reader.read< std::string >();
    _snap_search_radius = restart_reader.read< double >();
    _snap_maximum_displacement = restart_reader.read< double >();
    _snap_minimum_number_density = restart_reader.read< double >();
    _snap_minimum_density_contrast = restart_reader.read< double >();
    _snap_minimum_neutral_fraction = restart_reader.read< double >();
    _snap_maximum_temperature = restart_reader.read< double >();
    _snap_distance_scale = restart_reader.read< double >();
    _snap_centroid_radius = restart_reader.read< double >();
    _snap_log_filename = restart_reader.read< std::string >();
  }
  initialize_spectra();
  _supernova_handler = new SupernovaHandler(_supernova_energy);
  load_file();
  if (_snap_clusters_to_gas) {
    validate_snap_configuration();
    load_cluster_birth_file();
    validate_snap_clusters();
    const size_t number_of_saved_clusters = restart_reader.read< size_t >();
    if (number_of_saved_clusters != _snap_clusters.size()) {
      cmac_error("Swiggum restart contains %zu snap clusters, but %zu were "
                 "loaded from %s.",
                 number_of_saved_clusters, _snap_clusters.size(),
                 _cluster_birth_filename.c_str());
    }
    for (size_t i = 0; i < number_of_saved_clusters; ++i) {
      const uint_fast32_t id = restart_reader.read< uint_fast32_t >();
      SnapCluster *cluster = nullptr;
      for (SnapCluster &candidate : _snap_clusters) {
        if (candidate.id == id) {
          cluster = &candidate;
          break;
        }
      }
      if (cluster == nullptr) {
        cmac_error("Unknown snap cluster ID %" PRIuFAST32 " in restart.", id);
      }
      cluster->placement_attempted = restart_reader.read< bool >();
      cluster->placement_successful = restart_reader.read< bool >();
      cluster->offset = CoordinateVector<>(restart_reader);
      cluster->status = restart_reader.read< std::string >();
    }
    open_snap_log();
  } else if (version_two) {
    const size_t number_of_saved_clusters = restart_reader.read< size_t >();
    if (number_of_saved_clusters != 0) {
      cmac_error("Disabled Swiggum cluster snapping has unexpected restart "
                 "state.");
    }
  }
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
  const std::vector< std::string > lines =
      read_text_lines(_filename, "Swiggum history");

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
  if (_snap_clusters_to_gas && columns.count("snap_cluster_id") == 0) {
    cmac_error("Snap-enabled Swiggum history %s has no snap_cluster_id column. "
               "Run prepare_snap_cluster_inputs.py first.",
               _filename.c_str());
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
    const uint_fast32_t snap_cluster_id =
        _snap_clusters_to_gas
            ? static_cast< uint_fast32_t >(
                  std::stoul(fields[columns["snap_cluster_id"]]))
            : std::numeric_limits< uint_fast32_t >::max();
    if (!star.track.empty() && star.snap_cluster_id != snap_cluster_id) {
      cmac_error("Star %" PRIuFAST64
                 " has multiple snap_cluster_id values in %s.",
                 id, _filename.c_str());
    }
    star.id = id;
    star.snap_cluster_id = snap_cluster_id;
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

void SwiggumFilePhotonSourceDistribution::load_cluster_birth_file() {
  const std::vector< std::string > lines =
      read_text_lines(_cluster_birth_filename, "Swiggum cluster birth file");
  const std::vector< std::string > header = split_csv_line(lines[0]);
  std::map< std::string, size_t > columns;
  for (size_t i = 0; i < header.size(); ++i) {
    columns[header[i]] = i;
  }
  const char *required[] = {
      "snap_cluster_id", "cluster_name",       "birth_time_myr",
      "birth_x_pc",     "birth_y_pc",          "birth_z_pc",
      "cluster_radius_pc", "number_of_sampled_stars"};
  for (const char *name : required) {
    if (columns.count(name) == 0) {
      cmac_error("Required column '%s' is absent from %s", name,
                 _cluster_birth_filename.c_str());
    }
  }

  const double pc = UnitConverter::to_SI< QUANTITY_LENGTH >(1., "pc");
  const double Myr = UnitConverter::to_SI< QUANTITY_TIME >(1., "Myr");
  std::map< uint_fast32_t, SnapCluster > clusters;
  std::map< std::string, uint_fast32_t > names;
  for (size_t iline = 1; iline < lines.size(); ++iline) {
    if (lines[iline].empty()) {
      continue;
    }
    const std::vector< std::string > fields = split_csv_line(lines[iline]);
    if (fields.size() != header.size()) {
      cmac_error("Malformed CSV row %zu in %s", iline + 1,
                 _cluster_birth_filename.c_str());
    }
    const uint_fast32_t id = static_cast< uint_fast32_t >(
        std::stoul(fields[columns["snap_cluster_id"]]));
    const std::string name = fields[columns["cluster_name"]];
    if (clusters.count(id) != 0 || names.count(name) != 0) {
      cmac_error("Duplicate snap cluster ID or name in %s: %" PRIuFAST32
                 " / %s.",
                 _cluster_birth_filename.c_str(), id, name.c_str());
    }
    SnapCluster cluster;
    cluster.id = id;
    cluster.name = name;
    cluster.birth_time =
        std::stod(fields[columns["birth_time_myr"]]) * Myr;
    cluster.nominal_birth_position = CoordinateVector<>(
        std::stod(fields[columns["birth_x_pc"]]) * pc,
        std::stod(fields[columns["birth_y_pc"]]) * pc,
        std::stod(fields[columns["birth_z_pc"]]) * pc);
    cluster.radius =
        std::stod(fields[columns["cluster_radius_pc"]]) * pc;
    cluster.placement_attempted = false;
    cluster.placement_successful = false;
    cluster.offset = CoordinateVector<>(0.);
    cluster.status = "not_yet_born";
    clusters[id] = cluster;
    names[name] = id;
  }
  _snap_clusters.clear();
  _snap_clusters.reserve(clusters.size());
  for (const auto &entry : clusters) {
    _snap_clusters.push_back(entry.second);
  }
}

void SwiggumFilePhotonSourceDistribution::validate_snap_configuration() const {
  if (_snap_search_radius < 0. || _snap_maximum_displacement < 0. ||
      _snap_minimum_number_density < 0. ||
      _snap_minimum_density_contrast < 0. ||
      _snap_minimum_neutral_fraction < 0. ||
      _snap_minimum_neutral_fraction > 1. ||
      _snap_maximum_temperature < 0. || !(_snap_distance_scale > 0.) ||
      _snap_centroid_radius < 0.) {
    cmac_error("Invalid dynamic Swiggum cluster-snapping configuration.");
  }
}

void SwiggumFilePhotonSourceDistribution::validate_snap_clusters() const {
  for (const SnapCluster &cluster : _snap_clusters) {
    if (!(cluster.birth_time > _history_start_time)) {
      cmac_error("Snap cluster %s is born at or before the configured history "
                 "start time.",
                 cluster.name.c_str());
    }
  }
  for (const Star &star : _stars) {
    if (get_snap_cluster(star.snap_cluster_id) == nullptr) {
      cmac_error("Star %" PRIuFAST64 " refers to missing snap cluster ID %"
                 PRIuFAST32 ".",
                 star.id, star.snap_cluster_id);
    }
  }
}

void SwiggumFilePhotonSourceDistribution::open_snap_log() {
  std::ifstream existing(_snap_log_filename.c_str(),
                         std::ios::binary | std::ios::ate);
  const bool write_header = !existing || existing.tellg() == 0;
  existing.close();
  std::ofstream output(_snap_log_filename.c_str(), std::ios::app);
  if (!output) {
    cmac_error("Unable to open cluster snap log: %s",
               _snap_log_filename.c_str());
  }
  if (write_header) {
    output << "cluster_id,cluster_name,birth_time_myr,nominal_x_pc,"
              "nominal_y_pc,nominal_z_pc,target_x_pc,target_y_pc,target_z_pc,"
              "dx_pc,dy_pc,dz_pc,displacement_pc,nominal_nH_cm3,"
              "target_nH_cm3,target_temperature_K,target_neutral_fraction,"
              "status,simulation_time_myr\n";
    output.flush();
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

const SwiggumFilePhotonSourceDistribution::SnapCluster *
SwiggumFilePhotonSourceDistribution::get_snap_cluster(
    const uint_fast32_t id) const {
  const auto cluster = std::lower_bound(
      _snap_clusters.begin(), _snap_clusters.end(), id,
      [](const SnapCluster &entry, const uint_fast32_t value) {
        return entry.id < value;
      });
  if (cluster == _snap_clusters.end() || cluster->id != id) {
    return nullptr;
  }
  return &(*cluster);
}

CoordinateVector<> SwiggumFilePhotonSourceDistribution::get_shifted_position(
    const Star &star, const double time) const {
  const CoordinateVector<> position = interpolate_position(star, time);
  if (!_snap_clusters_to_gas) {
    return position;
  }
  const SnapCluster *cluster = get_snap_cluster(star.snap_cluster_id);
  cmac_assert(cluster != nullptr);
  return position + cluster->offset;
}

void SwiggumFilePhotonSourceDistribution::log_snap_event(
    const SnapCluster &cluster, const CoordinateVector<> &target,
    const double nominal_number_density, const double target_number_density,
    const double target_temperature,
    const double target_neutral_fraction) const {
  const double pc = UnitConverter::to_SI< QUANTITY_LENGTH >(1., "pc");
  const double Myr = UnitConverter::to_SI< QUANTITY_TIME >(1., "Myr");
  const double per_cm3 = UnitConverter::to_SI< QUANTITY_NUMBER_DENSITY >(
      1., "cm^-3");
  const CoordinateVector<> displacement =
      target - cluster.nominal_birth_position;
  std::ofstream output(_snap_log_filename.c_str(), std::ios::app);
  if (!output) {
    cmac_error("Unable to append to cluster snap log: %s",
               _snap_log_filename.c_str());
  }
  output << std::setprecision(17) << cluster.id << ','
         << quote_csv(cluster.name) << ','
         << cluster.birth_time / Myr << ','
         << cluster.nominal_birth_position.x() / pc << ','
         << cluster.nominal_birth_position.y() / pc << ','
         << cluster.nominal_birth_position.z() / pc << ',' << target.x() / pc
         << ',' << target.y() / pc << ',' << target.z() / pc << ','
         << displacement.x() / pc << ',' << displacement.y() / pc << ','
         << displacement.z() / pc << ',' << displacement.norm() / pc << ','
         << nominal_number_density / per_cm3 << ','
         << target_number_density / per_cm3 << ',' << target_temperature << ','
         << target_neutral_fraction << ',' << cluster.status << ','
         << _history_time / Myr << '\n';
  output.flush();
}

void SwiggumFilePhotonSourceDistribution::place_cluster_in_current_gas(
    SnapCluster &cluster,
    DensitySubGridCreator< HydroDensitySubGrid > &grid_creator) {
  struct Candidate {
    CoordinateVector<> position;
    double number_density;
    double temperature;
    double neutral_fraction;
    double distance;
  };

  cluster.placement_attempted = true;
  cluster.placement_successful = false;
  cluster.offset = CoordinateVector<>(0.);
  _last_snap_subgrids_inspected = 0;
  const CoordinateVector<> nominal = cluster.nominal_birth_position;
  if (!_box.inside(nominal)) {
    cluster.status = "outside_box_at_birth";
    log_snap_event(cluster, nominal, 0., 0., 0., 0.);
    return;
  }
  if (grid_creator.is_spherical()) {
    cmac_error("Dynamic Swiggum cluster snapping currently requires a "
               "Cartesian task-based hydro grid.");
  }

  HydroDensitySubGrid &nominal_subgrid = *grid_creator.get_subgrid(nominal);
  const auto nominal_cell = nominal_subgrid.get_cell(nominal);
  const IonizationVariables &nominal_variables =
      nominal_cell.get_ionization_variables();
  const double nominal_density = nominal_variables.get_number_density();
  const double nominal_temperature = nominal_variables.get_temperature();
  const double nominal_neutral =
      nominal_variables.get_ionic_fraction(ION_H_n);
  const auto is_thermally_eligible =
      [this](const double density, const double temperature,
             const double neutral_fraction) {
    return density >= _snap_minimum_number_density &&
           neutral_fraction >= _snap_minimum_neutral_fraction &&
           temperature <= _snap_maximum_temperature;
  };
  if (is_thermally_eligible(nominal_density, nominal_temperature,
                            nominal_neutral)) {
    cluster.placement_successful = true;
    cluster.status = "already_in_gas";
    log_snap_event(cluster, nominal, nominal_density, nominal_density,
                   nominal_temperature, nominal_neutral);
    return;
  }

  const Box<> grid_box = grid_creator.get_box();
  const CoordinateVector< int_fast32_t > layout =
      grid_creator.get_subgrid_layout();
  const CoordinateVector<> subgrid_sides(
      grid_box.get_sides().x() / layout.x(),
      grid_box.get_sides().y() / layout.y(),
      grid_box.get_sides().z() / layout.z());
  const auto subgrid_index = [](const double coordinate, const double anchor,
                                const double side,
                                const int_fast32_t count) {
    return std::max< int_fast32_t >(
        0, std::min< int_fast32_t >(
               count - 1,
               static_cast< int_fast32_t >(std::floor((coordinate - anchor) /
                                                       side))));
  };
  const CoordinateVector<> radius_vector(_snap_search_radius);
  const CoordinateVector<> lower = nominal - radius_vector;
  const CoordinateVector<> upper = nominal + radius_vector;
  const int_fast32_t ix_min =
      subgrid_index(lower.x(), grid_box.get_anchor().x(), subgrid_sides.x(),
                    layout.x());
  const int_fast32_t iy_min =
      subgrid_index(lower.y(), grid_box.get_anchor().y(), subgrid_sides.y(),
                    layout.y());
  const int_fast32_t iz_min =
      subgrid_index(lower.z(), grid_box.get_anchor().z(), subgrid_sides.z(),
                    layout.z());
  const int_fast32_t ix_max =
      subgrid_index(upper.x(), grid_box.get_anchor().x(), subgrid_sides.x(),
                    layout.x());
  const int_fast32_t iy_max =
      subgrid_index(upper.y(), grid_box.get_anchor().y(), subgrid_sides.y(),
                    layout.y());
  const int_fast32_t iz_max =
      subgrid_index(upper.z(), grid_box.get_anchor().z(), subgrid_sides.z(),
                    layout.z());

  std::vector< Candidate > eligible;
  for (int_fast32_t ix = ix_min; ix <= ix_max; ++ix) {
    for (int_fast32_t iy = iy_min; iy <= iy_max; ++iy) {
      for (int_fast32_t iz = iz_min; iz <= iz_max; ++iz) {
        const size_t igrid =
            ix * layout.y() * layout.z() + iy * layout.z() + iz;
        cmac_assert(igrid < grid_creator.number_of_original_subgrids());
        ++_last_snap_subgrids_inspected;
        HydroDensitySubGrid &subgrid = *grid_creator.get_subgrid(igrid);
        for (auto cell = subgrid.begin(); cell != subgrid.end(); ++cell) {
          const CoordinateVector<> position = cell.get_cell_midpoint();
          const double distance = (position - nominal).norm();
          if (distance > _snap_search_radius) {
            continue;
          }
          const IonizationVariables &variables =
              cell.get_ionization_variables();
          const double density = variables.get_number_density();
          const double temperature = variables.get_temperature();
          const double neutral = variables.get_ionic_fraction(ION_H_n);
          if (is_thermally_eligible(density, temperature, neutral)) {
            eligible.push_back(
                {position, density, temperature, neutral, distance});
          }
        }
      }
    }
  }
  if (eligible.empty()) {
    cluster.status = "no_eligible_target";
    log_snap_event(cluster, nominal, nominal_density, 0., 0., 0.);
    return;
  }

  const double required_density =
      std::max(_snap_minimum_number_density,
               _snap_minimum_density_contrast * nominal_density);
  const Candidate *best = nullptr;
  double best_score = -std::numeric_limits< double >::infinity();
  for (const Candidate &candidate : eligible) {
    if (candidate.number_density <= required_density) {
      continue;
    }
    const double scaled_distance = candidate.distance / _snap_distance_scale;
    const double score =
        std::log(candidate.number_density / required_density) -
        0.5 * scaled_distance * scaled_distance;
    if (score > best_score) {
      best_score = score;
      best = &candidate;
    }
  }
  if (best == nullptr) {
    cluster.status = "density_contrast_not_met";
    log_snap_event(cluster, nominal, nominal_density, 0., 0., 0.);
    return;
  }

  CoordinateVector<> weighted_position;
  double centroid_weight = 0.;
  for (const Candidate &candidate : eligible) {
    if (candidate.number_density > required_density &&
        (candidate.position - best->position).norm() <=
            _snap_centroid_radius) {
      weighted_position += candidate.number_density * candidate.position;
      centroid_weight += candidate.number_density;
    }
  }
  cmac_assert(centroid_weight > 0.);
  const CoordinateVector<> target = weighted_position / centroid_weight;
  const CoordinateVector<> offset = target - nominal;
  if (offset.norm() > _snap_maximum_displacement) {
    cluster.status = "maximum_displacement_exceeded";
    log_snap_event(cluster, target, nominal_density, best->number_density,
                   best->temperature, best->neutral_fraction);
    return;
  }
  cluster.offset = offset;
  cluster.placement_successful = true;
  cluster.status = "snapped";
  log_snap_event(cluster, target, nominal_density, best->number_density,
                 best->temperature, best->neutral_fraction);
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
    const CoordinateVector<> position =
        get_shifted_position(star, _history_time);
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

double SwiggumFilePhotonSourceDistribution::get_photon_frequency_weighted(
    RandomGenerator &random, photonsourcenumber_t index,
    double uniform_fraction, double &weight) {
  return _spectra[_spectrum_indices[index]]->get_random_frequency_weighted(
      random, uniform_fraction, weight);
}

bool SwiggumFilePhotonSourceDistribution::update(
    DensitySubGridCreator< HydroDensitySubGrid > *grid_creator,
    const double actual_timestep) {
  const double old_time = _history_time;
  const std::vector< uint_fast64_t > old_ids = _source_ids;
  const std::vector< CoordinateVector<> > old_positions = _source_positions;
  _history_time += actual_timestep;

  if (_snap_clusters_to_gas) {
    if (grid_creator == nullptr) {
      cmac_error("Snap-enabled Swiggum source update requires the live hydro "
                 "grid.");
    }
    for (SnapCluster &cluster : _snap_clusters) {
      if (!cluster.placement_attempted && old_time < cluster.birth_time &&
          cluster.birth_time <= _history_time) {
        place_cluster_in_current_gas(cluster, *grid_creator);
      }
    }
  }

  for (const Star &star : _stars) {
    if (old_time < star.death_time && star.death_time <= _history_time) {
      const CoordinateVector<> position =
          get_shifted_position(star, star.death_time);
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
  restart_writer.write(restart_prefix + _filename);
  _box.get_anchor().write_restart_file(restart_writer);
  _box.get_sides().write_restart_file(restart_writer);
  restart_writer.write(_history_start_time);
  restart_writer.write(_history_time);
  restart_writer.write(_luminosity_adjustment);
  restart_writer.write(_supernova_energy);
  restart_writer.write(_snap_clusters_to_gas);
  restart_writer.write(_cluster_birth_filename);
  restart_writer.write(_snap_search_radius);
  restart_writer.write(_snap_maximum_displacement);
  restart_writer.write(_snap_minimum_number_density);
  restart_writer.write(_snap_minimum_density_contrast);
  restart_writer.write(_snap_minimum_neutral_fraction);
  restart_writer.write(_snap_maximum_temperature);
  restart_writer.write(_snap_distance_scale);
  restart_writer.write(_snap_centroid_radius);
  restart_writer.write(_snap_log_filename);
  restart_writer.write(_snap_clusters.size());
  for (const SnapCluster &cluster : _snap_clusters) {
    restart_writer.write(cluster.id);
    restart_writer.write(cluster.placement_attempted);
    restart_writer.write(cluster.placement_successful);
    cluster.offset.write_restart_file(restart_writer);
    restart_writer.write(cluster.status);
  }
}
