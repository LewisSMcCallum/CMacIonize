/*******************************************************************************
 * This file is part of CMacIonize
 ******************************************************************************/

/**
 * @file TaskBasedRadiationRestart.hpp
 *
 * @brief Restart helpers for task-based radiation hydrodynamics.
 */
#ifndef TASKBASEDRADIATIONRESTART_HPP
#define TASKBASEDRADIATIONRESTART_HPP

#include "RestartReader.hpp"
#include "RestartWriter.hpp"

#include <algorithm>
#include <cstdint>

namespace TaskBasedRadiationRestart {

/*! @brief Marker for the optional restart appendix. */
static constexpr uint64_t LAST_RADIATION_TIME_MARKER =
    UINT64_C(0x434d414352414454);

/**
 * @brief Reconstruct the last radiation time from a legacy restart.
 *
 * Older task-based RHD restart files stored the radiation-step counter, but
 * not the actual time of the last radiation update. The returned value is
 * bounded by the end of the last completed hydro step and is accurate to
 * within one hydro timestep.
 *
 * @param current_time End time of the next hydro step.
 * @param actual_timestep Duration of the next hydro step.
 * @param hydro_lastrad Counter for the next scheduled radiation threshold.
 * @param hydro_radtime Radiation update interval.
 * @return Reconstructed time of the last radiation update.
 */
inline double reconstruct_last_radiation_time(
    const double current_time, const double actual_timestep,
    const uint_fast32_t hydro_lastrad, const double hydro_radtime) {
  if (hydro_lastrad == 0) {
    return 0.;
  }

  const double last_completed_time =
      std::max(0., current_time - actual_timestep);
  if (hydro_radtime < 0.) {
    return last_completed_time;
  }

  const double threshold_time =
      static_cast< double >(hydro_lastrad - 1) * hydro_radtime;
  return std::min(last_completed_time, threshold_time + actual_timestep);
}

/**
 * @brief Append the exact last radiation time to a restart file.
 *
 * Appending rather than inserting this value keeps new restart files readable
 * by older executables and keeps old restart files readable by new ones.
 *
 * @param restart_writer Restart writer.
 * @param last_radiation_time Exact time of the last radiation update.
 */
inline void write_last_radiation_time(RestartWriter &restart_writer,
                                      const double last_radiation_time) {
  restart_writer.write(LAST_RADIATION_TIME_MARKER);
  restart_writer.write(last_radiation_time);
}

/**
 * @brief Read the exact last radiation time, or reconstruct it for an old dump.
 *
 * @param restart_reader Restart reader positioned after the legacy fields.
 * @param current_time End time of the next hydro step.
 * @param actual_timestep Duration of the next hydro step.
 * @param hydro_lastrad Counter for the next scheduled radiation threshold.
 * @param hydro_radtime Radiation update interval.
 * @return Exact or reconstructed time of the last radiation update.
 */
inline double read_last_radiation_time(
    RestartReader &restart_reader, const double current_time,
    const double actual_timestep, const uint_fast32_t hydro_lastrad,
    const double hydro_radtime) {
  if (restart_reader.has_more_data()) {
    const uint64_t marker = restart_reader.read< uint64_t >();
    if (marker == LAST_RADIATION_TIME_MARKER &&
        restart_reader.has_more_data()) {
      return restart_reader.read< double >();
    }
  }
  return reconstruct_last_radiation_time(
      current_time, actual_timestep, hydro_lastrad, hydro_radtime);
}

} // namespace TaskBasedRadiationRestart

#endif // TASKBASEDRADIATIONRESTART_HPP
