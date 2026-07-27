/*******************************************************************************
 * This file is part of CMacIonize
 ******************************************************************************/

/**
 * @file testTaskBasedRadiationRestart.cpp
 *
 * @brief Regression test for the task-based RHD last-radiation restart state.
 */

#include "Assert.hpp"
#include "TaskBasedRadiationRestart.hpp"

#include <cstdio>

int main(int argc, char **argv) {
  const double radiation_interval = 1.577e12;
  const double hydro_timestep = 6.01578e9;
  const uint_fast32_t next_radiation = 82;
  const double restart_current_time = 1.27751e14;

  // This is the midsfr legacy-restart case. The old behaviour returned zero,
  // causing a roughly 4.1 Myr chemistry update instead of roughly 0.05 Myr.
  const double reconstructed =
      TaskBasedRadiationRestart::reconstruct_last_radiation_time(
          restart_current_time, hydro_timestep, next_radiation,
          radiation_interval);
  assert_condition(reconstructed > 1.277e14);
  assert_condition(restart_current_time - reconstructed <
                   1.01 * radiation_interval);

  // An old restart ends immediately after the legacy fields. Verify that an
  // empty appendix takes the reconstruction path used by the midsfr dump.
  const char *legacy_filename = "test_taskbased_radiation_legacy.dump";
  {
    RestartWriter writer(legacy_filename);
  }
  {
    RestartReader reader(legacy_filename);
    const double restored = TaskBasedRadiationRestart::read_last_radiation_time(
        reader, restart_current_time, hydro_timestep, next_radiation,
        radiation_interval);
    assert_condition(restored == reconstructed);
  }
  std::remove(legacy_filename);

  // Verify that an exact value appended by new versions takes precedence.
  const char *filename = "test_taskbased_radiation_restart.dump";
  const double exact_last_radiation_time = 1.277432345e14;
  {
    RestartWriter writer(filename);
    TaskBasedRadiationRestart::write_last_radiation_time(
        writer, exact_last_radiation_time);
  }
  {
    RestartReader reader(filename);
    const double restored = TaskBasedRadiationRestart::read_last_radiation_time(
        reader, restart_current_time, hydro_timestep, next_radiation,
        radiation_interval);
    assert_condition(restored == exact_last_radiation_time);
  }
  std::remove(filename);

  // Radiation on every hydro step should resume from the completed step.
  assert_condition(
      TaskBasedRadiationRestart::reconstruct_last_radiation_time(
          100., 2., 10, -1.) == 98.);

  return 0;
}
