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
 * @file TaskBasedIonizationSimulation.hpp
 *
 * @brief Ionization radiative transfer simulation using a task-based parallel
 * algorithm.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef TASKBASEDIONIZATIONSIMULATION_HPP
#define TASKBASEDIONIZATIONSIMULATION_HPP

#include "AbundanceModel.hpp"
#include "Abundances.hpp"
#include "ChargeTransferRates.hpp"
#include "CollisionalRates.hpp"
#include "LineCoolingData.hpp"
#include "MemoryLogger.hpp"
#include "ParameterFile.hpp"
#include "RandomGenerator.hpp"
#include "SimulationBox.hpp"
#include "Task.hpp"
#include "ThreadSafeVector.hpp"
#include "TimeLogger.hpp"

#include <vector>

class ContinuousPhotonSource;
class CrossSections;
class DensityFunction;
class DensityGridWriter;
class DensitySubGrid;
template < class _subgrid_type_ > class DensitySubGridCreator;
class DiffuseReemissionHandler;
class MemorySpace;
class PhotonSourceDistribution;
class PhotonSourceSpectrum;
class RecombinationRates;
class TaskQueue;
class TemperatureCalculator;
class TrackerManager;

/**
 * @brief Ionization radiative transfer simulation using a task-based parallel
 * algorithm.
 */
class TaskBasedIonizationSimulation {
private:
  /*! @brief PhotonPacket buffers. */
  MemorySpace *_buffers;

  /*! @brief Queues per thread. */
  std::vector< TaskQueue * > _queues;

  /*! @brief General shared queue. */
  TaskQueue *_shared_queue;

  /*! @brief Task space. */
  ThreadSafeVector< Task > *_tasks;

  /*! @brief Random number generator per thread. */
  std::vector< RandomGenerator > _random_generators;

  /*! @brief ParameterFile containing the run parameters. */
  ParameterFile _parameter_file;

  /*! @brief Number of iterations of the ray tracing loop. */
  const uint_fast32_t _number_of_iterations;

  /*! @brief Number of photons used in every iteration of the ray tracing
   *  loop. */
  uint_fast64_t _number_of_photons;

  /*! @brief Copy level for subgrids that contain a source. */
  const uint_fast8_t _source_copy_level;

  /*! @brief Simulation box (in m). */
  SimulationBox _simulation_box;

  /*! @brief Grid creator. */
  DensitySubGridCreator< DensitySubGrid > *_grid_creator;

  /*! @brief DensityFunction that sets the density field. */
  DensityFunction *_density_function;

  /*! @brief DensityGridWriter used for snapshots. */
  DensityGridWriter *_density_grid_writer;

  /*! @brief PhotonSourceDistribution specifying the positions of the
   *  sources. */
  PhotonSourceDistribution *_photon_source_distribution;

  /*! @brief Spectrum for the discrete UV sources. */
  PhotonSourceSpectrum *_photon_source_spectrum;

  /*! @brief Continuous source of UV light. */
  ContinuousPhotonSource *_continuous_photon_source;

  /*! @brief Spectrum for the continuous UV source. */
  PhotonSourceSpectrum *_continuous_photon_source_spectrum;

  /*! @brief Total ionizing luminosity of all sources (in s^-1). */
  double _total_luminosity;

  /*! @brief Object used to compute the combined ionization and temperature
   *  balance at the end of a ray tracing step. */
  TemperatureCalculator *_temperature_calculator;

  /*! @brief Data values for line cooling. */
  const LineCoolingData _line_cooling_data;

  /*! @brief Charge transfer rates. */
  const ChargeTransferRates _charge_transfer_rates;

  /*! @brief Abundance model. */
  const AbundanceModel *_abundance_model;

  /*! @brief Abundances. */
  const Abundances _abundances;

  /*! @brief Cross sections for photoionization. */
  CrossSections *_cross_sections;

  /*! @brief Recombination rates. */
  RecombinationRates *_recombination_rates;

  /*! @brief Reemission handler. */
  DiffuseReemissionHandler *_reemission_handler;

  /*! @brief Log to write logging info to. */
  Log *_log;

  /*! @brief Optional spectrum tracker manager. */
  TrackerManager *_trackers;

  /*! @brief Timer for the total simulation time. */
  Timer _total_timer;

  /*! @brief Timer for serial simulation time. */
  Timer _serial_timer;

  /*! @brief Timer for parallel simulation time. */
  Timer _parallel_timer;

  /*! @brief Timer for the time spent in photon propagations. */
  Timer _photon_propagation_timer;

  /*! @brief Timer for the time spent in cell updates. */
  Timer _cell_update_timer;

  /*! @brief Start time of the program (in CPU cycles). */
  uint_fast64_t _program_start;

  /*! @brief Memory log. */
  MemoryLogger _memory_log;

  /*! @brief Time log. */
  TimeLogger _time_log;

  /*! @brief Output task plot information? */
  const bool _task_plot;

  /*! @brief Output a snapshot before the initial iteration? */
  const bool _output_initial_snapshot;

  const bool _time_dependent_ionization;

  const double _time_dependent_timestep;

public:
  TaskBasedIonizationSimulation(const int_fast32_t num_thread,
                                const std::string parameterfile_name,
                                const bool task_plot = false,
                                const bool output_initial_snapshot = false,
                                Log *log = nullptr);
  ~TaskBasedIonizationSimulation();

  void initialize(DensityFunction *density_function = nullptr);
  void run(DensityGridWriter *density_grid_writer = nullptr);
};

#endif // TASKBASEDIONIZATIONSIMULATION_HPP
