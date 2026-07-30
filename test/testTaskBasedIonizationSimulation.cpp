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
 * @file testTaskBasedIonizationSimulation.cpp
 *
 * @brief Unit test for the TaskBasedIonizationSimulation library.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Assert.hpp"
#include "Log.hpp"
#include "TaskBasedIonizationSimulation.hpp"

#include <string>
#include <vector>

class CapturingLog : public Log {
private:
  std::vector< std::string > _messages;

protected:
  virtual void write_message(std::string message) {
    _messages.push_back(message);
  }

public:
  CapturingLog() : Log(LOGLEVEL_STATUS) {}

  bool contains(const std::string &text) const {
    for (const std::string &message : _messages) {
      if (message.find(text) != std::string::npos) {
        return true;
      }
    }
    return false;
  }
};

/**
 * @brief Unit test for the TaskBasedIonizationSimulation library.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  CapturingLog log;
  TaskBasedIonizationSimulation simulation(
      1, "test_taskbasedionizationsimulation.param", false, false, &log);
  simulation.initialize(nullptr);
  simulation.run(nullptr);

  assert_condition(log.contains("Starting loop 0 with 10000 photon packets."));
  assert_condition(log.contains("Starting loop 1 with 20000 photon packets."));
  assert_condition(log.contains("Starting loop 8 with 20000 photon packets."));
  assert_condition(log.contains("Starting loop 9 with 30000 photon packets."));

  return 0;
}
