
#include "CollisionalRates.hpp"
#include "CollisionalRatesDataLocation.hpp"
#include "ElementNames.hpp"
#include "UnitConverter.hpp"
#include <fstream>
#include <sstream>
#include <string>

/**
 * @brief Constructor.
 */
CollisionalRates::CollisionalRates() {

std::cout << "HAS " << NUMBER_OF_IONNAMES << " different ions." << std::endl;

  _collisional_rates.resize(NUMBER_OF_IONNAMES);
  _temperatures.resize(250);


  for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {

    std::cout << "Setting up data from file - " << get_ion_collisional_filename(i) << std::endl;

    _collisional_rates[i].resize(250);


    std::stringstream filenamestream;
    filenamestream << COLLISIONALDATALOCATION
                   << get_ion_collisional_filename(i);




    std::ifstream drfile(filenamestream.str());

    // skip the first two lines
    std::string line;
    std::getline(drfile, line);
    std::getline(drfile, line);
    // now parse the remaining lines
    for (uint_fast32_t j = 0; j < 250; ++j) {
      std::getline(drfile, line);

      std::istringstream linestream(line);

      //overwriting temperatures 181 times.... not great but fine

      linestream >> _temperatures[j] >> _collisional_rates[i][j];
      // temperature to K
      _temperatures[j] = std::pow(10.,_temperatures[j]);
      // cm^3/s to m^3/s

      _collisional_rates[i][j] *= 1.e-6;
  

      //std::cout << _temperatures[j] << "kelvin with " << _collisional_rates[i][j] << std::endl;
    }

  }

  _min_logT = std::log(_temperatures[0]);
  _inverse_avg_dlogT =
      250. /
      (std::log(_temperatures[250 - 1]) -
       _min_logT);

}
