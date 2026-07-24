/*******************************************************************************
 * This file is part of CMacIonize
 ******************************************************************************/

/**
 * @file SphericalHDF5DensityFunction.hpp
 *
 * @brief DensityFunction for a cell-centred cubed-sphere HDF5 map.
 */
#ifndef SPHERICALHDF5DENSITYFUNCTION_HPP
#define SPHERICALHDF5DENSITYFUNCTION_HPP

#include "DensityFunction.hpp"

#include <string>
#include <vector>

class Log;
class ParameterFile;

class SphericalHDF5DensityFunction : public DensityFunction {
private:
  CoordinateVector<> _centre;
  int_fast32_t _angular_cells;
  int_fast32_t _radial_cells;
  std::vector< double > _radial_edges;
  // Keep the prepared map in its native float format: the N=256 map already
  // contains over 200 million cells, so promoting the input would waste 0.8 GB.
  std::vector< float > _number_densities;
  double _temperature;
  double _dust_gas_ratio;
  double _fraction_silicates;

public:
  SphericalHDF5DensityFunction(ParameterFile &params, Log *log = nullptr);

  virtual DensityValues operator()(const Cell &cell);
};

#endif
