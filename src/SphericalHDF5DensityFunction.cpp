/*******************************************************************************
 * This file is part of CMacIonize
 ******************************************************************************/

/**
 * @file SphericalHDF5DensityFunction.cpp
 *
 * @brief DensityFunction for a cell-centred cubed-sphere HDF5 map.
 */
#include "SphericalHDF5DensityFunction.hpp"

#include "Error.hpp"
#include "HDF5Tools.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

#include <algorithm>
#include <cmath>

namespace {
inline std::vector< float > read_number_density(hid_t group,
                                                const std::string &name,
                                                const int_fast32_t angular,
                                                const int_fast32_t radial) {
  const hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);
  if (dataset < 0) {
    cmac_error("Failed to open dataset \"%s\".", name.c_str());
  }
  const hid_t dataspace = H5Dget_space(dataset);
  const int dimensions = H5Sget_simple_extent_ndims(dataspace);
  if (dimensions != 4) {
    cmac_error("Dataset \"%s\" should have four dimensions.", name.c_str());
  }
  std::vector< hsize_t > shape(dimensions);
  H5Sget_simple_extent_dims(dataspace, shape.data(), nullptr);
  if (shape[0] != 6 || shape[1] != static_cast< hsize_t >(angular) ||
      shape[2] != static_cast< hsize_t >(angular) ||
      shape[3] != static_cast< hsize_t >(radial)) {
    cmac_error("Dataset \"%s\" dimensions do not match the spherical grid.",
               name.c_str());
  }
  const size_t size = shape[0] * shape[1] * shape[2] * shape[3];
  std::vector< float > values(size);
  if (H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
              values.data()) < 0) {
    cmac_error("Failed to read dataset \"%s\".", name.c_str());
  }
  H5Sclose(dataspace);
  H5Dclose(dataset);
  return values;
}

inline int_fast8_t face_for_direction(const CoordinateVector<> &x) {
  int_fast8_t axis = 0;
  if (std::abs(x[1]) > std::abs(x[axis])) {
    axis = 1;
  }
  if (std::abs(x[2]) > std::abs(x[axis])) {
    axis = 2;
  }
  return 2 * axis + (x[axis] < 0.);
}

inline CoordinateVector<> face_normal(const int_fast8_t face) {
  static const double axes[6][3] = {{1, 0, 0},  {-1, 0, 0}, {0, 1, 0},
                                     {0, -1, 0}, {0, 0, 1},  {0, 0, -1}};
  return CoordinateVector<>(axes[face][0], axes[face][1], axes[face][2]);
}

inline CoordinateVector<> face_u(const int_fast8_t face) {
  static const double axes[6][3] = {{0, 1, 0}, {0, -1, 0}, {-1, 0, 0},
                                     {1, 0, 0}, {0, 1, 0},  {0, 1, 0}};
  return CoordinateVector<>(axes[face][0], axes[face][1], axes[face][2]);
}

inline CoordinateVector<> face_v(const int_fast8_t face) {
  static const double axes[6][3] = {{0, 0, 1}, {0, 0, 1}, {0, 0, 1},
                                     {0, 0, 1}, {-1, 0, 0}, {1, 0, 0}};
  return CoordinateVector<>(axes[face][0], axes[face][1], axes[face][2]);
}

inline double dot(const CoordinateVector<> &a, const CoordinateVector<> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
} // namespace

SphericalHDF5DensityFunction::SphericalHDF5DensityFunction(
    ParameterFile &params, Log *log)
    : _centre(params.get_physical_vector< QUANTITY_LENGTH >(
          "SphericalDensityGrid:centre", "[0. m, 0. m, 0. m]")),
      _angular_cells(params.get_value< int_fast32_t >(
          "SphericalDensityGrid:angular cells per face", 32)),
      _radial_cells(params.get_value< int_fast32_t >(
          "SphericalDensityGrid:number of radial cells", 64)),
      _temperature(params.get_physical_value< QUANTITY_TEMPERATURE >(
          "DensityFunction:temperature", "100. K")),
      _dust_gas_ratio(
          params.get_value< double >("DensityFunction:dust to gas", 0.01)),
      _fraction_silicates(params.get_value< double >(
          "DensityFunction:fraction silicates", 1.)) {
  const std::string filename = params.get_filename("DensityFunction:filename");
  HDF5Tools::initialize();
  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_READ);
  HDF5Tools::HDF5Group group =
      HDF5Tools::open_group(file, "/SphericalGrid");
  _number_densities = read_number_density(
      group, "NumberDensity", _angular_cells, _radial_cells);
  _radial_edges = HDF5Tools::read_dataset< double >(group, "RadialEdges");
  HDF5Tools::close_group(group);
  HDF5Tools::close_file(file);

  const size_t expected = static_cast< size_t >(6) * _angular_cells *
                          _angular_cells * _radial_cells;
  if (_number_densities.size() != expected ||
      _radial_edges.size() != static_cast< size_t >(_radial_cells + 1)) {
    cmac_error("Spherical HDF5 dimensions do not match the parameter file.");
  }
  if (log) {
    log->write_status("Read ", expected,
                      " spherical density cells from \"", filename, "\".");
  }
}

DensityValues
SphericalHDF5DensityFunction::operator()(const Cell &cell) {
  const CoordinateVector<> x = cell.get_cell_midpoint() - _centre;
  const double radius = x.norm();
  const int_fast8_t face = face_for_direction(x);
  const double denominator = dot(face_normal(face), x);
  const double u = dot(face_u(face), x) / denominator;
  const double v = dot(face_v(face), x) / denominator;
  const int_fast32_t iu = std::min(
      _angular_cells - 1,
      std::max< int_fast32_t >(
          0, std::floor(0.5 * (u + 1.) * _angular_cells)));
  const int_fast32_t iv = std::min(
      _angular_cells - 1,
      std::max< int_fast32_t >(
          0, std::floor(0.5 * (v + 1.) * _angular_cells)));
  const std::vector< double >::const_iterator edge =
      std::upper_bound(_radial_edges.begin(), _radial_edges.end(), radius);
  const int_fast32_t ir =
      std::min(_radial_cells - 1,
               std::max< int_fast32_t >(
                   0, edge - _radial_edges.begin() - 1));
  const size_t index =
      ((static_cast< size_t >(face) * _angular_cells + iu) *
           _angular_cells +
       iv) *
          _radial_cells +
      ir;

  DensityValues values;
  values.set_number_density(_number_densities[index]);
  values.set_temperature(_temperature);
  values.set_dust_gas_ratio(_dust_gas_ratio);
  values.set_fraction_silicates(_fraction_silicates);
  values.set_ionic_fraction(ION_H_n, 0.999);
  return values;
}
