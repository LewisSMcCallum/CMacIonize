/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file GadgetDensityGridWriter.cpp
 *
 * @brief GadgetDensityGridWriter implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "GadgetDensityGridWriter.hpp"
#include "Box.hpp"
#include "CompilerInfo.hpp"
#include "ConfigurationInfo.hpp"
#include "CoordinateVector.hpp"
#include "DensityGrid.hpp"
#include "DensitySubGridCreator.hpp"
#include "DensityValues.hpp"
#include "HDF5Tools.hpp"
#include "HydroDensitySubGrid.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "Utilities.hpp"

#include <vector>

namespace {
struct SphericalOutputDataset {
  hid_t dataset;
  int_fast32_t property;
  int_fast32_t component;
};

inline hid_t create_spherical_dataset(
    const hid_t group, const std::string &name, const hsize_t dimensions[4],
    const hsize_t chunks[4], const bool compression) {
  const hid_t dataspace = H5Screate_simple(4, dimensions, nullptr);
  const hid_t properties = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(properties, 4, chunks);
  if (compression) {
    H5Pset_fletcher32(properties);
    H5Pset_shuffle(properties);
    H5Pset_deflate(properties, 9);
  }
#ifdef HDF5_OLD_API
  const hid_t dataset =
      H5Dcreate(group, name.c_str(), H5T_NATIVE_DOUBLE, dataspace, properties);
#else
  const hid_t dataset =
      H5Dcreate(group, name.c_str(), H5T_NATIVE_DOUBLE, dataspace,
                H5P_DEFAULT, properties, H5P_DEFAULT);
#endif
  H5Pclose(properties);
  H5Sclose(dataspace);
  if (dataset < 0) {
    cmac_error("Failed to create spherical dataset \"%s\".", name.c_str());
  }
  return dataset;
}

inline void write_spherical_block(const hid_t dataset,
                                  const hsize_t offset[4],
                                  const hsize_t count[4],
                                  const std::vector< double > &values) {
  const hid_t filespace = H5Dget_space(dataset);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, nullptr, count,
                      nullptr);
  const hid_t memoryspace = H5Screate_simple(4, count, nullptr);
  if (H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memoryspace, filespace, H5P_DEFAULT,
               values.data()) < 0) {
    cmac_error("Failed to write spherical output block.");
  }
  H5Sclose(memoryspace);
  H5Sclose(filespace);
}

inline void write_spherical_grid(
    const hid_t file, DensitySubGridCreator< DensitySubGrid > &grid_creator,
    const DensityGridWriterFields &fields, const bool compression, Log *log) {
  const CoordinateVector< int_fast32_t > subgrids =
      grid_creator.get_subgrid_layout();
  const CoordinateVector< int_fast32_t > cells =
      grid_creator.get_subgrid_cell_layout();
  const int_fast32_t angular_subgrids = subgrids[1];
  const int_fast32_t angular_cells = angular_subgrids * cells[1];
  const int_fast32_t radial_cells = subgrids[2] * cells[2];
  const hsize_t dimensions[4] = {
      6, static_cast< hsize_t >(angular_cells),
      static_cast< hsize_t >(angular_cells),
      static_cast< hsize_t >(radial_cells)};
  const hsize_t chunks[4] = {
      1, static_cast< hsize_t >(cells[0]),
      static_cast< hsize_t >(cells[1]), static_cast< hsize_t >(cells[2])};

  HDF5Tools::HDF5Group geometry =
      HDF5Tools::create_group(file, "SphericalGrid");
  CoordinateVector<> centre = grid_creator.get_spherical_centre();
  HDF5Tools::write_attribute< CoordinateVector<> >(geometry, "Centre", centre);
  int32_t angular = angular_cells;
  int32_t radial = radial_cells;
  HDF5Tools::write_attribute< int32_t >(geometry, "AngularCellsPerFace",
                                        angular);
  HDF5Tools::write_attribute< int32_t >(geometry, "RadialCells", radial);
  std::string ordering = "face,u,v,r";
  HDF5Tools::write_attribute< std::string >(geometry, "ArrayOrdering",
                                            ordering);
  std::string faces = "+x,-x,+y,-y,+z,-z";
  HDF5Tools::write_attribute< std::string >(geometry, "FaceOrdering", faces);
  std::vector< double > radial_edges =
      grid_creator.get_spherical_radial_edges();
  HDF5Tools::write_dataset< double >(geometry, "RadialEdges", radial_edges);
  std::vector< double > angular_coordinates(angular_cells);
  for (int_fast32_t i = 0; i < angular_cells; ++i) {
    angular_coordinates[i] =
        -1. + 2. * (static_cast< double >(i) + 0.5) / angular_cells;
  }
  HDF5Tools::write_dataset< double >(geometry, "UCoordinates",
                                     angular_coordinates);
  HDF5Tools::write_dataset< double >(geometry, "VCoordinates",
                                     angular_coordinates);
  HDF5Tools::close_group(geometry);

  HDF5Tools::HDF5Group output =
      HDF5Tools::create_group(file, "PartType0");
  std::vector< SphericalOutputDataset > datasets;
  for (int_fast32_t property = 0; property < DENSITYGRIDFIELD_NUMBER;
       ++property) {
    if (!fields.field_present(property)) {
      continue;
    }
    const std::string name = DensityGridWriterFields::get_name(property);
    if (DensityGridWriterFields::get_type(property) ==
        DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE) {
      if (property != DENSITYGRIDFIELD_COORDINATES) {
        cmac_warning("Skipping vector field \"%s\" in spherical output.",
                     name.c_str());
      }
    } else if (DensityGridWriterFields::is_ion_property(property)) {
      for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
        if (fields.ion_present(property, ion)) {
          datasets.push_back(
              {create_spherical_dataset(output, name + get_ion_name(ion),
                                        dimensions, chunks, compression),
               property, ion});
        }
      }
    } else if (DensityGridWriterFields::is_heating_property(property)) {
      for (int_fast32_t heating = 0; heating < NUMBER_OF_HEATINGTERMS;
           ++heating) {
        if (fields.heatingterm_present(property, heating)) {
          datasets.push_back(
              {create_spherical_dataset(output, name + get_ion_name(heating),
                                        dimensions, chunks, compression),
               property, heating});
        }
      }
    } else {
      datasets.push_back(
          {create_spherical_dataset(output, name, dimensions, chunks,
                                    compression),
           property, -1});
    }
  }

  for (auto grid = grid_creator.begin();
       grid != grid_creator.original_end(); ++grid) {
    const CoordinateVector< int_fast32_t > position =
        grid_creator.get_grid_position(grid.get_index());
    const hsize_t offset[4] = {
        static_cast< hsize_t >(position[0] / angular_subgrids),
        static_cast< hsize_t >((position[0] % angular_subgrids) * cells[0]),
        static_cast< hsize_t >(position[1] * cells[1]),
        static_cast< hsize_t >(position[2] * cells[2])};
    const hsize_t count[4] = {
        1, static_cast< hsize_t >(cells[0]),
        static_cast< hsize_t >(cells[1]), static_cast< hsize_t >(cells[2])};
    std::vector< std::vector< double > > values(
        datasets.size(),
        std::vector< double >((*grid).get_number_of_cells()));
    size_t cell_index = 0;
    for (auto cell = (*grid).begin(); cell != (*grid).end(); ++cell) {
      for (size_t field = 0; field < datasets.size(); ++field) {
        const SphericalOutputDataset &description = datasets[field];
        if (DensityGridWriterFields::is_ion_property(description.property)) {
          values[field][cell_index] =
              DensityGridWriterFields::get_scalar_double_ion_value(
                  description.property, description.component, cell);
        } else if (DensityGridWriterFields::is_heating_property(
                       description.property)) {
          values[field][cell_index] =
              DensityGridWriterFields::get_scalar_double_heating_value(
                  description.property, description.component, cell);
        } else {
          values[field][cell_index] =
              DensityGridWriterFields::get_scalar_double_value(
                  description.property, cell);
        }
      }
      ++cell_index;
    }
    for (size_t field = 0; field < datasets.size(); ++field) {
      write_spherical_block(datasets[field].dataset, offset, count,
                            values[field]);
    }
  }
  for (size_t field = 0; field < datasets.size(); ++field) {
    H5Dclose(datasets[field].dataset);
  }
  HDF5Tools::close_group(output);
  if (log) {
    log->write_status("Wrote native spherical arrays with dimensions [6, ",
                      angular_cells, ", ", angular_cells, ", ", radial_cells,
                      "].");
  }
}
} // namespace

/**
 * @brief Constructor.
 *
 * @param prefix Prefix for the name of the file to write.
 * @param output_folder Name of the folder where output files should be placed.
 * @param hydro Flag specifying whether or not hydro is active.
 * @param fields DensityGridWriterFields containing information about which
 * output fields are active.
 * @param log Log to write logging information to.
 * @param padding Number of digits used for the counter in the filenames.
 * @param compression Compress the HDF5 output?
 */
GadgetDensityGridWriter::GadgetDensityGridWriter(
    std::string prefix, std::string output_folder, const bool hydro,
    const DensityGridWriterFields fields, Log *log, uint_fast8_t padding,
    const bool compression)
    : DensityGridWriter(output_folder, hydro, fields, log), _prefix(prefix),
      _padding(padding), _compression(compression) {

  // turn off default HDF5 error handling: we catch errors ourselves
  HDF5Tools::initialize();
  if (_log) {
    _log->write_status("Set up GadgetDensityGridWriter with prefix \"", _prefix,
                       "\".");
    if (_compression) {
      _log->write_status("Compression enabled.");
    } else {
      _log->write_status("Compression disabled.");
    }
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - prefix: Prefix to prepend to all snapshot file names (default: snapshot)
 *  - padding: Number of digits to use in the output file names (default: 3)
 *  - compression: Compress HDF5 datasets? (default: false)
 *
 * @param output_folder Name of the folder where output files should be placed.
 * @param params ParameterFile to read.
 * @param hydro Flag specifying whether or not hydro is active.
 * @param log Log to write logging information to.
 */
GadgetDensityGridWriter::GadgetDensityGridWriter(std::string output_folder,
                                                 ParameterFile &params,
                                                 const bool hydro, Log *log)
    : GadgetDensityGridWriter(
          params.get_value< std::string >("DensityGridWriter:prefix",
                                          "snapshot"),
          output_folder, hydro, DensityGridWriterFields(params, hydro), log,
          params.get_value< uint_fast8_t >("DensityGridWriter:padding", 3),
          params.get_value< bool >("DensityGridWriter:compression", false)) {}

/**
 * @brief Write the file.
 *
 * @param grid DensityGrid to write out.
 * @param iteration Value of the counter to append to the filename.
 * @param params ParameterFile containing the run parameters that should be
 * written to the file.
 * @param time Simulation time (in s).
 * @param hydro_units Internal unit system for the hydro.
 */
void GadgetDensityGridWriter::write(DensityGrid &grid, uint_fast32_t iteration,
                                    ParameterFile &params, double time,
                                    const InternalHydroUnits *hydro_units) {

  std::string filename = Utilities::compose_filename(
      _output_folder, _prefix, "hdf5", iteration, _padding);

  if (_log) {
    _log->write_status("Writing file \"", filename, "\".");
  }

  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_WRITE);

  // write header
  HDF5Tools::HDF5Group group = HDF5Tools::create_group(file, "Header");
  Box<> box = grid.get_box();
  CoordinateVector<> boxsize = box.get_sides();
  HDF5Tools::write_attribute< CoordinateVector<> >(group, "BoxSize", boxsize);
  int32_t dimension = 3;
  HDF5Tools::write_attribute< int32_t >(group, "Dimension", dimension);
  std::vector< uint32_t > flag_entropy(6, 0);
  HDF5Tools::write_attribute< std::vector< uint32_t > >(
      group, "Flag_Entropy_ICs", flag_entropy);
  std::vector< double > masstable(6, 0.);
  HDF5Tools::write_attribute< std::vector< double > >(group, "MassTable",
                                                      masstable);
  int32_t numfiles = 1;
  HDF5Tools::write_attribute< int32_t >(group, "NumFilesPerSnapshot", numfiles);
  std::vector< uint32_t > numpart(6, 0);
  numpart[0] = grid.get_number_of_cells();
  std::vector< uint32_t > numpart_high(6, 0);
  HDF5Tools::write_attribute< std::vector< uint32_t > >(
      group, "NumPart_ThisFile", numpart);
  HDF5Tools::write_attribute< std::vector< uint32_t > >(group, "NumPart_Total",
                                                        numpart);
  HDF5Tools::write_attribute< std::vector< uint32_t > >(
      group, "NumPart_Total_HighWord", numpart_high);
  HDF5Tools::write_attribute< double >(group, "Time", time);
  HDF5Tools::close_group(group);

  // write code info
  group = HDF5Tools::create_group(file, "Code");
  for (auto it = CompilerInfo::begin(); it != CompilerInfo::end(); ++it) {
    std::string key = it.get_key();
    std::string value = it.get_value();
    HDF5Tools::write_attribute< std::string >(group, key, value);
  }
  HDF5Tools::close_group(group);

  // write configuration info
  group = HDF5Tools::create_group(file, "Configuration");
  for (auto it = ConfigurationInfo::begin(); it != ConfigurationInfo::end();
       ++it) {
    std::string key = it.get_key();
    std::string value = it.get_value();
    HDF5Tools::write_attribute< std::string >(group, key, value);
  }
  HDF5Tools::close_group(group);

  // write parameters
  group = HDF5Tools::create_group(file, "Parameters");
  for (auto it = params.begin(); it != params.end(); ++it) {
    std::string key = it.get_key();
    std::string value = it.get_value();
    HDF5Tools::write_attribute< std::string >(group, key, value);
  }
  HDF5Tools::close_group(group);

  // write runtime parameters
  group = HDF5Tools::create_group(file, "RuntimePars");
  std::string timestamp = Utilities::get_timestamp();
  HDF5Tools::write_attribute< std::string >(group, "Creation time", timestamp);
  // an uint_fast32_t does not necessarily have the expected 32-bit size, while
  // we really need a 32-bit variable to write to the file
  uint32_t uint32_iteration = iteration;
  HDF5Tools::write_attribute< uint32_t >(group, "Iteration", uint32_iteration);
  HDF5Tools::close_group(group);

  // write units, we use SI units everywhere
  group = HDF5Tools::create_group(file, "Units");
  double unit_current_in_cgs = 1.;
  double unit_length_in_cgs = 100.;
  double unit_mass_in_cgs = 1000.;
  double unit_temperature_in_cgs = 1.;
  double unit_time_in_cgs = 1.;
  HDF5Tools::write_attribute< double >(group, "Unit current in cgs (U_I)",
                                       unit_current_in_cgs);
  HDF5Tools::write_attribute< double >(group, "Unit length in cgs (U_L)",
                                       unit_length_in_cgs);
  HDF5Tools::write_attribute< double >(group, "Unit mass in cgs (U_M)",
                                       unit_mass_in_cgs);
  HDF5Tools::write_attribute< double >(group, "Unit temperature in cgs (U_T)",
                                       unit_temperature_in_cgs);
  HDF5Tools::write_attribute< double >(group, "Unit time in cgs (U_t)",
                                       unit_time_in_cgs);
  HDF5Tools::close_group(group);

  // write particles
  // to limit memory usage, we first create all datasets, and then add the data
  // in small blocks
  group = HDF5Tools::create_group(file, "PartType0");
  for (int_fast32_t property = 0; property < DENSITYGRIDFIELD_NUMBER;
       ++property) {
    if (_fields.field_present(property)) {
      const std::string name = DensityGridWriterFields::get_name(property);
      if (DensityGridWriterFields::get_type(property) ==
          DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE) {
        HDF5Tools::create_dataset< CoordinateVector<> >(group, name, numpart[0],
                                                        _compression);
      } else {
        if (DensityGridWriterFields::is_ion_property(property)) {
          for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
            if (_fields.ion_present(property, ion)) {
              const std::string prop_name = name + get_ion_name(ion);
              HDF5Tools::create_dataset< double >(group, prop_name, numpart[0],
                                                  _compression);
            }
          }
        } else if (DensityGridWriterFields::is_heating_property(property)) {
          for (int_fast32_t heating = 0; heating < NUMBER_OF_HEATINGTERMS;
               ++heating) {
            if (_fields.heatingterm_present(property, heating)) {
              const std::string prop_name = name + get_ion_name(heating);
              HDF5Tools::create_dataset< double >(group, prop_name, numpart[0],
                                                  _compression);
            }
          }
        } else {
          HDF5Tools::create_dataset< double >(group, name, numpart[0],
                                              _compression);
        }
      }
    }
  }

  const uint_fast32_t blocksize = 10000;
  const uint_fast32_t numblock =
      numpart[0] / blocksize + (numpart[0] % blocksize > 0);
  for (uint_fast32_t iblock = 0; iblock < numblock; ++iblock) {
    const uint_fast32_t offset = iblock * blocksize;
    const uint_fast32_t upper_limit =
        std::min(offset + blocksize, uint_fast32_t(numpart[0]));
    const uint_fast32_t thisblocksize = upper_limit - offset;

    std::vector< std::vector< CoordinateVector<> > > vector_props(
        _fields.get_field_count(DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE),
        std::vector< CoordinateVector<> >(thisblocksize));
    std::vector< std::vector< double > > scalar_props(
        _fields.get_field_count(DENSITYGRIDFIELDTYPE_SCALAR_DOUBLE),
        std::vector< double >(thisblocksize));

    size_t index = 0;
    for (auto it = grid.begin() + offset; it != grid.begin() + upper_limit;
         ++it) {
      uint_fast8_t vector_index = 0;
      uint_fast8_t scalar_index = 0;
      for (int_fast32_t property = 0; property < DENSITYGRIDFIELD_NUMBER;
           ++property) {
        if (_fields.field_present(property)) {
          if (DensityGridWriterFields::get_type(property) ==
              DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE) {
            vector_props[vector_index][index] =
                DensityGridWriterFields::get_vector_double_value(
                    property, it, box.get_anchor(), hydro_units);
            ++vector_index;
          } else {
            if (DensityGridWriterFields::is_ion_property(property)) {
              for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
                if (_fields.ion_present(property, ion)) {
                  scalar_props[scalar_index][index] =
                      DensityGridWriterFields::get_scalar_double_ion_value(
                          property, ion, it);
                  ++scalar_index;
                }
              }
            } else if (DensityGridWriterFields::is_heating_property(property)) {
              for (int_fast32_t heating = 0; heating < NUMBER_OF_HEATINGTERMS;
                   ++heating) {
                if (_fields.heatingterm_present(property, heating)) {
                  scalar_props[scalar_index][index] =
                      DensityGridWriterFields::get_scalar_double_heating_value(
                          property, heating, it);
                  ++scalar_index;
                }
              }
            } else {
              scalar_props[scalar_index][index] =
                  DensityGridWriterFields::get_scalar_double_value(property, it,
                                                                   hydro_units);
              ++scalar_index;
            }
          }
        }
      }
      ++index;
    }

    uint_fast8_t vector_index = 0;
    uint_fast8_t scalar_index = 0;
    for (int_fast32_t property = 0; property < DENSITYGRIDFIELD_NUMBER;
         ++property) {
      if (_fields.field_present(property)) {
        const std::string name = DensityGridWriterFields::get_name(property);
        if (DensityGridWriterFields::get_type(property) ==
            DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE) {
          HDF5Tools::append_dataset< CoordinateVector<> >(
              group, name, offset, vector_props[vector_index]);
          ++vector_index;
        } else {
          if (DensityGridWriterFields::is_ion_property(property)) {
            for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
              if (_fields.ion_present(property, ion)) {
                const std::string prop_name = name + get_ion_name(ion);
                HDF5Tools::append_dataset< double >(group, prop_name, offset,
                                                    scalar_props[scalar_index]);
                ++scalar_index;
              }
            }
          } else if (DensityGridWriterFields::is_heating_property(property)) {
            for (int_fast32_t heating = 0; heating < NUMBER_OF_HEATINGTERMS;
                 ++heating) {
              if (_fields.heatingterm_present(property, heating)) {
                const std::string prop_name = name + get_ion_name(heating);
                HDF5Tools::append_dataset< double >(group, prop_name, offset,
                                                    scalar_props[scalar_index]);
                ++scalar_index;
              }
            }
          } else {
            HDF5Tools::append_dataset< double >(group, name, offset,
                                                scalar_props[scalar_index]);
            ++scalar_index;
          }
        }
      }
    }
  }
  HDF5Tools::close_group(group);

  // close file
  HDF5Tools::close_file(file);
}

/**
 * @brief Write a snapshot for a split grid.
 *
 * @param grid_creator Grid.
 * @param counter Counter value to add to the snapshot file name.
 * @param params ParameterFile containing the run parameters that should be
 * written to the file.
 * @param time Simulation time (in s).
 */
void GadgetDensityGridWriter::write(
    DensitySubGridCreator< DensitySubGrid > &grid_creator,
    const uint_fast32_t counter, ParameterFile &params, double time) {

  std::string filename = Utilities::compose_filename(_output_folder, _prefix,
                                                     "hdf5", counter, _padding);

  if (_log) {
    _log->write_status("Writing file \"", filename, "\".");
  }

  const Box<> box = grid_creator.get_box();

  // we force output of the required fields for now
  //  uint_fast32_t field_flags[DENSITYGRIDFIELD_NUMBER];
  //  for (uint_fast32_t i = 0; i < DENSITYGRIDFIELD_NUMBER; ++i) {
  //    field_flags[i] = 0;
  //  }
  //  field_flags[DENSITYGRIDFIELD_COORDINATES] = 1;
  //  field_flags[DENSITYGRIDFIELD_NUMBER_DENSITY] = 1;
  //  field_flags[DENSITYGRIDFIELD_NEUTRAL_FRACTION] = 1;
  //  const DensityGridWriterFields fields(field_flags);
  // this line is what we actually want...
  const DensityGridWriterFields &fields = _fields;

  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_WRITE);

  // write header
  HDF5Tools::HDF5Group group = HDF5Tools::create_group(file, "Header");
  CoordinateVector<> boxsize = box.get_sides();
  HDF5Tools::write_attribute< CoordinateVector<> >(group, "BoxSize", boxsize);
  int32_t dimension = 3;
  HDF5Tools::write_attribute< int32_t >(group, "Dimension", dimension);
  std::vector< uint32_t > flag_entropy(6, 0);
  HDF5Tools::write_attribute< std::vector< uint32_t > >(
      group, "Flag_Entropy_ICs", flag_entropy);
  std::vector< double > masstable(6, 0.);
  HDF5Tools::write_attribute< std::vector< double > >(group, "MassTable",
                                                      masstable);
  int32_t numfiles = 1;
  HDF5Tools::write_attribute< int32_t >(group, "NumFilesPerSnapshot", numfiles);
  const uint64_t number_of_cells = grid_creator.number_of_cells();
  std::vector< uint32_t > numpart(6, 0);
  numpart[0] = static_cast< uint32_t >(number_of_cells);
  std::vector< uint32_t > numpart_high(6, 0);
  numpart_high[0] = static_cast< uint32_t >(number_of_cells >> 32);
  HDF5Tools::write_attribute< std::vector< uint32_t > >(
      group, "NumPart_ThisFile", numpart);
  HDF5Tools::write_attribute< std::vector< uint32_t > >(group, "NumPart_Total",
                                                        numpart);
  HDF5Tools::write_attribute< std::vector< uint32_t > >(
      group, "NumPart_Total_HighWord", numpart_high);
  HDF5Tools::write_attribute< double >(group, "Time", time);
  HDF5Tools::close_group(group);

  // write code info
  group = HDF5Tools::create_group(file, "Code");
  for (auto it = CompilerInfo::begin(); it != CompilerInfo::end(); ++it) {
    std::string key = it.get_key();
    std::string value = it.get_value();
    HDF5Tools::write_attribute< std::string >(group, key, value);
  }
  HDF5Tools::close_group(group);

  // write configuration info
  group = HDF5Tools::create_group(file, "Configuration");
  for (auto it = ConfigurationInfo::begin(); it != ConfigurationInfo::end();
       ++it) {
    std::string key = it.get_key();
    std::string value = it.get_value();
    HDF5Tools::write_attribute< std::string >(group, key, value);
  }
  HDF5Tools::close_group(group);

  // write parameters
  group = HDF5Tools::create_group(file, "Parameters");
  for (auto it = params.begin(); it != params.end(); ++it) {
    std::string key = it.get_key();
    std::string value = it.get_value();
    HDF5Tools::write_attribute< std::string >(group, key, value);
  }
  HDF5Tools::close_group(group);

  // write runtime parameters
  group = HDF5Tools::create_group(file, "RuntimePars");
  std::string timestamp = Utilities::get_timestamp();
  HDF5Tools::write_attribute< std::string >(group, "Creation time", timestamp);
  // an uint_fast32_t does not necessarily have the expected 32-bit size, while
  // we really need a 32-bit variable to write to the file
  uint32_t uint32_iteration = counter;
  HDF5Tools::write_attribute< uint32_t >(group, "Iteration", uint32_iteration);
  HDF5Tools::close_group(group);

  // write units, we use SI units everywhere
  group = HDF5Tools::create_group(file, "Units");
  double unit_current_in_cgs = 1.;
  double unit_length_in_cgs = 100.;
  double unit_mass_in_cgs = 1000.;
  double unit_temperature_in_cgs = 1.;
  double unit_time_in_cgs = 1.;
  HDF5Tools::write_attribute< double >(group, "Unit current in cgs (U_I)",
                                       unit_current_in_cgs);
  HDF5Tools::write_attribute< double >(group, "Unit length in cgs (U_L)",
                                       unit_length_in_cgs);
  HDF5Tools::write_attribute< double >(group, "Unit mass in cgs (U_M)",
                                       unit_mass_in_cgs);
  HDF5Tools::write_attribute< double >(group, "Unit temperature in cgs (U_T)",
                                       unit_temperature_in_cgs);
  HDF5Tools::write_attribute< double >(group, "Unit time in cgs (U_t)",
                                       unit_time_in_cgs);
  HDF5Tools::close_group(group);

  if (grid_creator.is_spherical()) {
    write_spherical_grid(file, grid_creator, fields, _compression, _log);
    HDF5Tools::close_file(file);
    return;
  }

  // write particles
  // to limit memory usage, we first create all datasets, and then add the data
  // in small blocks
  group = HDF5Tools::create_group(file, "PartType0");
  for (int_fast32_t property = 0; property < DENSITYGRIDFIELD_NUMBER;
       ++property) {
    if (fields.field_present(property)) {
      const std::string name = DensityGridWriterFields::get_name(property);
      if (DensityGridWriterFields::get_type(property) ==
          DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE) {
        HDF5Tools::create_dataset< CoordinateVector<> >(group, name, numpart[0],
                                                        _compression);
      } else {
        if (DensityGridWriterFields::is_ion_property(property)) {
          for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
            if (fields.ion_present(property, ion)) {
              const std::string prop_name = name + get_ion_name(ion);
              HDF5Tools::create_dataset< double >(group, prop_name, numpart[0],
                                                  _compression);
            }
          }
        } else if (DensityGridWriterFields::is_heating_property(property)) {
          for (int_fast32_t heating = 0; heating < NUMBER_OF_HEATINGTERMS;
               ++heating) {
            if (fields.heatingterm_present(property, heating)) {
              const std::string prop_name = name + get_ion_name(heating);
              HDF5Tools::create_dataset< double >(group, prop_name, numpart[0],
                                                  _compression);
            }
          }
        } else {
          HDF5Tools::create_dataset< double >(group, name, numpart[0],
                                              _compression);
        }
      }
    }
  }

  const uint_fast32_t blocksize = 10000;
  uint_fast32_t block_offset = 0;
  for (auto gridit = grid_creator.begin();
       gridit != grid_creator.original_end(); ++gridit) {

    const uint_fast32_t numblock =
        (*gridit).get_number_of_cells() / blocksize +
        ((*gridit).get_number_of_cells() % blocksize > 0);
    for (uint_fast32_t iblock = 0; iblock < numblock; ++iblock) {
      const uint_fast32_t offset = iblock * blocksize;
      const uint_fast32_t upper_limit = std::min(
          offset + blocksize, uint_fast32_t((*gridit).get_number_of_cells()));
      const uint_fast32_t thisblocksize = upper_limit - offset;

      std::vector< std::vector< CoordinateVector<> > > vector_props(
          fields.get_field_count(DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE),
          std::vector< CoordinateVector<> >(thisblocksize));
      std::vector< std::vector< double > > scalar_props(
          fields.get_field_count(DENSITYGRIDFIELDTYPE_SCALAR_DOUBLE),
          std::vector< double >(thisblocksize));

      size_t index = 0;
      for (auto cellit = (*gridit).begin() + offset;
           cellit != (*gridit).begin() + upper_limit; ++cellit) {
        uint_fast8_t vector_index = 0;
        uint_fast8_t scalar_index = 0;
        for (int_fast32_t property = 0; property < DENSITYGRIDFIELD_NUMBER;
             ++property) {
          if (fields.field_present(property)) {
            if (DensityGridWriterFields::get_type(property) ==
                DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE) {
              vector_props[vector_index][index] =
                  DensityGridWriterFields::get_vector_double_value(
                      property, cellit, box.get_anchor());
              ++vector_index;
            } else {
              if (DensityGridWriterFields::is_ion_property(property)) {
                for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
                  if (fields.ion_present(property, ion)) {
                    scalar_props[scalar_index][index] =
                        DensityGridWriterFields::get_scalar_double_ion_value(
                            property, ion, cellit);
                    ++scalar_index;
                  }
                }
              } else if (DensityGridWriterFields::is_heating_property(
                             property)) {
                for (int_fast32_t heating = 0; heating < NUMBER_OF_HEATINGTERMS;
                     ++heating) {
                  if (fields.heatingterm_present(property, heating)) {
                    scalar_props[scalar_index][index] =
                        DensityGridWriterFields::
                            get_scalar_double_heating_value(property, heating,
                                                            cellit);
                    ++scalar_index;
                  }
                }
              } else {
                scalar_props[scalar_index][index] =
                    DensityGridWriterFields::get_scalar_double_value(property,
                                                                     cellit);
                ++scalar_index;
              }
            }
          }
        }
        ++index;
      }

      uint_fast8_t vector_index = 0;
      uint_fast8_t scalar_index = 0;
      for (int_fast32_t property = 0; property < DENSITYGRIDFIELD_NUMBER;
           ++property) {
        if (fields.field_present(property)) {
          const std::string name = DensityGridWriterFields::get_name(property);
          if (DensityGridWriterFields::get_type(property) ==
              DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE) {
            HDF5Tools::append_dataset< CoordinateVector<> >(
                group, name, block_offset + offset, vector_props[vector_index]);
            ++vector_index;
          } else {
            if (DensityGridWriterFields::is_ion_property(property)) {
              for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
                if (fields.ion_present(property, ion)) {
                  const std::string prop_name = name + get_ion_name(ion);
                  HDF5Tools::append_dataset< double >(
                      group, prop_name, block_offset + offset,
                      scalar_props[scalar_index]);
                  ++scalar_index;
                }
              }
            } else if (DensityGridWriterFields::is_heating_property(property)) {
              for (int_fast32_t heating = 0; heating < NUMBER_OF_HEATINGTERMS;
                   ++heating) {
                if (fields.heatingterm_present(property, heating)) {
                  const std::string prop_name = name + get_ion_name(heating);
                  HDF5Tools::append_dataset< double >(
                      group, prop_name, block_offset + offset,
                      scalar_props[scalar_index]);
                  ++scalar_index;
                }
              }
            } else {
              HDF5Tools::append_dataset< double >(group, name,
                                                  block_offset + offset,
                                                  scalar_props[scalar_index]);
              ++scalar_index;
            }
          }
        }
      }
    }
    block_offset += (*gridit).get_number_of_cells();
  }
  HDF5Tools::close_group(group);

  // close file
  HDF5Tools::close_file(file);
}

/**
 * @brief Write a snapshot for a split grid with hydro.
 *
 * @param grid_creator Grid.
 * @param counter Counter value to add to the snapshot file name.
 * @param params ParameterFile containing the run parameters that should be
 * written to the file.
 * @param time Simulation time (in s).
 */
void GadgetDensityGridWriter::write(
    DensitySubGridCreator< HydroDensitySubGrid > &grid_creator,
    const uint_fast32_t counter, ParameterFile &params, double time) {

  std::string filename = Utilities::compose_filename(_output_folder, _prefix,
                                                     "hdf5", counter, _padding);

  if (_log) {
    _log->write_status("Writing file \"", filename, "\".");
  }

  const Box<> box = grid_creator.get_box();

  // we force output of the required fields for now
  //  uint_fast32_t field_flags[DENSITYGRIDFIELD_NUMBER];
  //  for (uint_fast32_t i = 0; i < DENSITYGRIDFIELD_NUMBER; ++i) {
  //    field_flags[i] = 0;
  //  }
  //  field_flags[DENSITYGRIDFIELD_COORDINATES] = 1;
  //  field_flags[DENSITYGRIDFIELD_NUMBER_DENSITY] = 1;
  //  field_flags[DENSITYGRIDFIELD_NEUTRAL_FRACTION] = 1;
  //  const DensityGridWriterFields fields(field_flags);
  // this line is what we actually want...
  const DensityGridWriterFields &fields = _fields;

  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_WRITE);

  // write header
  HDF5Tools::HDF5Group group = HDF5Tools::create_group(file, "Header");
  CoordinateVector<> boxsize = box.get_sides();
  HDF5Tools::write_attribute< CoordinateVector<> >(group, "BoxSize", boxsize);
  int32_t dimension = 3;
  HDF5Tools::write_attribute< int32_t >(group, "Dimension", dimension);
  std::vector< uint32_t > flag_entropy(6, 0);
  HDF5Tools::write_attribute< std::vector< uint32_t > >(
      group, "Flag_Entropy_ICs", flag_entropy);
  std::vector< double > masstable(6, 0.);
  HDF5Tools::write_attribute< std::vector< double > >(group, "MassTable",
                                                      masstable);
  int32_t numfiles = 1;
  HDF5Tools::write_attribute< int32_t >(group, "NumFilesPerSnapshot", numfiles);
  const uint64_t number_of_cells = grid_creator.number_of_cells();
  std::vector< uint32_t > numpart(6, 0);
  numpart[0] = static_cast< uint32_t >(number_of_cells);
  std::vector< uint32_t > numpart_high(6, 0);
  numpart_high[0] = static_cast< uint32_t >(number_of_cells >> 32);
  HDF5Tools::write_attribute< std::vector< uint32_t > >(
      group, "NumPart_ThisFile", numpart);
  HDF5Tools::write_attribute< std::vector< uint32_t > >(group, "NumPart_Total",
                                                        numpart);
  HDF5Tools::write_attribute< std::vector< uint32_t > >(
      group, "NumPart_Total_HighWord", numpart_high);
  HDF5Tools::write_attribute< double >(group, "Time", time);
  HDF5Tools::close_group(group);

  // write code info
  group = HDF5Tools::create_group(file, "Code");
  for (auto it = CompilerInfo::begin(); it != CompilerInfo::end(); ++it) {
    std::string key = it.get_key();
    std::string value = it.get_value();
    HDF5Tools::write_attribute< std::string >(group, key, value);
  }
  HDF5Tools::close_group(group);

  // write configuration info
  group = HDF5Tools::create_group(file, "Configuration");
  for (auto it = ConfigurationInfo::begin(); it != ConfigurationInfo::end();
       ++it) {
    std::string key = it.get_key();
    std::string value = it.get_value();
    HDF5Tools::write_attribute< std::string >(group, key, value);
  }
  HDF5Tools::close_group(group);

  // write parameters
  group = HDF5Tools::create_group(file, "Parameters");
  for (auto it = params.begin(); it != params.end(); ++it) {
    std::string key = it.get_key();
    std::string value = it.get_value();
    HDF5Tools::write_attribute< std::string >(group, key, value);
  }
  HDF5Tools::close_group(group);

  // write runtime parameters
  group = HDF5Tools::create_group(file, "RuntimePars");
  std::string timestamp = Utilities::get_timestamp();
  HDF5Tools::write_attribute< std::string >(group, "Creation time", timestamp);
  // an uint_fast32_t does not necessarily have the expected 32-bit size, while
  // we really need a 32-bit variable to write to the file
  uint32_t uint32_iteration = counter;
  HDF5Tools::write_attribute< uint32_t >(group, "Iteration", uint32_iteration);
  HDF5Tools::close_group(group);

  // write units, we use SI units everywhere
  group = HDF5Tools::create_group(file, "Units");
  double unit_current_in_cgs = 1.;
  double unit_length_in_cgs = 100.;
  double unit_mass_in_cgs = 1000.;
  double unit_temperature_in_cgs = 1.;
  double unit_time_in_cgs = 1.;
  HDF5Tools::write_attribute< double >(group, "Unit current in cgs (U_I)",
                                       unit_current_in_cgs);
  HDF5Tools::write_attribute< double >(group, "Unit length in cgs (U_L)",
                                       unit_length_in_cgs);
  HDF5Tools::write_attribute< double >(group, "Unit mass in cgs (U_M)",
                                       unit_mass_in_cgs);
  HDF5Tools::write_attribute< double >(group, "Unit temperature in cgs (U_T)",
                                       unit_temperature_in_cgs);
  HDF5Tools::write_attribute< double >(group, "Unit time in cgs (U_t)",
                                       unit_time_in_cgs);
  HDF5Tools::close_group(group);

  // write particles
  // to limit memory usage, we first create all datasets, and then add the data
  // in small blocks
  group = HDF5Tools::create_group(file, "PartType0");
  for (int_fast32_t property = 0; property < DENSITYGRIDFIELD_NUMBER;
       ++property) {
    if (fields.field_present(property)) {
      const std::string name = DensityGridWriterFields::get_name(property);
      if (DensityGridWriterFields::get_type(property) ==
          DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE) {
        HDF5Tools::create_dataset< CoordinateVector<> >(group, name, numpart[0],
                                                        _compression);
      } else {
        if (DensityGridWriterFields::is_ion_property(property)) {
          for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
            if (fields.ion_present(property, ion)) {
              const std::string prop_name = name + get_ion_name(ion);
              HDF5Tools::create_dataset< double >(group, prop_name, numpart[0],
                                                  _compression);
            }
          }
        } else if (DensityGridWriterFields::is_heating_property(property)) {
          for (int_fast32_t heating = 0; heating < NUMBER_OF_HEATINGTERMS;
               ++heating) {
            if (fields.heatingterm_present(property, heating)) {
              const std::string prop_name = name + get_ion_name(heating);
              HDF5Tools::create_dataset< double >(group, prop_name, numpart[0],
                                                  _compression);
            }
          }
        } else {
          HDF5Tools::create_dataset< double >(group, name, numpart[0],
                                              _compression);
        }
      }
    }
  }

  const uint_fast32_t blocksize = 10000;
  uint_fast32_t block_offset = 0;
  for (auto gridit = grid_creator.begin();
       gridit != grid_creator.original_end(); ++gridit) {

    const uint_fast32_t numblock =
        (*gridit).get_number_of_cells() / blocksize +
        ((*gridit).get_number_of_cells() % blocksize > 0);
    for (uint_fast32_t iblock = 0; iblock < numblock; ++iblock) {
      const uint_fast32_t offset = iblock * blocksize;
      const uint_fast32_t upper_limit = std::min(
          offset + blocksize, uint_fast32_t((*gridit).get_number_of_cells()));
      const uint_fast32_t thisblocksize = upper_limit - offset;

      std::vector< std::vector< CoordinateVector<> > > vector_props(
          fields.get_field_count(DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE),
          std::vector< CoordinateVector<> >(thisblocksize));
      std::vector< std::vector< double > > scalar_props(
          fields.get_field_count(DENSITYGRIDFIELDTYPE_SCALAR_DOUBLE),
          std::vector< double >(thisblocksize));

      size_t index = 0;
      for (auto cellit = (*gridit).hydro_begin() + offset;
           cellit != (*gridit).hydro_begin() + upper_limit; ++cellit) {
        uint_fast8_t vector_index = 0;
        uint_fast8_t scalar_index = 0;
        for (int_fast32_t property = 0; property < DENSITYGRIDFIELD_NUMBER;
             ++property) {
          if (fields.field_present(property)) {
            if (DensityGridWriterFields::get_type(property) ==
                DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE) {
              vector_props[vector_index][index] =
                  DensityGridWriterFields::get_vector_double_value(
                      property, cellit, box.get_anchor());
              ++vector_index;
            } else {
              if (DensityGridWriterFields::is_ion_property(property)) {
                for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
                  if (fields.ion_present(property, ion)) {
                    scalar_props[scalar_index][index] =
                        DensityGridWriterFields::get_scalar_double_ion_value(
                            property, ion, cellit);
                    ++scalar_index;
                  }
                }
              } else if (DensityGridWriterFields::is_heating_property(
                             property)) {
                for (int_fast32_t heating = 0; heating < NUMBER_OF_HEATINGTERMS;
                     ++heating) {
                  if (fields.heatingterm_present(property, heating)) {
                    scalar_props[scalar_index][index] =
                        DensityGridWriterFields::
                            get_scalar_double_heating_value(property, heating,
                                                            cellit);
                    ++scalar_index;
                  }
                }
              } else {
                scalar_props[scalar_index][index] =
                    DensityGridWriterFields::get_scalar_double_value(property,
                                                                     cellit);
                ++scalar_index;
              }
            }
          }
        }
        ++index;
      }

      uint_fast8_t vector_index = 0;
      uint_fast8_t scalar_index = 0;
      for (int_fast32_t property = 0; property < DENSITYGRIDFIELD_NUMBER;
           ++property) {
        if (fields.field_present(property)) {
          const std::string name = DensityGridWriterFields::get_name(property);
          if (DensityGridWriterFields::get_type(property) ==
              DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE) {
            HDF5Tools::append_dataset< CoordinateVector<> >(
                group, name, block_offset + offset, vector_props[vector_index]);
            ++vector_index;
          } else {
            if (DensityGridWriterFields::is_ion_property(property)) {
              for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
                if (fields.ion_present(property, ion)) {
                  const std::string prop_name = name + get_ion_name(ion);
                  HDF5Tools::append_dataset< double >(
                      group, prop_name, block_offset + offset,
                      scalar_props[scalar_index]);
                  ++scalar_index;
                }
              }
            } else if (DensityGridWriterFields::is_heating_property(property)) {
              for (int_fast32_t heating = 0; heating < NUMBER_OF_HEATINGTERMS;
                   ++heating) {
                if (fields.heatingterm_present(property, heating)) {
                  const std::string prop_name = name + get_ion_name(heating);
                  HDF5Tools::append_dataset< double >(
                      group, prop_name, block_offset + offset,
                      scalar_props[scalar_index]);
                  ++scalar_index;
                }
              }
            } else {
              HDF5Tools::append_dataset< double >(group, name,
                                                  block_offset + offset,
                                                  scalar_props[scalar_index]);
              ++scalar_index;
            }
          }
        }
      }
    }
    block_offset += (*gridit).get_number_of_cells();
  }
  HDF5Tools::close_group(group);

  // close file
  HDF5Tools::close_file(file);
}
std::string GadgetDensityGridWriter::get_snapshot_filename(
    const uint_fast32_t counter) const {
  return Utilities::compose_filename(_output_folder, _prefix, "hdf5", counter,
                                     _padding);
}
