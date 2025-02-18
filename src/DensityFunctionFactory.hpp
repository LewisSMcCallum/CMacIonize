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
 * @file DensityFunctionFactory.hpp
 *
 * @brief Factory for DensityFunction implementations.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYFUNCTIONFACTORY_HPP
#define DENSITYFUNCTIONFACTORY_HPP

#include "Configuration.hpp"
#include "DensityFunction.hpp"
#include "Error.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

// non library dependent implementations
#include "AsciiFileDensityFunction.hpp"
#include "BlockSyntaxDensityFunction.hpp"
#include "BondiProfileDensityFunction.hpp"
#include "CoredDMProfileDensityFunction.hpp"
#include "DiscICDensityFunction.hpp"
#include "DiscPatchDensityFunction.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "InterpolatedDensityFunction.hpp"
#include "PhantomSnapshotDensityFunction.hpp"
#include "SPHNGSnapshotDensityFunction.hpp"
#include "SpiralGalaxyDensityFunction.hpp"

// HDF5 dependent implementations
#ifdef HAVE_HDF5
#include "AmunSnapshotDensityFunction.hpp"
#include "ArepoSnapshotDensityFunction.hpp"
#include "BufferedCMacIonizeSnapshotDensityFunction.hpp"
#include "CMacIonizeSnapshotDensityFunction.hpp"
#include "FLASHSnapshotDensityFunction.hpp"
#include "GadgetSnapshotDensityFunction.hpp"
#include "HDF5DensityFunction.hpp"
#endif

#include <string>

/**
 * @brief Factory for DensityFunction implementations.
 */
class DensityFunctionFactory {
public:
  /**
   * @brief Method that checks if the requested DensityFunction implementation
   * requires HDF5.
   *
   * @param type Requested DensityFunction type.
   * @param log Log to write logging info to.
   */
  static void check_hdf5(std::string type, Log *log = nullptr) {
    if (type == "AmunSnapshot" || type == "BufferedCMacIonizeSnapshot" ||
        type == "CMacIonizeSnapshot" || type == "FLASHSnapshot" ||
        type == "GadgetSnapshot" || type == "HDF5" || type == "ArepoSnapshot") {
      if (log) {
        log->write_error("Cannot create an instance of ", type,
                         "DensityFunction, since the code was "
                         "compiled without HDF5 support.");
      }
      cmac_error("A %sDensityFunction requires HDF5. However, the code "
                 "was compiled without HDF5 support!",
                 type.c_str());
    }
  }

  /**
   * @brief Generate a DensityFunction based on the type chosen in the parameter
   * file.
   *
   * Supported types are (default: Homogeneous):
   *  - AmunSnapshot: Implementation that reads a density grid from an Amun
   *    snapshot file
   *  - AsciiFile: Implementation that reads a density grid from an ASCII text
   *    file
   *  - BlockSyntax: Implementation that reads a geometrically constructed
   *    density field from a text file containing block syntax
   *  - BondiProfile: Bondi accretion profile.
   *  - CoredDMProfile: Density profile in hydrostatic equilibrium with a cored
   *    DM density profile.
   *  - DiscIC: Constant value density field with a tangential velocity profile
   *  - DiscPatch: Hydrostatic disc patch (Creasey, Theuns & Bower, 2013).
   *  - Homogeneous: Constant value density field.
   *  - Interpolated: Implementation that reads a density field from a text file
   *    and interpolates on it
   *  - SPHNGSnapshot: Implementation that reads a density field from a snapshot
   *    file of the SPH code SPHNG
   *  - SpiralGalaxy: Implementation that sets up a diffuse galactic density
   *    field
   *  - CMacIonizeSnapshot: Implementation that reads a density field from a
   *    snapshot of another CMacIonize run
   *  - FLASHSnapshot: Implementation that reads a density field from the AMR
   *    simulation code FLASH
   *  - GadgetSnapshot: Implementation that reads a density field from the HDF5
   *    file format of the SPH simulation code Gadget2 (also supported by SWIFT,
   *    AREPO, GIZMO and Shadowfax; the CMacIonize snapshot format is a variant
   *    of this format)
   *  - PhantomSnapshot: Implementation that reads a density field from the
   *    binary file format used by the SPH simulation code Phantom
   *
   * @param params ParameterFile containing the parameters used by the specific
   * implementation.
   * @param log Log to write logging information to.
   * @return Pointer to a newly created DensityFunction implementation. Memory
   * management for the pointer needs to be done by the calling routine.
   */
  static DensityFunction *generate(ParameterFile &params, Log *log = nullptr) {
    std::string type =
        params.get_value< std::string >("DensityFunction:type", "Homogeneous");
    if (log) {
      log->write_info("Requested DensityFunction type: ", type);
    }
#ifndef HAVE_HDF5
    check_hdf5(type, log);
#endif
    // there is some order here: first the non-library dependent
    // implementations, then the library dependent ones (sorted alphabetically
    // on library name). Each group is sorted alphabetically as well.
    if (type == "AsciiFile") {
      return new AsciiFileDensityFunction(params, log);
    } else if (type == "BlockSyntax") {
      return new BlockSyntaxDensityFunction(params, log);
    } else if (type == "BondiProfile") {
      return new BondiProfileDensityFunction(params);
    } else if (type == "CoredDMProfile") {
      return new CoredDMProfileDensityFunction(params);
    } else if (type == "DiscIC") {
      return new DiscICDensityFunction(params, log);
    } else if (type == "DiscPatch") {
      return new DiscPatchDensityFunction(params);
    } else if (type == "Homogeneous") {
      return new HomogeneousDensityFunction(params, log);
    } else if (type == "HDF5") {
      return new HDF5DensityFunction(params, log);
    } else if (type == "Interpolated") {
      return new InterpolatedDensityFunction(params, log);
    } else if (type == "PhantomSnapshot") {
      return new PhantomSnapshotDensityFunction(params, log);
    } else if (type == "SPHNGSnapshot") {
      return new SPHNGSnapshotDensityFunction(params, log);
    } else if (type == "SpiralGalaxy") {
      return new SpiralGalaxyDensityFunction(params, log);
#ifdef HAVE_HDF5
    } else if (type == "AmunSnapshot") {
      return new AmunSnapshotDensityFunction(params, log);
    } else if (type == "ArepoSnapshot") {
      return new ArepoSnapshotDensityFunction(params, log);
    } else if (type == "BufferedCMacIonizeSnapshot") {
      return new BufferedCMacIonizeSnapshotDensityFunction(params, log);
    } else if (type == "CMacIonizeSnapshot") {
      return new CMacIonizeSnapshotDensityFunction(params, log);
    } else if (type == "HDF5File") {
        return new HDF5DensityFunction(params,log);
    } else if (type == "FLASHSnapshot") {
      return new FLASHSnapshotDensityFunction(params, log);
    } else if (type == "GadgetSnapshot") {
      return new GadgetSnapshotDensityFunction(params, log);
#endif
    } else if (type == "None") {
      return nullptr;
    } else {
      cmac_error("Unknown DensityFunction type: \"%s\".", type.c_str());
      return nullptr;
    }
  }
};

#endif // DENSITYFUNCTIONFACTORY_HPP
