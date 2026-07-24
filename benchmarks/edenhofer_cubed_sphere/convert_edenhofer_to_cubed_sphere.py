#!/usr/bin/env python3
"""Resample the Edenhofer et al. (2023) mean map onto CMacIonize grids."""

import argparse
from pathlib import Path

import h5py
import healpy as hp
import numpy as np
from astropy.io import fits

PC_IN_M = 3.085677581491367e16
# User-supplied conversion: one native differential-extinction unit per pc
# corresponds to 2700 cm^-3.  CMacIonize uses SI number densities.
EDENHOFER_TO_M3 = 2700.0e6
DEFAULT_RESOLUTIONS = (64, 128, 256)
LOG_FLOOR = 1.0e-30

NORMAL = np.array(((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                   (0, -1, 0), (0, 0, 1), (0, 0, -1)), dtype=float)
U_AXIS = np.array(((0, 1, 0), (0, -1, 0), (-1, 0, 0),
                   (1, 0, 0), (0, 1, 0), (0, 1, 0)), dtype=float)
V_AXIS = np.array(((0, 0, 1), (0, 0, 1), (0, 0, 1),
                   (0, 0, 1), (-1, 0, 0), (1, 0, 0)), dtype=float)


def cubed_sphere_directions(n):
    """Return directions in CMac's flattened (face,u,v) cell order."""
    q = -1.0 + (np.arange(n) + 0.5) * 2.0 / n
    u, v = np.meshgrid(q, q, indexing="ij")
    directions = []
    for face in range(6):
        xyz = (NORMAL[face, :, None, None]
               + U_AXIS[face, :, None, None] * u
               + V_AXIS[face, :, None, None] * v)
        xyz /= np.linalg.norm(xyz, axis=0)
        directions.append(xyz.reshape(3, -1).T)
    return np.concatenate(directions)


def log_angular_interpolation(values, pixels, weights):
    """Match dustmaps' positive-field interpolation in logarithmic space."""
    logged = np.log(np.maximum(np.asarray(values, dtype=np.float64), LOG_FLOOR))
    return np.exp(np.sum(logged[..., pixels] * weights, axis=-2))


def convert(fits_path, output_dir, angular_n, radial_chunk):
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"edenhofer_cubed_n{angular_n}_r517.hdf5"
    directions = cubed_sphere_directions(angular_n)
    theta, phi = hp.vec2ang(directions)
    pixels, weights = hp.get_interp_weights(
        256, theta, phi, nest=True)

    with fits.open(fits_path, memmap=True) as hdus:
        source = hdus["MEAN"].data
        native_edges_pc = np.asarray(
            hdus["RADIAL PIXEL BOUNDARIES"].data.field(0), dtype=np.float64)
        inner = hdus["MEAN OF INTEGRATED INNER 68.8 PC"].data
        # One average inner cell preserves the supplied integrated inner column;
        # all remaining cells use the native radial boundaries unchanged.
        radial_edges_pc = np.concatenate(([0.1], native_edges_pc))

        with h5py.File(output, "w") as hdf:
            group = hdf.create_group("SphericalGrid")
            density = group.create_dataset(
                "NumberDensity",
                shape=(6, angular_n, angular_n, source.shape[0] + 1),
                dtype="f4", chunks=(1, min(32, angular_n),
                                    min(32, angular_n), 16),
                compression="gzip", compression_opts=1, shuffle=True)
            group.create_dataset("RadialEdges",
                                 data=radial_edges_pc * PC_IN_M)

            inner_differential = log_angular_interpolation(
                inner, pixels, weights) / native_edges_pc[0]
            density[..., 0] = (
                inner_differential.reshape(6, angular_n, angular_n)
                * EDENHOFER_TO_M3).astype(np.float32)

            for start in range(0, source.shape[0], radial_chunk):
                stop = min(start + radial_chunk, source.shape[0])
                interpolated = log_angular_interpolation(
                    source[start:stop], pixels, weights)
                block = interpolated.T.reshape(
                    6, angular_n, angular_n, stop - start)
                density[..., start + 1:stop + 1] = (
                    block * EDENHOFER_TO_M3).astype(np.float32)
                print(f"{output.name}: radial layers {start + 1}-{stop}"
                      f"/{source.shape[0]}", flush=True)

            group.attrs["AngularCellsPerFace"] = angular_n
            group.attrs["RadialCells"] = source.shape[0] + 1
            group.attrs["NumberDensityUnit"] = "m^-3"
            group.attrs["NativeConversionToCm-3"] = 2700.0
            group.attrs["AppliedConversionToM-3"] = EDENHOFER_TO_M3
            group.attrs["SourceNSIDE"] = int(hdus["MEAN"].header["NSIDE"])
            group.attrs["SourceOrdering"] = hdus["MEAN"].header["ORDERING"]
            group.attrs["AngularInterpolation"] = "bilinear in log(density)"
            group.attrs["Coordinates"] = (
                "Sun-centred Galactic Cartesian: +x GC, +y rotation, +z north")
            group.attrs["SourceFile"] = str(fits_path)
    print(f"Wrote {output}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fits_file", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--angular-resolution", type=int, nargs="+",
                        default=DEFAULT_RESOLUTIONS, metavar="N")
    parser.add_argument("--radial-chunk", type=int, default=8)
    args = parser.parse_args()
    for angular_n in args.angular_resolution:
        convert(args.fits_file, args.output_directory, angular_n,
                args.radial_chunk)


if __name__ == "__main__":
    main()
