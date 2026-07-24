#!/usr/bin/env python3
"""Plot the radial hydrogen column of a prepared cubed-sphere map."""

import argparse
from pathlib import Path

import h5py
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


def cubed_indices(directions, n):
    absolute = np.abs(directions)
    axis = np.argmax(absolute, axis=0)
    face = 2 * axis + (directions[axis, np.arange(directions.shape[1])] < 0)
    normal = np.array(((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                       (0, -1, 0), (0, 0, 1), (0, 0, -1)), dtype=float)
    u_axis = np.array(((0, 1, 0), (0, -1, 0), (-1, 0, 0),
                       (1, 0, 0), (0, 1, 0), (0, 1, 0)), dtype=float)
    v_axis = np.array(((0, 0, 1), (0, 0, 1), (0, 0, 1),
                       (0, 0, 1), (-1, 0, 0), (1, 0, 0)), dtype=float)
    selected_n = normal[face]
    denominator = np.sum(selected_n * directions.T, axis=1)
    u = np.sum(u_axis[face] * directions.T, axis=1) / denominator
    v = np.sum(v_axis[face] * directions.T, axis=1) / denominator
    iu = np.clip(((u + 1.0) * 0.5 * n).astype(int), 0, n - 1)
    iv = np.clip(((v + 1.0) * 0.5 * n).astype(int), 0, n - 1)
    return face, iu, iv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("map", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--display-nside", type=int, default=256)
    args = parser.parse_args()
    output = args.output or args.map.with_suffix(".sky.png")

    with h5py.File(args.map, "r") as hdf:
        group = hdf["SphericalGrid"]
        edges = group["RadialEdges"][:]
        density = group["NumberDensity"]
        n = density.shape[1]
        column = np.zeros(density.shape[:3], dtype=np.float64)
        widths = np.diff(edges)
        for start in range(0, density.shape[-1], 32):
            stop = min(start + 32, density.shape[-1])
            column += np.sum(
                density[..., start:stop] * widths[start:stop], axis=-1)

    pixels = np.arange(hp.nside2npix(args.display_nside))
    directions = np.asarray(
        hp.pix2vec(args.display_nside, pixels, nest=False))
    face, iu, iv = cubed_indices(directions, n)
    sky = np.log10(np.maximum(column[face, iu, iv] * 1.0e-4, 1.0))
    hp.mollview(sky, coord="G", unit=r"$\log_{10} N_H\ [{\rm cm}^{-2}]$",
                title=f"Edenhofer cubed sphere: N={n}", cmap="magma")
    hp.graticule()
    plt.savefig(output, dpi=180, bbox_inches="tight")
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
