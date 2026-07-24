#!/usr/bin/env python3
"""Plot an HI or HII column from a spherical CMacIonize Gadget output."""

import argparse
from pathlib import Path

import h5py
import healpy as hp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def cubed_coordinates(points):
    axis = np.argmax(np.abs(points), axis=1)
    face = 2 * axis + (points[np.arange(len(points)), axis] < 0)
    normal = np.array(((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                       (0, -1, 0), (0, 0, 1), (0, 0, -1)), dtype=float)
    u_axis = np.array(((0, 1, 0), (0, -1, 0), (-1, 0, 0),
                       (1, 0, 0), (0, 1, 0), (0, 1, 0)), dtype=float)
    v_axis = np.array(((0, 0, 1), (0, 0, 1), (0, 0, 1),
                       (0, 0, 1), (-1, 0, 0), (1, 0, 0)), dtype=float)
    denominator = np.sum(normal[face] * points, axis=1)
    u = np.sum(u_axis[face] * points, axis=1) / denominator
    v = np.sum(v_axis[face] * points, axis=1) / denominator
    return face, u, v


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("snapshot", type=Path,
                        help="CMacIonize Gadget HDF5 output")
    parser.add_argument("input_map", type=Path,
                        help="Edenhofer cubed-sphere input used by the run")
    parser.add_argument("--component", choices=("ionized", "neutral", "total"),
                        default="ionized")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--display-nside", type=int, default=256)
    parser.add_argument("--chunk", type=int, default=1_000_000)
    args = parser.parse_args()
    output = args.output or args.snapshot.with_name(
        args.snapshot.stem + f".{args.component}_sky.png")

    with h5py.File(args.input_map, "r") as source:
        edges = source["SphericalGrid/RadialEdges"][:]
        n = source["SphericalGrid/NumberDensity"].shape[1]
    widths = np.diff(edges)
    column = np.zeros((6, n, n), dtype=np.float64)

    with h5py.File(args.snapshot, "r") as snapshot:
        cells = snapshot["PartType0"]
        coordinates = cells["Coordinates"]
        density = cells["NumberDensity"]
        neutral_fraction = cells.get("NeutralFractionH")
        if args.component != "total" and neutral_fraction is None:
            raise KeyError("Snapshot does not contain PartType0/NeutralFractionH")
        # Gadget output coordinates are relative to the box anchor. The
        # spherical centre is the centre of the supporting Cartesian box.
        centre = 0.5 * np.asarray(snapshot["Header"].attrs["BoxSize"])
        for start in range(0, len(coordinates), args.chunk):
            stop = min(start + args.chunk, len(coordinates))
            points = np.asarray(coordinates[start:stop]) - centre
            radius = np.linalg.norm(points, axis=1)
            face, u, v = cubed_coordinates(points)
            iu = np.clip(((u + 1.0) * 0.5 * n).astype(int), 0, n - 1)
            iv = np.clip(((v + 1.0) * 0.5 * n).astype(int), 0, n - 1)
            ir = np.clip(np.searchsorted(edges, radius, side="right") - 1,
                         0, len(widths) - 1)
            values = np.asarray(density[start:stop])
            if args.component != "total":
                neutral = np.asarray(neutral_fraction[start:stop])
                values *= neutral if args.component == "neutral" else 1.0-neutral
            np.add.at(column, (face, iu, iv), values * widths[ir])
            print(f"Read {stop}/{len(coordinates)} cells", flush=True)

    pixels = np.arange(hp.nside2npix(args.display_nside))
    directions = np.asarray(hp.pix2vec(args.display_nside, pixels)).T
    face, u, v = cubed_coordinates(directions)
    iu = np.clip(((u + 1.0) * 0.5 * n).astype(int), 0, n - 1)
    iv = np.clip(((v + 1.0) * 0.5 * n).astype(int), 0, n - 1)
    sky = np.log10(np.maximum(column[face, iu, iv] * 1.0e-4, 1.0))
    label = {"ionized": "HII", "neutral": "HI", "total": "H"}[args.component]
    hp.mollview(sky, coord="G", cmap="magma",
                unit=rf"$\log_{{10}} N_{{\rm {label}}}\ "
                     r"[{\rm cm}^{-2}]$",
                title=f"CMacIonize {label} column")
    hp.graticule()
    plt.savefig(output, dpi=180, bbox_inches="tight")
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
