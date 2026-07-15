#!/usr/bin/env python3
"""Create a top-down gas surface-density GIF from CMacIonize snapshots."""

import argparse
import glob

import h5py
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np


def get_tb_grid(grid, subgriddim, gridsize):
    """Restore task-based subgrid output ordering to a Cartesian array."""
    result = np.zeros((gridsize[0], gridsize[1], gridsize[2]))
    cx = int(gridsize[0] / subgriddim[0])
    cy = int(gridsize[1] / subgriddim[1])
    cz = int(gridsize[2] / subgriddim[2])
    startchunk = 0
    endchunk = cx * cy * cz
    ix = 0
    iy = 0
    iz = 0
    while endchunk <= gridsize[0] * gridsize[1] * gridsize[2]:
        chunk = np.asarray(grid[startchunk:endchunk])
        result[ix : ix + cx, iy : iy + cy, iz : iz + cz] = chunk.reshape(
            cx, cy, cz
        )
        startchunk += cx * cy * cz
        endchunk += cx * cy * cz
        iz += cz
        if iz == gridsize[2]:
            iz = 0
            iy += cy
            if iy == gridsize[1]:
                iy = 0
                ix += cx
    return result


def read_surface_density(filename, subgrids, grid_size):
    with h5py.File(filename, "r") as snapshot:
        density = snapshot["PartType0/Density"][:]
        box_size = np.asarray(snapshot["Header"].attrs["BoxSize"])
        time_myr = float(snapshot["Header"].attrs["Time"]) / 3.15576e13
    density_grid = get_tb_grid(density, subgrids, grid_size)
    dz = box_size[2] / grid_size[2]
    return density_grid.sum(axis=2) * dz, time_myr, box_size


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pattern", help="quoted snapshot glob")
    parser.add_argument("output", help="output GIF filename")
    parser.add_argument("--grid-size", nargs=3, type=int, required=True)
    parser.add_argument("--subgrids", nargs=3, type=int, required=True)
    parser.add_argument("--fps", type=int, default=8)
    parser.add_argument("--vmin", type=float, default=None)
    parser.add_argument("--vmax", type=float, default=None)
    args = parser.parse_args()

    filenames = sorted(glob.glob(args.pattern))
    if not filenames:
        raise SystemExit(f"No snapshots match {args.pattern!r}")

    grid_size = tuple(args.grid_size)
    subgrids = tuple(args.subgrids)
    if any(grid_size[i] % subgrids[i] for i in range(3)):
        raise SystemExit("Grid dimensions must be divisible by subgrid dimensions")

    if args.vmin is None or args.vmax is None:
        low, high = np.inf, -np.inf
        for filename in filenames:
            surface_density, _, _ = read_surface_density(
                filename, subgrids, grid_size
            )
            positive = surface_density[surface_density > 0.0]
            if positive.size:
                low = min(low, np.percentile(positive, 1.0))
                high = max(high, np.percentile(positive, 99.5))
        vmin = np.log10(args.vmin if args.vmin is not None else low)
        vmax = np.log10(args.vmax if args.vmax is not None else high)
    else:
        vmin, vmax = np.log10(args.vmin), np.log10(args.vmax)

    first, first_time, box_size = read_surface_density(
        filenames[0], subgrids, grid_size
    )
    kpc = 3.0856775814913673e19
    extent = np.array([-0.5 * box_size[0], 0.5 * box_size[0],
                       -0.5 * box_size[1], 0.5 * box_size[1]]) / kpc

    figure, axis = plt.subplots(figsize=(7, 6))
    image = axis.imshow(
        np.log10(np.maximum(first.T, np.finfo(float).tiny)),
        origin="lower",
        extent=extent,
        vmin=vmin,
        vmax=vmax,
        cmap="hot",
        interpolation="nearest",
    )
    title = axis.set_title(f"history time = {first_time - 100.0:.2f} Myr")
    axis.set_xlabel("X toward Galactic centre [kpc]")
    axis.set_ylabel("Y along Galactic rotation [kpc]")
    figure.colorbar(image, ax=axis, label=r"$\log_{10}\,\Sigma$ [kg m$^{-2}$]")
    figure.tight_layout()

    def update(frame):
        surface_density, simulation_time, _ = read_surface_density(
            filenames[frame], subgrids, grid_size
        )
        image.set_data(
            np.log10(np.maximum(surface_density.T, np.finfo(float).tiny))
        )
        title.set_text(f"history time = {simulation_time - 100.0:.2f} Myr")
        return image, title

    animation = FuncAnimation(
        figure, update, frames=len(filenames), interval=1000.0 / args.fps,
        blit=False
    )
    animation.save(args.output, writer=PillowWriter(fps=args.fps), dpi=120)


if __name__ == "__main__":
    main()
