#!/usr/bin/env python3
"""Create a top-down gas surface-density GIF from CMacIonize snapshots."""

import argparse
import csv
import glob
import gzip
import hashlib
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np


PC_IN_M = 3.0856775814913673e16
MYR_IN_S = 3.15576e13
CACHE_VERSION = 1

# Display controls for supernova circles.
SN_CIRCLE_SPEED_KMS = 10.0
SN_CIRCLE_FADE_TIME_MYR = 1.0


def luminosity_from_mass(mass):
    """Return Q_H0 using the same table and interpolation as CMacIonize."""
    masses = np.array(
        [57.95, 46.94, 38.08, 34.39, 30.98, 28.0,
         25.29, 22.90, 20.76, 18.80, 17.08, 15.55]
    )
    log_luminosities = np.array(
        [49.64, 49.44, 49.22, 49.10, 48.99, 48.88,
         48.75, 48.61, 48.44, 48.27, 48.06, 47.88]
    )
    if mass < masses[-1]:
        return 0.0
    if mass >= masses[0]:
        return 10.0 ** log_luminosities[0]
    return 10.0 ** np.interp(mass, masses[::-1], log_luminosities[::-1])


def load_stellar_history(filename):
    """Read the Swiggum tracks and retain one trajectory per star."""
    opener = gzip.open if str(filename).endswith(".gz") else open
    stars = {}
    with opener(filename, "rt", newline="") as source:
        for row in csv.DictReader(source):
            star_id = int(row["star_id"])
            star = stars.setdefault(
                star_id,
                {
                    "birth": float(row["stellar_birth_time_myr"]),
                    "death": float(row["SNTime"]),
                    "mass": float(row["Stellar_Mass"]),
                    "track": [],
                },
            )
            star["track"].append(
                (float(row["trajectory_time_myr"]),
                 float(row["X"]), float(row["Y"]), float(row["Z"]))
            )
    for star in stars.values():
        track = sorted(set(star["track"]))
        star["times"] = np.array([point[0] for point in track])
        star["positions"] = np.array([point[1:] for point in track])
        star["luminosity"] = luminosity_from_mass(star["mass"])
        del star["track"]
    return list(stars.values())


def interpolate_position(star, history_time):
    """Linearly interpolate a catalogue track as the C++ reader does."""
    times = star["times"]
    positions = star["positions"]
    if history_time <= times[0]:
        return positions[0]
    if history_time >= times[-1]:
        return positions[-1]
    right = np.searchsorted(times, history_time, side="right")
    fraction = ((history_time - times[right - 1]) /
                (times[right] - times[right - 1]))
    return positions[right - 1] + fraction * (
        positions[right] - positions[right - 1]
    )


def point_in_box(position, half_sides_pc):
    return np.all(position >= -half_sides_pc) and np.all(position <= half_sides_pc)


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


def read_surface_density(filename, subgrids, grid_size, cache_directory):
    """Read or build a cached top-down projection for one snapshot."""
    snapshot_path = Path(filename).resolve()
    stat = snapshot_path.stat()
    cache_key = repr(
        (CACHE_VERSION, str(snapshot_path), stat.st_size, stat.st_mtime_ns,
         subgrids, grid_size)
    ).encode("utf-8")
    cache_file = cache_directory / (hashlib.sha256(cache_key).hexdigest() + ".npz")
    if cache_file.exists():
        try:
            with np.load(cache_file) as cached:
                return (cached["surface_density"], float(cached["time_myr"]),
                        cached["box_size"])
        except (OSError, ValueError, KeyError):
            cache_file.unlink(missing_ok=True)

    with h5py.File(filename, "r") as snapshot:
        density = snapshot["PartType0/Density"][:]
        box_size = np.asarray(snapshot["Header"].attrs["BoxSize"])
        time_myr = float(snapshot["Header"].attrs["Time"]) / MYR_IN_S
    density_grid = get_tb_grid(density, subgrids, grid_size)
    dz = box_size[2] / grid_size[2]
    surface_density = density_grid.sum(axis=2) * dz
    np.savez_compressed(
        cache_file, surface_density=surface_density, time_myr=time_myr,
        box_size=box_size
    )
    return surface_density, time_myr, box_size


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pattern", help="quoted snapshot glob")
    parser.add_argument("output", help="output GIF filename")
    parser.add_argument("--grid-size", nargs=3, type=int, required=True)
    parser.add_argument("--subgrids", nargs=3, type=int, required=True)
    parser.add_argument("--fps", type=int, default=8)
    parser.add_argument("--vmin", type=float, default=None)
    parser.add_argument("--vmax", type=float, default=None)
    parser.add_argument(
        "--track-stars", action="store_true",
        help="overlay living stars and supernova remnants from the Swiggum history",
    )
    parser.add_argument(
        "--star-history", type=Path,
        default=Path(__file__).resolve().parents[2] / "data" /
        "sn_history_july9_imf_sample0_minus100_to_plus20myr.csv.gz",
        help="stellar-history CSV or CSV.GZ used with --track-stars",
    )
    parser.add_argument(
        "--cache-dir", type=Path, default=Path(".swiggum_plot_cache"),
        help="directory for reusable 2-D projections",
    )
    args = parser.parse_args()

    filenames = sorted(glob.glob(args.pattern))
    if not filenames:
        raise SystemExit(f"No snapshots match {args.pattern!r}")

    grid_size = tuple(args.grid_size)
    subgrids = tuple(args.subgrids)
    if any(grid_size[i] % subgrids[i] for i in range(3)):
        raise SystemExit("Grid dimensions must be divisible by subgrid dimensions")
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    stars = load_stellar_history(args.star_history) if args.track_stars else []

    if args.vmin is None or args.vmax is None:
        low, high = np.inf, -np.inf
        for filename in filenames:
            surface_density, _, _ = read_surface_density(
                filename, subgrids, grid_size, args.cache_dir
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
        filenames[0], subgrids, grid_size, args.cache_dir
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

    star_artist = axis.scatter(
        [], [], marker="x", c="deepskyblue", linewidths=1.2, zorder=3
    )
    remnant_artists = []

    def update_stars(history_time):
        nonlocal remnant_artists
        positions = []
        sizes = []
        half_sides_pc = 0.5 * box_size / PC_IN_M
        for star in stars:
            if (star["birth"] <= history_time < star["death"] and
                    star["luminosity"] > 0.0):
                position = interpolate_position(star, history_time)
                if point_in_box(position, half_sides_pc):
                    positions.append(position[:2] / 1000.0)
                    # Scatter sizes are areas: this makes marker diameter scale
                    # approximately as Q_H0^(1/4), without huge O-star crosses.
                    sizes.append(24.0 * (star["luminosity"] / 1.0e48) ** 0.5)
        star_artist.set_offsets(np.asarray(positions).reshape(-1, 2))
        star_artist.set_sizes(sizes)

        for artist in remnant_artists:
            artist.remove()
        remnant_artists = []
        for star in stars:
            age_myr = history_time - star["death"]
            if 0.0 <= age_myr < SN_CIRCLE_FADE_TIME_MYR:
                position = interpolate_position(star, star["death"])
                if point_in_box(position, half_sides_pc):
                    radius_kpc = (SN_CIRCLE_SPEED_KMS * 1.0e3 * age_myr *
                                  MYR_IN_S / (1000.0 * PC_IN_M))
                    remnant = plt.Circle(
                        position[:2] / 1000.0, radius_kpc, fill=False,
                        edgecolor="cyan", linewidth=1.2,
                        alpha=1.0 - age_myr / SN_CIRCLE_FADE_TIME_MYR,
                        zorder=3,
                    )
                    axis.add_patch(remnant)
                    remnant_artists.append(remnant)

    def update(frame):
        surface_density, simulation_time, _ = read_surface_density(
            filenames[frame], subgrids, grid_size, args.cache_dir
        )
        image.set_data(
            np.log10(np.maximum(surface_density.T, np.finfo(float).tiny))
        )
        title.set_text(f"history time = {simulation_time - 100.0:.2f} Myr")
        if args.track_stars:
            update_stars(simulation_time - 100.0)
        return image, title, star_artist, *remnant_artists

    animation = FuncAnimation(
        figure, update, frames=len(filenames), interval=1000.0 / args.fps,
        blit=False
    )
    animation.save(args.output, writer=PillowWriter(fps=args.fps), dpi=120)


if __name__ == "__main__":
    main()
