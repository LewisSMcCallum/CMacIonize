#!/usr/bin/env python3
"""Make Cartesian and all-sky diagnostics for one Swiggum snapshot.

The observer is at the box centre.  HEALPix longitude zero points along +X
(towards the Galactic centre), longitude 90 degrees points along +Y (Galactic
rotation), and +Z is the north pole.
"""

import argparse
from pathlib import Path

import h5py
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Ellipse
import numpy as np

try:
    from scipy.ndimage import map_coordinates
except ImportError:  # pragma: no cover - nearest-neighbour fallback is tested below
    map_coordinates = None


PC_IN_M = 3.0856775814913673e16
KPC_IN_M = 1000.0 * PC_IN_M
M_IN_CM = 100.0
MYR_IN_S = 3.15576e13
KEV_IN_J = 1.602176634e-16
BOLTZMANN = 1.380649e-23

# Approximate eROSITA bands in keV.  The emission model below includes only
# hydrogen free-free continuum, so soft-band metal-line emission is absent.
XRAY_BANDS_KEV = ((0.2, 0.6), (0.6, 2.3), (2.3, 5.0))

# Plot-only choices.  The percentile clipping prevents a few cells from
# stretching the colour bars over many unhelpful decades.  The maximum log
# span can be changed without changing the data or the ray integration.
LOG_LOWER_PERCENTILE = 2.0
LOG_UPPER_PERCENTILE = 98.0
LOG_MAX_DECADES = 5.0
LINEAR_LOWER_PERCENTILE = 1.0
LINEAR_UPPER_PERCENTILE = 99.0

# Approximate guide regions in Galactic longitude/latitude.  These are labels
# for orientation, not extra data in the simulation; edit them if a preferred
# literature footprint is desired.
SKY_REGIONS = (
    {"name": "Orion–Eridanus", "lon": 205.0, "lat": -18.0,
     "width": 28.0, "height": 18.0, "angle": -20.0,
     "colour": "cyan"},
    {"name": "Gum", "lon": 258.0, "lat": -2.0,
     "width": 24.0, "height": 18.0, "angle": 0.0,
     "colour": "deepskyblue"},
)

# Approximate 3-D guide regions for the Cartesian panels.  Distances and
# radii are deliberately exposed here because the plotted regions are only
# contextual overlays, not inputs to CMacIonize.
CARTESIAN_REGIONS = (
    {"name": "Orion–Eridanus", "lon": 205.0, "lat": -18.0,
     "distance_pc": 400.0, "radius_pc": 250.0, "colour": "cyan"},
    {"name": "Gum", "lon": 258.0, "lat": -2.0,
     "distance_pc": 450.0, "radius_pc": 140.0, "colour": "deepskyblue"},
)


def get_tb_grid(values, subgrids, grid_size):
    """Restore task-based subgrid ordering to an X,Y,Z Cartesian array."""
    result = np.empty(grid_size, dtype=np.float32)
    chunk_size = tuple(grid_size[i] // subgrids[i] for i in range(3))
    cells_per_chunk = int(np.prod(chunk_size))
    offset = 0
    for ix in range(0, grid_size[0], chunk_size[0]):
        for iy in range(0, grid_size[1], chunk_size[1]):
            for iz in range(0, grid_size[2], chunk_size[2]):
                result[
                    ix : ix + chunk_size[0],
                    iy : iy + chunk_size[1],
                    iz : iz + chunk_size[2],
                ] = np.asarray(values[offset : offset + cells_per_chunk]).reshape(
                    chunk_size
                )
                offset += cells_per_chunk
    return result


def read_grid(snapshot, name, subgrids, grid_size):
    path = f"PartType0/{name}"
    if path not in snapshot:
        raise SystemExit(
            f"Snapshot does not contain {path}; enable {name} in "
            "DensityGridWriterFields."
        )
    return get_tb_grid(snapshot[path], subgrids, grid_size)


def positive_limits(values, lower=LOG_LOWER_PERCENTILE,
                   upper=LOG_UPPER_PERCENTILE, max_decades=LOG_MAX_DECADES):
    array = np.asarray(values)
    finite = array[np.isfinite(array) & (array > 0.0)]
    if finite.size == 0:
        return 1.0e-30, 1.0
    low, high = np.percentile(finite, (lower, upper))
    low = max(low, high / 10.0 ** max_decades)
    if high <= low:
        high = low * 10.0
    return low, high


def linear_limits(values):
    array = np.asarray(values)
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return 0.0, 1.0
    low, high = np.percentile(
        finite, (LINEAR_LOWER_PERCENTILE, LINEAR_UPPER_PERCENTILE)
    )
    if high <= low:
        high = low + (abs(low) if low != 0.0 else 1.0)
    return low, high


def show_map(axis, values, extent, title, label, cmap="magma", linear=False,
             xlabel="X toward Galactic centre [kpc]",
             ylabel="Y along Galactic rotation [kpc]", aspect=None):
    if linear:
        low, high = linear_limits(values)
        norm = Normalize(low, high)
    else:
        low, high = positive_limits(values)
        norm = LogNorm(low, high)
    image = axis.imshow(np.ma.masked_invalid(np.asarray(values).T),
                        origin="lower", extent=extent, cmap=cmap, norm=norm,
                        interpolation="nearest", aspect=aspect)
    axis.set_title(title)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    plt.colorbar(image, ax=axis, label=label)


def add_cartesian_regions(axis, plane, box_size_m):
    """Add approximate Orion–Eridanus and Gum footprints to an XY/XZ panel."""
    half = 0.5 * np.asarray(box_size_m) / KPC_IN_M
    for region in CARTESIAN_REGIONS:
        lon = np.radians(region["lon"])
        lat = np.radians(region["lat"])
        distance = region["distance_pc"] / 1000.0
        centre = np.array([
            distance * np.cos(lat) * np.cos(lon),
            distance * np.cos(lat) * np.sin(lon),
            distance * np.sin(lat),
        ])
        radius = 2.0 * region["radius_pc"] / 1000.0
        if plane == "xy":
            x, y = centre[0], centre[1]
            width, height = radius, radius
        else:
            x, y = centre[0], centre[2]
            width, height = radius, radius
        if (x + width / 2.0 < -half[0] or x - width / 2.0 > half[0] or
                y + height / 2.0 < -half[1 if plane == "xy" else 2] or
                y - height / 2.0 > half[1 if plane == "xy" else 2]):
            continue
        axis.add_patch(Ellipse((x, y), width, height, fill=False,
                               linewidth=1.5, edgecolor=region["colour"]))
        axis.text(x, y, region["name"], color=region["colour"],
                  fontsize=10, ha="center", va="center",
                  bbox={"facecolor": "black", "alpha": 0.35, "pad": 1.5})


def free_free_emissivity(ne_cm3, ni_cm3, temperature, band):
    """Hydrogen free-free band emissivity in erg s^-1 cm^-3.

    This integrates a constant-Gaunt-factor bremsstrahlung spectrum over the
    requested energy band.  It deliberately excludes metal lines and absorption.
    """
    safe_temperature = np.maximum(temperature, 1.0)
    kT = BOLTZMANN * safe_temperature
    lower = np.exp(-band[0] * KEV_IN_J / kT)
    upper = np.exp(-band[1] * KEV_IN_J / kT)
    return (1.426e-27 * 1.2 * ne_cm3 * ni_cm3 *
            np.sqrt(safe_temperature) * (lower - upper))


def _sample_cartesian(field, coordinates):
    """Sample a Cartesian field at fractional cell-centre coordinates."""
    if map_coordinates is not None:
        return map_coordinates(field, coordinates, order=1, mode="nearest",
                               cval=0.0, prefilter=False)
    # Keep the plotting script usable in a minimal Python environment.  This
    # is the same trilinear interpolation as map_coordinates(order=1), with
    # constant zero outside the box.
    # A cell is centred on index i but occupies i +/- 0.5, so coordinates
    # near a box face should use the boundary cell rather than be discarded.
    valid = np.all((coordinates >= -0.5) &
                   (coordinates <= (np.asarray(field.shape)[:, None] - 0.5)),
                   axis=0)
    clipped = np.clip(coordinates, 0.0,
                      np.asarray(field.shape)[:, None] - 1.0)
    lower = np.floor(clipped).astype(np.int64)
    upper = np.minimum(lower + 1, np.asarray(field.shape)[:, None] - 1)
    fraction = clipped - lower
    x0, y0, z0 = lower
    x1, y1, z1 = upper
    fx, fy, fz = fraction
    c000 = field[x0, y0, z0]
    c100 = field[x1, y0, z0]
    c010 = field[x0, y1, z0]
    c110 = field[x1, y1, z0]
    c001 = field[x0, y0, z1]
    c101 = field[x1, y0, z1]
    c011 = field[x0, y1, z1]
    c111 = field[x1, y1, z1]
    sampled = ((1.0 - fx) * (1.0 - fy) * (1.0 - fz) * c000 +
               fx * (1.0 - fy) * (1.0 - fz) * c100 +
               (1.0 - fx) * fy * (1.0 - fz) * c010 +
               fx * fy * (1.0 - fz) * c110 +
               (1.0 - fx) * (1.0 - fy) * fz * c001 +
               fx * (1.0 - fy) * fz * c101 +
               (1.0 - fx) * fy * fz * c011 +
               fx * fy * fz * c111)
    sampled = np.asarray(sampled)
    sampled[~valid] = 0.0
    return sampled


def make_healpix_maps(hi_cm3, halpha_proxy, box_size, nside,
                      radial_oversampling=2, pixel_chunk=2048):
    """Ray-integrate the Cartesian fields onto a centre-observer HEALPix sky.

    Each HEALPix pixel is represented by a ray from the box centre.  The
    fields are linearly sampled at sub-cell radial intervals and integrated to
    the cube boundary.  This produces H I columns (cm^-2) and an Halpha
    emission-measure proxy (cm^-6 pc K^-0.9), rather than the old cell-volume
    / r^2 deposition which was not a line-of-sight sky map.
    """
    grid_size = np.asarray(hi_cm3.shape, dtype=int)
    box_size = np.asarray(box_size, dtype=float)
    cell_size_m = box_size / grid_size
    dr = np.min(cell_size_m) / max(int(radial_oversampling), 1)
    npix = hp.nside2npix(nside)
    directions = np.column_stack(hp.pix2vec(nside, np.arange(npix)))
    half_box = 0.5 * box_size
    # Distance to the first cube face along each ray.
    ray_limits = np.full(npix, np.inf)
    for axis in range(3):
        absolute_direction = np.abs(directions[:, axis])
        candidate = np.divide(half_box[axis], absolute_direction,
                              out=np.full(npix, np.inf),
                              where=absolute_direction > 0.0)
        ray_limits = np.minimum(ray_limits, candidate)
    nsteps = int(np.ceil(np.max(ray_limits) / dr))
    radii = (np.arange(nsteps, dtype=float) + 0.5) * dr
    hi_map = np.zeros(npix, dtype=float)
    halpha_map = np.zeros(npix, dtype=float)
    chunk = max(int(pixel_chunk), 1)
    for start in range(0, npix, chunk):
        stop = min(start + chunk, npix)
        direction_chunk = directions[start:stop]
        valid = radii[None, :] < ray_limits[start:stop, None]
        positions = direction_chunk[:, :, None] * radii[None, None, :]
        coordinates = ((positions + half_box[None, :, None]) /
                       cell_size_m[None, :, None] - 0.5)
        coordinates = coordinates.transpose(1, 0, 2).reshape(3, -1)
        hi_values = _sample_cartesian(hi_cm3, coordinates).reshape(stop - start,
                                                                    nsteps)
        halpha_values = _sample_cartesian(halpha_proxy, coordinates).reshape(
            stop - start, nsteps
        )
        hi_values[~valid] = 0.0
        halpha_values[~valid] = 0.0
        hi_map[start:stop] = hi_values.sum(axis=1) * dr * M_IN_CM
        halpha_map[start:stop] = halpha_values.sum(axis=1) * dr / PC_IN_M
    return halpha_map, hi_map


def _plot_cartesian_grid(prefix, panels, box_size, time_myr, linear, suffix):
    mode = "linear" if linear else "log"
    figure, axes = plt.subplots(3, 4, figsize=(20, 15), constrained_layout=True)
    for row, row_panels in enumerate(panels):
        for column, panel in enumerate(row_panels):
            values, title, label, cmap, extent, xlabel, ylabel, plane = panel
            show_map(axes[row, column], values, extent, title, label, cmap,
                     linear=linear, xlabel=xlabel, ylabel=ylabel)
            if plane is not None:
                add_cartesian_regions(axes[row, column], plane, box_size)
    figure.suptitle(f"Simulation time {time_myr:.3f} Myr — {mode} colour scaling")
    figure.savefig(f"{prefix}_{suffix}.png", dpi=160)
    plt.close(figure)


def plot_cartesian(prefix, density, hplus_cm3, halpha, temperature, box_size,
                   time_myr):
    box_size = np.asarray(box_size)
    cell_size = box_size / np.asarray(density.shape)
    dz_m, dy_m = cell_size[2], cell_size[1]
    dz_cm, dy_cm = dz_m * M_IN_CM, dy_m * M_IN_CM
    dz_pc, dy_pc = dz_m / PC_IN_M, dy_m / PC_IN_M
    extent_xy = (-0.5 * box_size[0] / KPC_IN_M, 0.5 * box_size[0] / KPC_IN_M,
                 -0.5 * box_size[1] / KPC_IN_M, 0.5 * box_size[1] / KPC_IN_M)
    extent_xz = (-0.5 * box_size[0] / KPC_IN_M, 0.5 * box_size[0] / KPC_IN_M,
                 -0.5 * box_size[2] / KPC_IN_M, 0.5 * box_size[2] / KPC_IN_M)
    midplane = density.shape[2] // 2

    panels = [
        [
            (density.sum(axis=2) * dz_m, "Face-on gas surface density",
             r"$\Sigma$ [kg m$^{-2}$]", "magma", extent_xy,
             "X toward Galactic centre [kpc]", "Y along Galactic rotation [kpc]", "xy"),
            (hplus_cm3.sum(axis=2) * dz_cm, "Face-on ionized-H column",
             r"$N_{H^+}$ [cm$^{-2}$]", "viridis", extent_xy,
             "X [kpc]", "Y [kpc]", "xy"),
            (halpha.sum(axis=2) * dz_pc, r"Face-on H$\alpha$ proxy",
             r"$\int n_e^2 T^{-0.9} dz$ [cm$^{-6}$ pc K$^{-0.9}$]",
             "inferno", extent_xy, "X [kpc]", "Y [kpc]", "xy"),
            (np.max(temperature, axis=2), "Face-on maximum temperature",
             "T [K]", "plasma", extent_xy, "X [kpc]", "Y [kpc]", "xy"),
        ],
        [
            (density.sum(axis=1) * dy_m, "Edge-on gas surface density",
             r"$\Sigma_y$ [kg m$^{-2}$]", "magma", extent_xz,
             "X toward Galactic centre [kpc]", "Z north [kpc]", "xz"),
            (hplus_cm3.sum(axis=1) * dy_cm, "Edge-on ionized-H column",
             r"$N_{H^+,y}$ [cm$^{-2}$]", "viridis", extent_xz,
             "X [kpc]", "Z [kpc]", "xz"),
            (halpha.sum(axis=1) * dy_pc, r"Edge-on H$\alpha$ proxy",
             r"$\int n_e^2 T^{-0.9} dy$ [cm$^{-6}$ pc K$^{-0.9}$]",
             "inferno", extent_xz, "X [kpc]", "Z [kpc]", "xz"),
            (np.max(temperature, axis=1), "Edge-on maximum temperature",
             "T [K]", "plasma", extent_xz, "X [kpc]", "Z [kpc]", "xz"),
        ],
        [
            (density[:, :, midplane] * 1.0e-3, "Midplane gas density",
             r"$\rho$ [g cm$^{-3}$]", "magma", extent_xy,
             "X [kpc]", "Y [kpc]", "xy"),
            (hplus_cm3[:, :, midplane], "Midplane ionized-H density",
             r"$n_{H^+}$ [cm$^{-3}$]", "viridis", extent_xy,
             "X [kpc]", "Y [kpc]", "xy"),
            (halpha[:, :, midplane], r"Midplane H$\alpha$ proxy",
             r"$n_e^2 T^{-0.9}$", "inferno", extent_xy,
             "X [kpc]", "Y [kpc]", "xy"),
            (temperature[:, :, midplane], "Midplane temperature", "T [K]",
             "plasma", extent_xy, "X [kpc]", "Y [kpc]", "xy"),
        ],
    ]
    _plot_cartesian_grid(prefix, panels, box_size, time_myr, False, "cartesian")
    _plot_cartesian_grid(prefix, panels, box_size, time_myr, True,
                         "cartesian_linear")


def add_sky_regions():
    theta = np.linspace(0.0, 2.0 * np.pi, 361)
    for region in SKY_REGIONS:
        angle = np.radians(region["angle"])
        x = 0.5 * region["width"] * np.cos(theta)
        y = 0.5 * region["height"] * np.sin(theta)
        rotated_x = x * np.cos(angle) - y * np.sin(angle)
        rotated_y = x * np.sin(angle) + y * np.cos(angle)
        lon = region["lon"] + rotated_x / np.cos(np.radians(region["lat"]))
        lat = region["lat"] + rotated_y
        hp.projplot(lon, lat, lonlat=True, color=region["colour"],
                    linewidth=1.5)
        hp.projtext(region["lon"], region["lat"], region["name"],
                    lonlat=True, color=region["colour"], fontsize=12,
                    ha="center", va="center")


def plot_sky(prefix, halpha_map, hi_map, nside, time_myr):
    figure = plt.figure(figsize=(14, 9))
    for panel, (values, title, unit) in enumerate(
        ((halpha_map, r"H$\alpha$ proxy sky", r"cm$^{-5}$ K$^{-0.9}$"),
         (hi_map, "H I column sky", r"cm$^{-2}$")), start=1
    ):
        low, high = positive_limits(values)
        displayed = values.copy()
        displayed[displayed <= 0.0] = hp.UNSEEN
        hp.mollview(
            displayed, fig=figure.number, sub=(2, 1, panel), title=title,
            unit=unit, norm="log", min=low, max=high, cmap="inferno",
            coord=None, notext=False, rot=180
        )
        add_sky_regions()
        hp.graticule(verbose=False)
    figure.suptitle(
        f"Observer at box centre; simulation time {time_myr:.3f} Myr\n"
        r"longitude 0$^\circ$ = +X (Galactic centre), 90$^\circ$ = +Y"
    )
    figure.savefig(f"{prefix}_sky.png", dpi=160)
    plt.close(figure)


def plot_xrays(prefix, ne_cm3, ni_cm3, temperature, box_size, time_myr):
    dz_cm = box_size[2] / temperature.shape[2] * M_IN_CM
    extent = (-0.5 * box_size[0] / KPC_IN_M, 0.5 * box_size[0] / KPC_IN_M,
              -0.5 * box_size[1] / KPC_IN_M, 0.5 * box_size[1] / KPC_IN_M)
    figure, axes = plt.subplots(1, 3, figsize=(18, 5.5), constrained_layout=True)
    for axis, band in zip(axes, XRAY_BANDS_KEV):
        emissivity = free_free_emissivity(ne_cm3, ni_cm3, temperature, band)
        projection = emissivity.sum(axis=2) * dz_cm
        show_map(
            axis, projection, extent,
            f"{band[0]:g}–{band[1]:g} keV free–free estimate",
            r"$\int\epsilon_{ff} dz$ [erg s$^{-1}$ cm$^{-2}$]", "magma",
        )
    figure.suptitle(
        f"Simulation time {time_myr:.3f} Myr — intrinsic hydrogen continuum; "
        "no metal lines or absorption"
    )
    figure.savefig(f"{prefix}_xray.png", dpi=160)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path)
    parser.add_argument(
        "--output-prefix", type=Path, default=None,
        help="output path prefix (default: snapshot filename without suffix)",
    )
    parser.add_argument("--grid-size", nargs=3, type=int, required=True)
    parser.add_argument("--subgrids", nargs=3, type=int, required=True)
    parser.add_argument(
        "--nside", type=int, default=64,
        help="HEALPix NSIDE (default: 64; use 128 for a finer sky)",
    )
    parser.add_argument(
        "--sky-radial-oversampling", type=int, default=2,
        help="samples per smallest Cartesian cell along each sky ray (default: 2)",
    )
    parser.add_argument(
        "--sky-pixel-chunk", type=int, default=2048,
        help="HEALPix rays integrated together (lower this to reduce memory)",
    )
    args = parser.parse_args()

    grid_size = tuple(args.grid_size)
    subgrids = tuple(args.subgrids)
    if any(grid_size[i] % subgrids[i] for i in range(3)):
        raise SystemExit("Grid dimensions must be divisible by subgrid dimensions")
    if not hp.isnsideok(args.nside):
        raise SystemExit("--nside must be a valid HEALPix NSIDE")
    prefix = args.output_prefix or args.snapshot.with_suffix("")
    prefix.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(args.snapshot, "r") as snapshot:
        box_size = np.asarray(snapshot["Header"].attrs["BoxSize"], dtype=float)
        time_myr = float(snapshot["Header"].attrs["Time"]) / MYR_IN_S
        density = read_grid(snapshot, "Density", subgrids, grid_size)
        number_density_cm3 = (
            read_grid(snapshot, "NumberDensity", subgrids, grid_size) / 1.0e6
        )
        neutral_fraction = read_grid(
            snapshot, "NeutralFractionH", subgrids, grid_size
        )
        temperature = read_grid(snapshot, "Temperature", subgrids, grid_size)

    np.clip(neutral_fraction, 0.0, 1.0, out=neutral_fraction)
    hplus_cm3 = number_density_cm3 * (1.0 - neutral_fraction)
    # Reuse the total-density array for H I to keep 256^3 memory use moderate.
    hi_cm3 = number_density_cm3
    hi_cm3 *= neutral_fraction
    del neutral_fraction
    # Hydrogen-only run: n_e = n_H+.  This is the requested relative Halpha
    # tracer, not a calibrated recombination-line luminosity.
    halpha = np.maximum(temperature, 1.0)
    np.power(halpha, -0.9, out=halpha)
    halpha *= hplus_cm3
    halpha *= hplus_cm3

    plot_cartesian(
        prefix, density, hplus_cm3, halpha, temperature, box_size, time_myr
    )
    del density
    halpha_sky, hi_sky = make_healpix_maps(
        hi_cm3, halpha, box_size, args.nside,
        radial_oversampling=args.sky_radial_oversampling,
        pixel_chunk=args.sky_pixel_chunk,
    )
    plot_sky(prefix, halpha_sky, hi_sky, args.nside, time_myr)
    plot_xrays(prefix, hplus_cm3, hplus_cm3, temperature, box_size, time_myr)
    print(
        f"Wrote {prefix}_cartesian.png, {prefix}_cartesian_linear.png, "
        f"{prefix}_sky.png and {prefix}_xray.png"
    )


if __name__ == "__main__":
    main()
