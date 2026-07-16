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
from matplotlib.colors import LogNorm
import numpy as np


PC_IN_M = 3.0856775814913673e16
KPC_IN_M = 1000.0 * PC_IN_M
M_IN_CM = 100.0
MYR_IN_S = 3.15576e13
KEV_IN_J = 1.602176634e-16
BOLTZMANN = 1.380649e-23

# Approximate eROSITA bands in keV.  The emission model below includes only
# hydrogen free-free continuum, so soft-band metal-line emission is absent.
XRAY_BANDS_KEV = ((0.2, 0.6), (0.6, 2.3), (2.3, 5.0))


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


def positive_limits(values, lower=1.0, upper=99.5):
    finite = np.asarray(values)[np.isfinite(values) & (np.asarray(values) > 0.0)]
    if finite.size == 0:
        return 1.0e-30, 1.0
    low, high = np.percentile(finite, (lower, upper))
    if high <= low:
        high = low * 10.0
    return low, high


def show_log(axis, values, extent, title, label, cmap="magma"):
    low, high = positive_limits(values)
    image = axis.imshow(
        np.asarray(values).T,
        origin="lower",
        extent=extent,
        cmap=cmap,
        norm=LogNorm(low, high),
        interpolation="nearest",
    )
    axis.set_title(title)
    axis.set_xlabel("X toward Galactic centre [kpc]")
    axis.set_ylabel("Y along Galactic rotation [kpc]")
    plt.colorbar(image, ax=axis, label=label)


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


def make_healpix_maps(hi_cm3, halpha_proxy, box_size, nside):
    """Deposit cell-volume/r^2 contributions into centre-observer sky pixels."""
    grid_size = np.asarray(hi_cm3.shape)
    cell_size_m = box_size / grid_size
    centres = [
        (np.arange(grid_size[i]) + 0.5) * cell_size_m[i] - 0.5 * box_size[i]
        for i in range(3)
    ]
    x, y = np.meshgrid(centres[0], centres[1], indexing="ij")
    npix = hp.nside2npix(nside)
    solid_angle = hp.nside2pixarea(nside)
    hi_map = np.zeros(npix)
    halpha_map = np.zeros(npix)
    cell_volume_cm3 = np.prod(cell_size_m) * M_IN_CM ** 3

    # Work one z plane at a time to avoid allocating full 3-D position arrays.
    for iz, z in enumerate(centres[2]):
        radius2_m2 = x * x + y * y + z * z
        valid = radius2_m2 > 0.0
        pixels = hp.vec2pix(nside, x[valid], y[valid], np.full(valid.sum(), z))
        geometric_weight = (
            cell_volume_cm3 / (radius2_m2[valid] * M_IN_CM ** 2 * solid_angle)
        )
        np.add.at(hi_map, pixels, hi_cm3[:, :, iz][valid] * geometric_weight)
        np.add.at(
            halpha_map,
            pixels,
            halpha_proxy[:, :, iz][valid] * geometric_weight,
        )
    return halpha_map, hi_map


def plot_cartesian(prefix, density, hplus_cm3, halpha, temperature, box_size,
                   time_myr):
    dz_m = box_size[2] / density.shape[2]
    dz_cm = dz_m * M_IN_CM
    dz_pc = dz_m / PC_IN_M
    extent = (-0.5 * box_size[0] / KPC_IN_M, 0.5 * box_size[0] / KPC_IN_M,
              -0.5 * box_size[1] / KPC_IN_M, 0.5 * box_size[1] / KPC_IN_M)
    midplane = density.shape[2] // 2
    figure, axes = plt.subplots(2, 4, figsize=(20, 10), constrained_layout=True)

    show_log(axes[0, 0], density.sum(axis=2) * dz_m, extent,
             "Gas surface density", r"$\Sigma$ [kg m$^{-2}$]")
    show_log(axes[0, 1], hplus_cm3.sum(axis=2) * dz_cm, extent,
             "Ionized-hydrogen column", r"$N_{H^+}$ [cm$^{-2}$]", "viridis")
    show_log(axes[0, 2], halpha.sum(axis=2) * dz_pc, extent,
             r"H$\alpha$ proxy projection",
             r"$\int n_e^2 T^{-0.9} dz$ [cm$^{-6}$ pc K$^{-0.9}$]", "inferno")
    axes[0, 3].axis("off")

    show_log(axes[1, 0], density[:, :, midplane] * 1.0e-3, extent,
             "Midplane gas density", r"$\rho$ [g cm$^{-3}$]")
    show_log(axes[1, 1], hplus_cm3[:, :, midplane], extent,
             "Midplane ionized-hydrogen density", r"$n_{H^+}$ [cm$^{-3}$]",
             "viridis")
    show_log(axes[1, 2], halpha[:, :, midplane], extent,
             r"Midplane H$\alpha$ proxy", r"$n_e^2 T^{-0.9}$", "inferno")
    show_log(axes[1, 3], temperature[:, :, midplane], extent,
             "Midplane temperature", "T [K]", "plasma")
    figure.suptitle(f"Simulation time {time_myr:.3f} Myr")
    figure.savefig(f"{prefix}_cartesian.png", dpi=160)
    plt.close(figure)


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
            coord=None, notext=False,
        )
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
        show_log(
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
        hi_cm3, halpha, box_size, args.nside
    )
    plot_sky(prefix, halpha_sky, hi_sky, args.nside, time_myr)
    plot_xrays(prefix, hplus_cm3, hplus_cm3, temperature, box_size, time_myr)
    print(f"Wrote {prefix}_cartesian.png, {prefix}_sky.png and {prefix}_xray.png")


if __name__ == "__main__":
    main()
