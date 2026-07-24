#!/usr/bin/env python3
"""Compare Cartesian and cubed-sphere task-based ionization outputs."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import map_coordinates

PC = 3.085677581491367e16
ROOT = Path(__file__).resolve().parent
NCART = 32
NANG = 16
NRAD = 32
RMIN = 0.05
RMAX = 8.0


def read_output(pattern):
    filename = sorted(ROOT.glob(pattern))[-1]
    data = np.loadtxt(filename)
    return data[:, :3] / PC, data[:, 3] * (1.0 - data[:, 5])


def face_coordinates(points):
    axis = np.argmax(np.abs(points), axis=-1)
    sign = np.take_along_axis(points, axis[..., None], axis=-1)[..., 0] < 0
    face = 2 * axis + sign
    n = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0],
                  [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)
    au = np.array([[0, 1, 0], [0, -1, 0], [-1, 0, 0],
                   [1, 0, 0], [0, 1, 0], [0, 1, 0]], float)
    av = np.array([[0, 0, 1], [0, 0, 1], [0, 0, 1],
                   [0, 0, 1], [-1, 0, 0], [1, 0, 0]], float)
    denominator = np.sum(n[face] * points, axis=-1)
    u = np.sum(au[face] * points, axis=-1) / denominator
    v = np.sum(av[face] * points, axis=-1) / denominator
    return face, u, v


cart_xyz, cart_values = read_output("output_cartesian/cartesian_*.txt")
cart = np.empty((NCART, NCART, NCART))
cart_i = np.floor((cart_xyz + 6.0) * NCART / 12.0).astype(int)
cart[cart_i[:, 0], cart_i[:, 1], cart_i[:, 2]] = cart_values

sph_xyz, sph_values = read_output("output_spherical/spherical_*.txt")
sph = np.empty((6, NANG, NANG, NRAD))
sface, su, sv = face_coordinates(sph_xyz)
si = np.floor(0.5 * (su + 1.0) * NANG).astype(int)
sj = np.floor(0.5 * (sv + 1.0) * NANG).astype(int)
sr = np.linalg.norm(sph_xyz, axis=1)
sk = np.floor((sr - RMIN) * NRAD / (RMAX - RMIN)).astype(int)
sph[sface, np.clip(si, 0, NANG - 1), np.clip(sj, 0, NANG - 1),
    np.clip(sk, 0, NRAD - 1)] = sph_values


def sample_cartesian(points):
    shape = points.shape[:-1]
    coordinates = (
        (points.reshape(-1, 3) + 6.0) * NCART / 12.0 - 0.5).T
    values = map_coordinates(cart, coordinates, order=1, mode="constant",
                            cval=0.0)
    return values.reshape(shape)


def sample_spherical(points):
    shape = points.shape[:-1]
    flat = points.reshape(-1, 3)
    radius = np.linalg.norm(flat, axis=1)
    face, u, v = face_coordinates(flat)
    result = np.zeros(len(flat))
    for this_face in range(6):
        select = face == this_face
        coordinates = np.vstack((
            0.5 * (u[select] + 1.0) * NANG - 0.5,
            0.5 * (v[select] + 1.0) * NANG - 0.5,
            (radius[select] - RMIN) * NRAD / (RMAX - RMIN) - 0.5,
        ))
        result[select] = map_coordinates(
            sph[this_face], coordinates, order=1, mode="nearest")
    result[(radius < RMIN) | (radius > RMAX)] = 0.0
    return result.reshape(shape)


# Integrate both solutions on exactly the same Cartesian sampling volume.
nplot = 128
edges = np.linspace(-6.0, 6.0, nplot + 1)
centres = 0.5 * (edges[1:] + edges[:-1])
x, y, z = np.meshgrid(centres, centres, centres, indexing="ij")
points = np.stack((x, y, z), axis=-1)
common = np.linalg.norm(points, axis=-1) <= RMAX
cart_column = np.sum(sample_cartesian(points) * common, axis=2) * (
    12.0 * PC / nplot)
sph_column = np.sum(sample_spherical(points) * common, axis=2) * (
    12.0 * PC / nplot)

# All-sky columns from the origin, using common rays and radial samples.
nlon, nlat, nray = 240, 120, 160
lon = np.linspace(-np.pi, np.pi, nlon)
lat = np.linspace(-0.5 * np.pi, 0.5 * np.pi, nlat)
llon, llat = np.meshgrid(lon, lat)
directions = np.stack((np.cos(llat) * np.cos(llon),
                       np.cos(llat) * np.sin(llon), np.sin(llat)), axis=-1)
radii = np.linspace(RMIN, 6.0, nray)
ray_points = directions[..., None, :] * radii[None, None, :, None]
cart_sky = np.trapezoid(sample_cartesian(ray_points), radii * PC, axis=2)
sph_sky = np.trapezoid(sample_spherical(ray_points), radii * PC, axis=2)

positive = np.concatenate((cart_column[cart_column > 0],
                           sph_column[sph_column > 0]))
vmin, vmax = np.percentile(positive, [2, 99.5])
fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
for ax, image, title in zip(
        axes[0], (cart_column, sph_column),
        ("Cartesian: integrated HII density",
         "Cubed sphere: integrated HII density")):
    plotted = ax.imshow(image.T, origin="lower", extent=(-6, 6, -6, 6),
                        norm="log", vmin=vmin, vmax=vmax, cmap="magma")
    ax.plot([0.25, 2.5], [0.0, 1.0], "c+", ms=9, mew=1.5)
    ax.set(title=title, xlabel="x (pc)", ylabel="y (pc)")
fig.colorbar(plotted, ax=axes[0], label=r"$N_{\rm HII}$ (m$^{-2}$)")

sky_positive = np.concatenate((cart_sky[cart_sky > 0],
                               sph_sky[sph_sky > 0]))
svmin, svmax = np.percentile(sky_positive, [2, 99.5])
for ax, image, title in zip(
        axes[1], (cart_sky, sph_sky),
        ("Cartesian viewed from (0,0,0)",
         "Cubed sphere viewed from (0,0,0)")):
    plotted_sky = ax.imshow(
        image, origin="lower", extent=(-180, 180, -90, 90), aspect="auto",
        norm="log", vmin=svmin, vmax=svmax, cmap="viridis")
    ax.set(title=title, xlabel="longitude (deg)", ylabel="latitude (deg)")
fig.colorbar(plotted_sky, ax=axes[1], label=r"$N_{\rm HII}$ (m$^{-2}$)")

output = ROOT / "two_source_comparison.png"
fig.savefig(output, dpi=180)
print(output)
