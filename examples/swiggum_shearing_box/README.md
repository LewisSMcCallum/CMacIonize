# Swiggum stellar-history shearing-box development run

This setup maps simulation time `0` to history time `-100 Myr` and runs to
`+20 Myr`. Stellar positions are linearly interpolated between trajectory rows.
Stars radiate only between their tabulated birth and death and while inside the
1 kpc box. A death crossed during a hydro step injects one supernova at the
tabulated death position when that position lies inside the box.

The supplied July 9 export contains only reconstructed ccSN progenitors. Its
own README states that direct-collapse outcomes are absent, so this run cannot
inject or record those missing events.

## Galactic rotation boundaries

`GalacticShearingBox` initializes the solar-circle linear shear and advances
the local Coriolis and radial tidal source terms. With
`shearing periodic boundaries: true`, the hydrodynamic radial faces use a
time-dependent, linearly interpolated y remap. Boundary fluxes are deposited
conservatively, with y momentum and total energy transformed between the two
background-shear frames. This removes the artificial edge-to-edge velocity
jump produced by an ordinary x-periodic boundary.

The photon transport still uses the ordinary geometric periodic neighbour
links. The new remap therefore fixes the hydrodynamic boundary discontinuity;
it does not yet remap a packet's y position when that packet crosses a radial
boundary.

## Build and run

The new source file and the hydrogen-only configuration require rerunning CMake
once:

```bash
cd build
cmake -DHYDROGEN_ONLY=ON ..
make -j
cd rundir
./CMacIonize --task-based-rhd \
  --params ../../examples/swiggum_shearing_box/swiggum_shearing_box.param
```

The parameter file expects the compressed history at
`build/data/sn_history_july9_imf_sample0_minus100_to_plus20myr.csv.gz`. CMake
copies the repository version there during configuration.

The default 128-cubed, one-million-packet setup is a development resolution.
Use short `--number-of-steps` runs first to establish memory use and runtime on
the target cluster before launching 120 Myr.

## Density movie

From the directory containing the snapshots:

```bash
python ../../examples/swiggum_shearing_box/plot_density_gif.py \
  'swiggum_*.hdf5' swiggum_density.gif \
  --grid-size 128 128 128 --subgrids 8 8 8 --track-stars
```

The script plots the top-down gas surface density and uses `h5py`, `numpy`,
`matplotlib` and Pillow. With `--track-stars`, crosses show living ionizing
stars with size set by Q_H0; circles show supernova remnants expanding at
10 km/s and fading for 1 Myr. Projected grids are cached in the ignored
`.swiggum_plot_cache` directory, so subsequent runs only read new or changed
snapshots. Omit `--track-stars` for a density-only animation. The constants
`SN_CIRCLE_SPEED_KMS` and `SN_CIRCLE_FADE_TIME_MYR` near the top of the script
control the displayed circle expansion and lifetime.

For detailed diagnostics of one snapshot (including HEALPix maps, requiring
`healpy`):

```bash
python ../../examples/swiggum_shearing_box/plot_snapshot_diagnostics.py \
  swiggum_0020.hdf5 --grid-size 256 256 256 --subgrids 8 8 8 \
  --nside 64
```

This writes `_cartesian.png`, `_sky.png`, and `_xray.png` figures beside the
snapshot. The Halpha panels use the requested relative tracer
`n_e^2 T^-0.9`. The three X-ray panels are intrinsic hydrogen free-free
estimates for 0.2–0.6, 0.6–2.3, and 2.3–5.0 keV; they omit metal lines and
foreground absorption and should not be interpreted as calibrated eROSITA
count-rate predictions.
