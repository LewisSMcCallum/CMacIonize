# Swiggum stellar-history shearing-box development run

This setup maps simulation time `0` to history time `-100 Myr` and runs to
`+20 Myr`. Stellar positions are linearly interpolated between trajectory rows.
Stars radiate only between their tabulated birth and death and while inside the
1 kpc box. A death crossed during a hydro step injects one supernova at the
tabulated death position when that position lies inside the box.

The supplied July 9 export contains only reconstructed ccSN progenitors. Its
own README states that direct-collapse outcomes are absent, so this run cannot
inject or record those missing events.

## Important rotation limitation

`GalacticShearingBox` initializes the solar-circle linear shear and advances
the local Coriolis and radial tidal source terms. The horizontal boundaries in
this development setup are ordinary periodic boundaries. They do **not** yet
perform the time-dependent azimuthal remap required by a complete radial
shearing-periodic boundary. This is suitable for developing and checking the
source history, feedback and local rotation terms, but not yet for a final
production shearing-box result over the full 120 Myr interval.

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
  --grid-size 128 128 128 --subgrids 8 8 8
```

The script plots the top-down gas surface density and uses `h5py`, `numpy`,
`matplotlib` and Pillow.
