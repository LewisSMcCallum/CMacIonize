# Edenhofer cubed-sphere input

`convert_edenhofer_to_cubed_sphere.py` converts the mean Edenhofer et al.
(2023) HEALPix map to the cell ordering used by CMacIonize's task-based
cubed-sphere grid. It stores hydrogen number density in SI units. The
user-supplied conversion is applied exactly:

```
n_H [m^-3] = Edenhofer differential-extinction density * 2700 [cm^-3] * 10^6
```

The first radial cell contains the integrated inner map divided by its outer
radius, so its radial column is preserved. The other 516 cells retain the
native Edenhofer radial boundaries. Angular interpolation follows the
`dustmaps` convention and is bilinear in log density.

The standard outputs have `N=64`, `128`, and `256` cells along each edge of
each cubed-sphere face. `N=256` is the conservative maximum: it does not claim
finer linear face sampling than the native `NSIDE=256` data.

Example:

```
python convert_edenhofer_to_cubed_sphere.py \
  /path/to/mean_and_std_healpix.fits /path/to/output
python plot_edenhofer_sky.py /path/to/output/edenhofer_cubed_n128_r517.hdf5
```

CMacIonize must be compiled with HDF5. Use `DensityFunction:type:
SphericalHDF5`, and set `SphericalDensityGrid:radial spacing: file` plus
`SphericalDensityGrid:radial edges filename` to the same HDF5 file.

Use the Gadget density-grid writer for simulation output (its default fields
already include coordinates, number density and neutral hydrogen fraction).
The result can then be plotted without loading the whole snapshot into memory:

```
python plot_ionization_sky.py snapshot_010.hdf5 \
  edenhofer_cubed_n128_r517.hdf5 --component ionized
```
