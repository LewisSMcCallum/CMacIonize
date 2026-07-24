# Task-based cubed-sphere comparison

This is a deliberately small Cartesian/cubed-sphere comparison using the same
homogeneous gas and two off-centre sources.

From this directory:

```sh
mkdir -p output_cartesian output_spherical
CMacIonize --task-based --params cartesian.param --threads 4 --dirty --no-initial-output
CMacIonize --task-based --params spherical.param --threads 4 --dirty --no-initial-output
python3 plot_comparison.py
```

The spherical grid parameters are all in `spherical.param`. Set `radial
spacing` to `logarithmic` for log-R cells. `HDF5DensityFunction` can be selected
without a spherical-specific file format: its Cartesian input field is sampled
at each spherical cell midpoint.
