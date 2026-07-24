# Task-based cubed-sphere comparison

This is a deliberately small Cartesian/cubed-sphere comparison using the same
homogeneous gas and two separated off-centre sources.

From this directory:

```sh
mkdir -p output_cartesian output_spherical
CMacIonize --task-based --params cartesian.param --threads 4 --dirty --no-initial-output
CMacIonize --task-based --params spherical.param --threads 4 --dirty --no-initial-output
python3 plot_comparison.py
```

The coarse spherical grid is intentionally useful as a failure/control case:
16 angular cells per face under-resolves an off-centre source at about 5 pc.
The recommended comparison uses nearly equal cell counts:

```sh
CMacIonize --task-based --params cartesian_highres.param --threads 4 --dirty --no-initial-output
CMacIonize --task-based --params spherical_highres.param --threads 4 --dirty --no-initial-output
python3 plot_comparison.py --high-resolution
```

The Cartesian run has 262,144 cells and 0.25 pc cubic cells. The spherical run
has 294,912 cells, 0.166 pc radial spacing, and approximately 0.31 pc
tangential spacing at the distant source.

The spherical grid parameters are all in `spherical.param`. Set `radial
spacing` to `logarithmic` for log-R cells. `HDF5DensityFunction` can be selected
without a spherical-specific file format: its Cartesian input field is sampled
at each spherical cell midpoint.
