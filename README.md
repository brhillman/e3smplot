# e3smplot
Collection of scripts and notebooks to perform quicklook analysis of E3SM output

## Prequisites
The following external libraries are needed:
    - xarray
    - netcdf4
    - dask
    - numpy
    - scipy
    - pyngl
    - matplotlib
    - plac

The easiest way to install is to create a conda environment:
```
    conda create --name=e3smplot -c conda-forge python=3 xarray netcdf4 dask numpy scipy pyngl matplotlib plac
    conda activate e3smplot
```

## Installing
You can install e3smplot into your environment by running the `setup.py` script which uses `distutils`. To be able to continue to receive updates via `git pulls`, you may want to install in `develop` mode:
```
    python setup.py develop
```
You should then be able to run the analysis scripts, and also import modules from the library from any directory on your system.

## Running
e3smplot includes standalone tools to make individual plot types, or you can import the codes directly into a script. `scripts/make_analysis.py` provides an example of how to generate many plots at once and compare against observations. You will want to edit the paths and test and cntl names for your local environment. You can also run each of the plotting scripts individually from the command line (thanks to the command line interfaces provided by the `plac` library). See code in `e3smplot/pyngl` and `e3smplot/mpl` (proper documentation on this is a TODO item).
