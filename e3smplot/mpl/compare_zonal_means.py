#!/usr/bin/env python3
import xarray
import numpy
import dask
from glob import glob
from ..e3sm_utils import get_data, get_area_weights, area_average, calculate_zonal_mean
from ..e3sm_utils import get_grid_name, get_mapping_file, get_scrip_grid_ds, create_scrip_grid, mask_all_zero
from ..e3sm_utils import open_dataset
from matplotlib import pyplot
import cftime
from ..utils import apply_map
import sys
import subprocess


def convert_time(ds):

    # Convert cftime coordinate
    if isinstance(ds.time.values[0], cftime._cftime.DatetimeNoLeap):
        try:
            ds['time'] = ds.indexes['time'].to_datetimeindex()
        except:
            print('Could not convert times to datetimeindex. But proceeding anyway...')

    return ds


def plot_zonal_mean(lat, data, **kwargs):
    ax = pyplot.gca()
    pl = ax.plot(lat, data, **kwargs)
    ax.set_xlabel('Latitude')
    ax.set_ylabel(f'{data.long_name} ({data.units})')
    return pl


def zonal_mean(data, weights=None):
    # Compute mean
    if weights is not None:
        weights, data = xarray.broadcast(weights, data)
        zonal_mean = (weights * data).sum(dim='lon', keep_attrs=True) / weights.sum(dim='lon', keep_attrs=True)
    else:
        zonal_mean = data.mean(dim='lon')

    # Copy attributes
    zonal_mean.attrs = data.attrs

    return zonal_mean


def to_latlon(da, x, y, map_file=None):
    if len(da.shape) == 1:
        if map_file is None:
            # Create latlon file
            print('Creating a default grid file...')

            # Get mapping file

        # Apply map
        da, x, y = apply_map(da)

        return da, x, y

def main(files, names, vname, fig_name, maps=None,
        time_offsets=None, t1=None, t2=None, verbose=False, dpi=800, **kwargs):

    if time_offsets is None: time_offsets = [None for x in files]
    if maps is None: maps = [None for x in files]

    # Load datasets
    if verbose: print('Load data...'); sys.stdout.flush()
    #with dask.config.set(**{'array.slicing.split_large_chunks': True}):
    datasets = [open_dataset(f, time_offset=dt, chunks={'time': 1}) for (f, dt) in zip(files, time_offsets)]

    # Subset for overlapping time periods, and compute time average
    if verbose: print('Get overlapping time range...'); sys.stdout.flush()
    if t1 is None: t1 = max([ds.time[0].values for ds in datasets])
    if t2 is None: t2 = min([ds.time[-1].values for ds in datasets])
    if verbose: print('Comparing period {} to {}'.format(str(t1), str(t2))); sys.stdout.flush()
    datasets = [ds.sel(time=slice(str(t1), str(t2))) for ds in datasets]
    if verbose: print('Compute time averages...'); sys.stdout.flush()
    datasets = [ds.mean(dim='time', keep_attrs=True) for ds in datasets]

    # Select (and remask) data
    if verbose: print('Select data...'); sys.stdout.flush()
    data_arrays = [mask_all_zero(get_data(ds, vname)) for ds in datasets]
    weights = [get_area_weights(ds) for ds in datasets]

    # Load data into memory
    if verbose: print('Force dask to execute and load into memory...'); sys.stdout.flush()
    for i in range(len(data_arrays)):
        data_arrays[i].load()
        weights[i].load()

    # Remap data to a lat/lon grid to compute zonal means. 
    if verbose: print('Remap to lat/lon grid if needed...'); sys.stdout.flush()
    maps = [xarray.open_dataset(f) if f is not None else None for f in maps]
    weights = [apply_map(m, f)[0] if f is not None else m for (m, f) in zip(weights, maps)]
    means = [apply_map(m, f)[0] if f is not None else m for (m, f) in zip(data_arrays, maps)]
    lats = [m.lat for m in means]

    # Compute *zonal* average. Note that this is tricky for unstructured data.
    # Our approach is to compute means over latitude bins, and return a new
    # coordinate for latitude of the zonal mean. The function defined in
    # e3sm_utils takes the data, a set of latitude weights, and the input
    # latitudes, and optionally the number of new latitude bands within which to
    # average the data, and returns the binned/averaged data and the new
    # latitudes.
    #if map_file is not None:
    #    if verbose: print('Apply map...')
    #    means[0] = apply_map(means[0], map_file, template=means[1])
    #    weights[0] = apply_map(weights[0], map_file, template=means[1])
    #    lats[0] = lats[1]
    #    weights, *__ = zip(*[xarray.broadcast(w, d) for (w, d) in zip(weights, means)])
    #    means = [zonal_mean(d, weights=w) for (d, w) in zip(means, weights)]
    #else:
    #    if verbose: print('Try using our slow zonal mean routine...')
    #    means, lats = zip(*[calculate_zonal_mean(d, w, y) for (d, w, y) in zip(means, weights, lats)])
    if verbose: print('Compute zonal means...'); sys.stdout.flush()
    weights, *__ = zip(*[xarray.broadcast(w, d) for (w, d) in zip(weights, means)])
    means = [zonal_mean(d, weights=w) for (d, w) in zip(means, weights)]
    #means = [d.mean(dim='lon', keep_attrs=True) for d in means]

    # Make line plots of zonal averages
    if verbose: print('Make line plots of zonal means...'); sys.stdout.flush()
    figure, ax = pyplot.subplots(1, 1)
    plots = [plot_zonal_mean(y, m, label=l, **kwargs) for (y, m, l) in zip(lats, means, names)]
    pyplot.legend()
    figure.savefig(fig_name, bbox_inches='tight', dpi=dpi)

    # Finally, trim whitespace from our figures
    if verbose: print('Trimming whitespace from figure...'); sys.stdout.flush()
    subprocess.call(f'convert -trim {fig_name} {fig_name}'.split(' '))


if __name__ == '__main__':
    import plac; plac.call(main)
