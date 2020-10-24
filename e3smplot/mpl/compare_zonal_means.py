#!/usr/bin/env python3
import xarray
import numpy
import dask
from glob import glob
from ..e3sm_utils import get_data, get_coords, get_area_weights, area_average, calculate_zonal_mean, get_grid_name, get_mapping_file, get_scrip_grid_ds, create_scrip_grid
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


def open_dataset(files, time_offset=None, **kwargs):
    # Open dataset as a dask array
    with dask.config.set(**{'array.slicing.split_large_chunks': True}):
        ds = xarray.open_mfdataset(
            sorted(files), combine='by_coords', drop_variables=('P3_output_dim', 'P3_input_dim'), **kwargs) 

    if time_offset is not None:
        print('Adding year offset...')
        ds['time'] = ds['time'] + time_offset

    # Rename coordinate variables
    if 'latitude' in ds.dims: ds = ds.rename({'latitude': 'lat'})
    if 'longitude' in ds.dims: ds = ds.rename({'longitude': 'lon'})

    if 'msshf' in ds.variables.keys():
        ds['msshf'].values = -ds['msshf'].values

    # Fix times so we can subset later
    ds = convert_time(ds)

    # Return dataset
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


def main(test_files, cntl_files, vname, fig_name, test_map=None, cntl_map=None,
        time_offsets=(None, None), test_name='Model', cntl_name='Obs', **kwargs):

    # Load datasets
    print('Load data...'); sys.stdout.flush()
    with dask.config.set(**{'array.slicing.split_large_chunks': True}):
        datasets = [open_dataset(f, time_offset=dt, chunks={'time': 1}) for (f, dt) in zip((test_files, cntl_files), time_offsets)]

    # Subset for overlapping time periods
    print('Get overlapping time range...'); sys.stdout.flush()
    t1 = max([ds.time[0].values for ds in datasets])
    t2 = min([ds.time[-1].values for ds in datasets])
    datasets = [ds.sel(time=slice(str(t1), str(t2))) for ds in datasets]

    # Select data
    print('Select data...'); sys.stdout.flush()
    data_arrays = [get_data(ds, vname) for ds in datasets]
    coords = [get_coords(ds) for ds in datasets]
    weights = [get_area_weights(ds) for ds in datasets]

    # Compute time averages
    print('Compute time averages...'); sys.stdout.flush()
    means = [da.mean(dim='time', keep_attrs=True) for da in data_arrays]
    weights = [wgt.mean(dim='time', keep_attrs=True) if 'time' in wgt.dims else wgt for wgt in weights]

    # Remap data to a lat/lon grid to compute zonal means. 
    print('Remap to lat/lon grid if needed...'); sys.stdout.flush()
    with dask.config.set(**{'array.slicing.split_large_chunks': True}):
        maps = [xarray.open_mfdataset(f) if f is not None else None for f in (test_map, cntl_map)]
    print('Remap weights...'); sys.stdout.flush()
    weights, lons, lats = zip(*[apply_map(m, f) if f is not None else m for (m, f) in zip(weights, maps)])
    print('Remap means...'); sys.stdout.flush()
    means, lons, lats = zip(*[apply_map(m, f) if f is not None else (m, m.lon, m.lat) for (m, f) in zip(means, maps)])
    print('Get lat...'); sys.stdout.flush()
    lats = [m.lat for m in means]
    print('done.')

    # Compute *zonal* average. Note that this is tricky for unstructured data.
    # Our approach is to compute means over latitude bins, and return a new
    # coordinate for latitude of the zonal mean. The function defined in
    # e3sm_utils takes the data, a set of latitude weights, and the input
    # latitudes, and optionally the number of new latitude bands within which to
    # average the data, and returns the binned/averaged data and the new
    # latitudes.
    #if map_file is not None:
    #    print('Apply map...')
    #    means[0] = apply_map(means[0], map_file, template=means[1])
    #    weights[0] = apply_map(weights[0], map_file, template=means[1])
    #    lats[0] = lats[1]
    #    weights, *__ = zip(*[xarray.broadcast(w, d) for (w, d) in zip(weights, means)])
    #    means = [zonal_mean(d, weights=w) for (d, w) in zip(means, weights)]
    #else:
    #    print('Try using our slow zonal mean routine...')
    #    means, lats = zip(*[calculate_zonal_mean(d, w, y) for (d, w, y) in zip(means, weights, lats)])
    print('Compute zonal means...', end=''); sys.stdout.flush()
    weights, *__ = zip(*[xarray.broadcast(w, d) for (w, d) in zip(weights, means)])
    means = [zonal_mean(d, weights=w) for (d, w) in zip(means, weights)]
    print('done.')

    # Make line plots of zonal averages
    print('Make line plots of zonal means...', end=''); sys.stdout.flush()
    figure, ax = pyplot.subplots(1, 1)
    plots = [plot_zonal_mean(y, m, label=l, **kwargs) for (y, m, l) in zip(lats, means, (test_name, cntl_name))]
    pyplot.legend()
    figure.savefig(fig_name, bbox_inches='tight')
    print('done.')

    # Finally, trim whitespace from our figures
    print('Trimming whitespace from figure...', end=''); sys.stdout.flush()
    subprocess.call(f'convert -trim {fig_name} {fig_name}'.split(' '))
    print('done.')


if __name__ == '__main__':
    import plac; plac.call(main)
