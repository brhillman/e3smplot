#!/usr/bin/env python3
import xarray
import numpy
import dask
from glob import glob
from ..e3sm_utils import get_data, get_coords, get_area_weights, area_average, calculate_zonal_mean, get_grid_name, get_mapping_file, get_scrip_grid_ds, create_scrip_grid
from pyngl import ngl
import cftime
from ..utils import apply_map
import sys


def open_dataset(files, **kwargs):
    # Open dataset as a dask array
    ds = xarray.open_mfdataset(sorted(files), combine='by_coords', drop_variables='P3_output_dim', **kwargs) 

    # Rename coordinate variables for consistency
    if 'latitude' in ds.dims: ds = ds.rename({'latitude': 'lat'})
    if 'longitude' in ds.dims: ds = ds.rename({'longitude': 'lon'})

    # Convert cftime DatetimeNoLeap to a datetimeindex
    if isinstance(ds.time.values[0], cftime._cftime.DatetimeNoLeap):
        try:
            ds['time'] = ds.indexes['time'].to_datetimeindex()
        except:
            print('Could not convert times to datetimeindex. But proceeding anyway...')

    # Return dataset
    return ds


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


def plot_lines(wks, lats, means, **kwargs):

    res = ngl.Resources()
    for key, val in kwargs.items():
        setattr(res, key, val)

    plot = ngl.xy(wks, lats, means, res)
    return plot


def main(test_files, cntl_files, vname, fig_name, test_map=None, cntl_map=None, test_name='Model', cntl_name='Obs', **kwargs):

    # Load datasets
    print('Load data...'); sys.stdout.flush()
    with dask.config.set(**{'array.slicing.split_large_chunks': True}):
        datasets = [open_dataset(f) for f in (test_files, cntl_files)]

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
    means = [apply_map(m, f) if f is not None else m for (m, f) in zip(means, (test_map, cntl_map))]
    weights = [apply_map(m, f) if f is not None else m for (m, f) in zip(weights, (test_map, cntl_map))]
    lats = [m.lat for m in means]

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

    # Make line plots of zonal averages
    print('Make line plots of zonal means...', end=''); sys.stdout.flush()
    plot_format = fig_name.split('.')[-1]
    wks = ngl.open_wks(plot_format, fig_name)
    pl = plot_lines(wks, lats, means, **kwargs)


if __name__ == '__main__':
    import plac; plac.call(main)
