#!/usr/bin/env python3

import ngl
import sys
import numpy
import scipy, scipy.sparse
import xarray
import dask
from .plot_map import plot_map
from ..e3sm_utils import get_data, area_average, regrid_data, get_coords, get_area_weights
from ..utils import nice_cntr_levels, apply_map
import cftime
import subprocess
import datetime


def convert_time(ds):

    # Convert cftime coordinate
    if isinstance(ds.time.values[0], cftime._cftime.DatetimeNoLeap):
        try:
            ds['time'] = ds.indexes['time'].to_datetimeindex()
        except:
            print('Could not convert times to datetimeindex. But proceeding anyway...')

    return ds


def subset_time(datasets):

    # Subset based on date range of first dataset
    t1 = max([ds.time[0].values for ds in datasets])
    t2 = min([ds.time[-1].values for ds in datasets])
    datasets = [ds.sel(time=slice(str(t1), str(t2))) for ds in datasets]
    return datasets


def compute_differences(da1, da2, x1, y1, x2, y2, map_file=None):

    # If coordinates are not the same, then we need to regrid
    if x1.shape != x2.shape or y1.shape != y2.shape:
        if map_file is None:
            print("                                                            ")
            print("Warning: no map file specified, will use griddata to regrid!")
            print("         This operation is SSSLLLLOOOOWWWWW for large grids!")
            print("         You should specify an offline map via the map_file ")
            print("         argument, or else refactor this code to calculate  ")
            print("         and cache an offline map itself.                   ")
            print("                                                            ")
            da1 = regrid_data(x1, y1, x2, y2, da1)
        else:
            # apply map
            with xarray.open_mfdataset(map_file) as ds_map:
                da1, *__ = apply_map(da1, ds_map, template=da2)

    # Now compute difference
    da_diff = da1 - da2
    da_diff.attrs = da1.attrs

    return da_diff


def compute_area_average(ds, vname):
    # Get weights; either use pre-computed or cosine(latitude) weights
    if 'area' in ds.variables.keys():
        wgt = ds.area
    else:
        wgt = numpy.cos(ds.lat * numpy.pi / 180.0)

    # Compute means
    return area_average(get_data(ds, vname), wgt)


def plot_panel(wks, data, x, y, **kwargs):

    # Make sure coordinate variables do not have a time dimension
    if 'time' in x.dims: x = x.isel(time=0)
    if 'time' in y.dims: y = y.isel(time=0)

    # Add cyclic point if we need to
    if len(data.shape) == 2 and x.max().values < 360:
        data_cyclic, x_cyclic = ngl.add_cyclic(data.values, x.values)
    else:
        data_cyclic = data.values
        x_cyclic = x.values

    # Select a colormap for this plot; if doing a difference plot that spans 0,
    # then choose a divergent colormap. Otherwise, stick with a perceptually
    # uniform colormap.
    if 'cnFillPalette' not in kwargs:
        if (data.min().values < 0 and data.max().values > 0):
            kwargs['cnFillPalette'] = 'MPL_coolwarm'
        else:
            kwargs['cnFillPalette'] = 'MPL_viridis'

    # Difference plots ought to have symmetric contour levels about 0
    if 'cnLevels' not in kwargs:
        if (data.min().values < 0 and data.max().values > 0):
            *__, levels = nice_cntr_levels(-abs(data).max().values,
                    abs(data).max().values, returnLevels=True, aboutZero=True)
            kwargs['cnLevelSelectionMode'] = 'ExplicitLevels'
            kwargs['cnLevels'] = levels
    else:
        # Make sure explicit levels is set if we passed cnLevels
        kwargs['cnLevelSelectionMode'] = 'ExplicitLevels'

    # Plot
    pl = plot_map(
        wks, x_cyclic, y.values, data_cyclic,
        lbOrientation='horizontal', 
        lbTitleString='%s (%s)'%(data.long_name, data.units),
        cnFillMode='RasterFill',
        cnLineLabelsOn=False, cnLinesOn=False, 
        nglDraw=False, nglFrame=False,
        **kwargs
    )

    return pl


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


def get_contour_levels(data_arrays, percentile=2, **kwargs):
    # Try to get robust contour intervals
    try:
        cmin = min([numpy.nanpercentile(da.values, percentile) for da in data_arrays])
        cmax = max([numpy.nanpercentile(da.values, 100-percentile) for da in data_arrays])
        if 'aboutZero' not in kwargs.keys(): kwargs['aboutZero'] = (cmin < 0 and cmax > 0)
        *__, clevels = nice_cntr_levels(cmin, cmax, returnLevels=True, max_steps=13, **kwargs)
    # Fall back to just using the min and max to set the limits
    except:
        cmin = min([da.min().values for da in data_arrays])
        cmax = max([da.max().values for da in data_arrays])
        if 'aboutZero' not in kwargs.keys(): kwargs['aboutZero'] = (cmin < 0 and cmax > 0)
        *__, clevels = nice_cntr_levels(cmin, cmax, returnLevels=True, max_steps=13, **kwargs)
    return clevels


def main(test_files, cntl_files, varname, plotname, 
        test_name="Model", cntl_name="Obs", time_offsets=[None, None], map_file=None, test_grid=None,
        cntl_grid=None, **kwargs):

    # Read data
    print('Open datasets...', end=''); sys.stdout.flush()
    datasets = [open_dataset(files, time_offset) for files, time_offset in zip((test_files,
        cntl_files), time_offsets)]
    grids = [xarray.open_mfdataset(f) if f is not None else None for f in (test_grid, cntl_grid)]
    print('done.')

    # Subset and time average
    print('Subset...', end=''); sys.stdout.flush()
    datasets = subset_time(datasets)
    datasets = [ds.mean(dim='time', keep_attrs=True) for ds in datasets]
    print('done.')

    # Select data
    print('Select data...', end=''); sys.stdout.flush()
    data_arrays = [get_data(ds, varname) for ds in datasets]
    coords = [get_coords(ds, ds_grid) for ds, ds_grid in zip(datasets, grids)]
    weights = [get_area_weights(ds) for ds in datasets]
    print('done.')

    # Compute differences
    print('Compute differences...', end=''); sys.stdout.flush()
    #datasets.append(compute_differences(datasets[0], datasets[1], varname, map_file))
    data_arrays.append(compute_differences(data_arrays[0], data_arrays[1], *coords[0], *coords[1], map_file))
    coords.append(coords[1])
    weights.append(weights[1])
    print('done.')

    # Compute global statistics for plot labels
    print('Compute statistics...', end=''); sys.stdout.flush()
    #means = [compute_area_average(ds, varname).values for ds in datasets]
    means = [area_average(d, w) for (d, w) in zip(data_arrays, weights)]
    print('done.')

    # Get common contour levels
    print('Find good clevels for plotting...', end=''); sys.stdout.flush()
    clevels = get_contour_levels([d for d in data_arrays[:-1]])
    dlevels = get_contour_levels([data_arrays[-1],], aboutZero=True)
    levels_list = [clevels, clevels, dlevels]
    print('done.')

    # Make plots
    print('Make plots...', end=''); sys.stdout.flush()
    plot_format = plotname.split('.')[-1]
    wks = ngl.open_wks(plot_format, plotname)
    diff_name = f'{test_name} minus {cntl_name}'
    labels = [f'{label} ({m.values:.1f})' for (label, m) in zip((test_name, cntl_name, diff_name), means)]
    plots = [
        plot_panel(wks, d, *c, tiMainString=label, cnLevels=levels, **kwargs)
        for (d, c, label, levels) in zip(data_arrays, coords, labels, levels_list)
    ]
    print('done.')

    # Panel plots
    print('Closing plot workspace and writing figures...', end=''); sys.stdout.flush()
    ngl.panel(wks, plots, [1, len(plots)])
    print('done.')

    # Finally, trim whitespace from our figures
    print('Trimming whitespace from figure...', end=''); sys.stdout.flush()
    subprocess.call(f'convert -trim {plotname} {plotname}'.split(' '))
    print('done.')

if __name__ == '__main__':
    import plac; plac.call(main)
