#!/usr/bin/env python3

import ngl
import sys
import numpy
import scipy, scipy.sparse
import xarray
from .plot_map import plot_map
from ..e3sm_utils import get_data, area_average, regrid_data, get_coords
from ..utils import nice_cntr_levels, apply_map
import cftime


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
    t1, t2 = datasets[0].time[0].values, datasets[0].time[-1].values
    for i in range(1, len(datasets)):
        datasets[i] = datasets[i].sel(time=slice(str(t1), str(t2)))
    return datasets


#def regrid(ds1, ds2, vname, method='nco'):
def compute_differences(ds1, ds2, vname, map_file=None):

    # Get data
    da1 = get_data(ds1, vname)
    da2 = get_data(ds2, vname)
    x1, y1 = get_coords(ds1)
    x2, y2 = get_coords(ds2) 

    # If coordinates are not the same, then we need to regrid
    if ds1.lon.shape != ds2.lon.shape or ds1.lat.shape != ds2.lat.shape:
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
            if True:
                da1 = apply_map(da1, map_file, template=da2)
            else:
                ds_map = xarray.open_dataset(map_file)
                weights = scipy.sparse.coo_matrix((ds_map['S'].values, (ds_map['row'].values-1, ds_map['col'].values-1)))
                da1_regrid = weights.dot(da1)
                da1_regrid = xarray.DataArray(
                    da1_regrid.reshape(da2.shape), dims=da2.dims, coords=da2.coords,
                    attrs=da1.attrs
                )
                da1 = da1_regrid

    # Now compute difference
    da_diff = da1 - da2
    da_diff.attrs = da1.attrs

    # Return new dataset
    ds_diff = xarray.Dataset({vname: da_diff, 'lon': x2, 'lat': y2})

    return ds_diff


def compute_area_average(ds, vname):
    # Get weights; either use pre-computed or cosine(latitude) weights
    if 'area' in ds.variables.keys():
        wgt = ds.area
    else:
        wgt = numpy.cos(ds.lat * numpy.pi / 180.0)

    # Compute means
    return area_average(get_data(ds, vname), wgt)


def plot_panel(wks, ds, vname, ds_grid=None, **kwargs):

    # Select data for plot
    data = get_data(ds, vname)
    x, y = get_coords(ds, ds_grid)

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


def open_dataset(files, **kwargs):
    # Open dataset as a dask array
    ds = xarray.open_mfdataset(
        sorted(files), combine='by_coords', drop_variables=('P3_output_dim', 'P3_input_dim'), **kwargs) 

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
        test_name="Model", cntl_name="Obs", map_file=None, **kwargs):

    # Read data
    print('Open datasets...', end=''); sys.stdout.flush()
    datasets = [open_dataset(files) for files in (test_files, cntl_files)]
    print('done.')

    # Subset and time average
    print('Subset...', end=''); sys.stdout.flush()
    datasets = subset_time(datasets)
    datasets = [ds.mean(dim='time', keep_attrs=True) for ds in datasets]
    print('done.')

    # Compute differences
    print('Compute differences...', end=''); sys.stdout.flush()
    datasets.append(compute_differences(datasets[0], datasets[1], varname, map_file))
    print('done.')

    # Compute global statistics for plot labels
    print('Compute statistics...', end=''); sys.stdout.flush()
    means = [compute_area_average(ds, varname).values for ds in datasets]
    print('done.')

    # Get common contour levels
    print('Find good clevels for plotting...', end=''); sys.stdout.flush()
    clevels = get_contour_levels([get_data(ds, varname) for ds in datasets[:-1]])
    dlevels = get_contour_levels([get_data(datasets[-1], varname),], aboutZero=True)
    levels_list = [clevels, clevels, dlevels]
    print('done.')

    # Make plots
    print('Make plots...', end=''); sys.stdout.flush()
    plot_format = plotname.split('.')[-1]
    wks = ngl.open_wks(plot_format, plotname)
    diff_name = f'{test_name} minus {cntl_name}'
    labels = [f'{label} ({m:.1f})' for (label, m) in zip((test_name, cntl_name, diff_name), means)]
    plots = [
        plot_panel(wks, ds, varname, tiMainString=label, cnLevels=levels, **kwargs)
        for (ds, label, levels) in zip(datasets, labels, levels_list)
    ]
    print('done.')

    # Panel plots
    print('Closing plot workspace and writing figures...', end=''); sys.stdout.flush()
    ngl.panel(wks, plots, [1, len(plots)])
    print('done.')


if __name__ == '__main__':
    import plac; plac.call(main)
