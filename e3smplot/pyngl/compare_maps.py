#!/usr/bin/env python3

import ngl
import sys
import numpy
import scipy, scipy.sparse
import xarray
import dask
from .plot_map import plot_map
from ..e3sm_utils import get_data, area_average, regrid_data, get_coords, get_area_weights, mask_all_zero
from ..e3sm_utils import get_mapping_file
from ..e3sm_utils import open_dataset
from ..utils import nice_cntr_levels, apply_map
import cftime
import subprocess
import datetime


def common_time_limits(datasets):
    t1 = max([ds.time[0].values for ds in datasets])
    t2 = min([ds.time[-1].values for ds in datasets])
    return t1, t2

def compute_differences(da1, da2, x1, y1, x2, y2, map_file=None):

    # If coordinates are not the same, then we need to regrid
    if x1.shape != x2.shape or y1.shape != y2.shape:
        if map_file is None:
            map_file = get_mapping_file(da1, da2)
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
    print('Open datasets...'); sys.stdout.flush()
    datasets = [open_dataset(files, time_offset) for files, time_offset in zip((test_files,
        cntl_files), time_offsets)]
    grids = [xarray.open_mfdataset(f) if f is not None else None for f in (test_grid, cntl_grid)]

    # Subset data
    print('Subset consistent time periods...'); sys.stdout.flush()
    t1, t2 = common_time_limits(datasets)
    datasets = [ds.sel(time=slice(str(t1), str(t2))) for ds in datasets]

    # Select data
    print('Select data...'); sys.stdout.flush()
    data_arrays = [mask_all_zero(get_data(ds, varname)) for ds in datasets]
    coords = [get_coords(ds, ds_grid) for ds, ds_grid in zip(datasets, grids)]
    weights = [get_area_weights(ds) for ds in datasets]

    # Time average
    print('Compute time averages...'); sys.stdout.flush()
    data_arrays = [da.mean(dim='time', keep_attrs=True) for da in data_arrays]
    coords = [[c.mean(dim='time', keep_attrs=True) for c in coord] if 'time' in
            coord[0].dims else coord for coord in coords]
    weights = [da.mean(dim='time', keep_attrs=True)  if 'time'
            in da.dims else da for da in weights]

    # Compute differences
    print('Compute differences...'); sys.stdout.flush()
    data_arrays.append(compute_differences(data_arrays[0], data_arrays[1], *coords[0], *coords[1], map_file))
    coords.append(coords[1])
    weights.append(weights[1])

    # Compute global statistics for plot labels
    print('Compute statistics...'); sys.stdout.flush()
    means = [area_average(d, w) for (d, w) in zip(data_arrays, weights)]

    # Get common contour levels
    print('Find good clevels for plotting...'); sys.stdout.flush()
    clevels = get_contour_levels([d for d in data_arrays[:-1]])
    dlevels = get_contour_levels([data_arrays[-1],], aboutZero=True)
    levels_list = [clevels, clevels, dlevels]

    # Make plots
    print('Make plots...'); sys.stdout.flush()
    plot_format = plotname.split('.')[-1]
    wks = ngl.open_wks(plot_format, plotname)
    diff_name = f'{test_name} minus {cntl_name}'
    labels = [f'{label} ({m.values:.1f})' for (label, m) in zip((test_name, cntl_name, diff_name), means)]
    plots = [
        plot_panel(wks, d, *c, tiMainString=label, cnLevels=levels, **kwargs)
        for (d, c, label, levels) in zip(data_arrays, coords, labels, levels_list)
    ]

    # Panel plots
    print('Closing plot workspace and writing figures...'); sys.stdout.flush()
    ngl.panel(wks, plots, [1, len(plots)])
    ngl.destroy(wks)

    # Finally, trim whitespace from our figures
    print('Trimming whitespace from figure...'); sys.stdout.flush()
    subprocess.call(f'convert -trim {plotname} {plotname}'.split(' '))

if __name__ == '__main__':
    import plac; plac.call(main)
