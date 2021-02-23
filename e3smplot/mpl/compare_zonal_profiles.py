#!/usr/bin/env python3
import xarray
import numpy
import dask
from glob import glob
from e3smplot.e3sm_utils import get_data, get_coords, get_area_weights, area_average, calculate_zonal_mean
from e3smplot.e3sm_utils import get_grid_name, get_mapping_file, get_scrip_grid_ds, create_scrip_grid, mask_all_zero
from e3smplot.e3sm_utils import open_dataset
from matplotlib import pyplot
import cftime
from e3smplot.utils import interp_along_axis
from e3smplot.smm import apply_weights_wrap
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


def plot_zonal_profile(d, label, **kwargs):
    ax = pyplot.gca()
    pl = ax.pcolormesh(
        d.lat, d.plev, d.transpose('plev', 'lat'),
        **kwargs
    )
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.set_xlabel('Latitude (degrees north)')
    ax.set_ylabel('Pressure (hPa)')
    return pl
    

def make_zonal_profile_figure(means, cases, **kwargs):

    # Find good contour levels
    vmin = min([d.min().values for d in means])
    vmax = max([d.max().values for d in means])

    # Plot zonal profiles
    from matplotlib import pyplot
    figure, axes = pyplot.subplots(3, 1, figsize=(10, 15))
    for icase, (case, d) in enumerate(zip(cases, means)):
        ax = figure.add_axes(axes[icase])
        pl = plot_zonal_profile(d, case, vmin=vmin, vmax=vmax, shading='nearest')
        cb = pyplot.colorbar(
            pl, orientation='horizontal', label=f'{d.long_name} ({d.units})',
            shrink=0.8,
        )
        ax.set_title(case)

    # Plot difference
    ax = figure.add_axes(axes[-1])
    d = means[0] - means[1]
    d.attrs = means[0].attrs
    vmax =  abs(d).max().values
    vmin = -abs(d).max().values
    pl = plot_zonal_profile(d, case, vmin=vmin, vmax=vmax, shading='nearest', cmap='RdBu_r')
    cb = pyplot.colorbar(
        pl, orientation='horizontal', label=f'{d.long_name} ({d.units})',
        shrink=0.8,
    )
    ax.set_title(f'{cases[0]} - {cases[1]}')

    return figure


def zonal_mean(data, weights=None, xname='lon'):
    # Compute mean
    if weights is not None:
        weights, data = xarray.broadcast(weights, data)
        zonal_mean = (weights * data).sum(dim=xname, keep_attrs=True) / weights.sum(dim=xname, keep_attrs=True)
    else:
        zonal_mean = data.mean(dim=xname)

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

def get_height(data):
    if 'lev' in data.dims:
        height = data.lev
    if 'plev' in data.dims:
        height = data.plev
    elif 'level' in data.dims:
        height = data.level
    else:
        raise RuntimeError('No valid height coordinate.')
    return height

def main(test_files, cntl_files, names, vname, fig_name, maps=None,
        time_offsets=None, t1=None, t2=None, dpi=400, verbose=False, **kwargs):

    if time_offsets is None: time_offsets = [None for x in files]
    if maps is None: maps = [None for x in files]

    # Load datasets
    print('Load data...'); sys.stdout.flush()
    #with dask.config.set(**{'array.slicing.split_large_chunks': True}):
    datasets = [open_dataset(f, time_offset=dt, chunks={'time': 1}) for (f, dt) in zip((test_files, cntl_files), time_offsets)]

    # Subset for overlapping time periods
    if verbose: print('Get overlapping time range...'); sys.stdout.flush()
    if t1 is None: t1 = max([ds.time[0].values for ds in datasets])
    if t2 is None: t2 = min([ds.time[-1].values for ds in datasets])
    datasets = [ds.sel(time=slice(str(t1), str(t2))) for ds in datasets]

    # Rename
    datasets = [ds.rename({'level': 'plev'}) if 'level' in ds.dims else ds for ds in datasets]

    # Select (and remask) data
    if verbose: print('Select data...'); sys.stdout.flush()
    data_arrays = [mask_all_zero(get_data(ds, vname)) for ds in datasets]
    coords = [get_coords(ds) for ds in datasets]
    weights = [get_area_weights(ds) for ds in datasets]

    for ii in range(len(weights)):
        *__, weights[ii] = xarray.broadcast(data_arrays[ii], weights[ii])

    # TODO: we need to interpolate to common pressure levels here

    # Compute time averages
    if verbose: print('Compute time averages...'); sys.stdout.flush()
    means = [da.mean(dim='time', keep_attrs=True) for da in data_arrays]
    weights = [wgt.mean(dim='time', keep_attrs=True) if 'time' in wgt.dims else wgt for wgt in weights]

    # Remap data to a common grid (since we will both compute zonal mean and
    # compute diffs)
    if verbose: print('Remap to common grid...'); sys.stdout.flush()
    dims = means[-1].dims[-2:]
    coords = {d: means[-1].coords[d] for d in dims}
    means = [
        apply_weights_wrap(f, m, x=coords['lon'], y=coords['lat']) 
        if f is not None else m for (f, m) in zip(maps, means)
    ]
    weights = [
        apply_weights_wrap(f, m, x=coords['lon'], y=coords['lat']) 
        if f is not None else m for (f, m) in zip(maps, weights)
    ]

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
    if verbose: print('Compute zonal means...'); sys.stdout.flush()
    weights, *__ = zip(*[xarray.broadcast(w, d) for (w, d) in zip(weights, means)])
    means = [zonal_mean(d, weights=w) for (d, w) in zip(means, weights)]

    # Make plots of zonal averages
    if verbose: print('Make pcolor plots of zonal means...'); sys.stdout.flush()
    figure = make_zonal_profile_figure(means, names)
    figure.savefig(fig_name, bbox_inches='tight', dpi=dpi)

    # Finally, trim whitespace from our figures
    if verbose: print('Trimming whitespace from figure...'); sys.stdout.flush()
    subprocess.call(f'convert -trim {fig_name} {fig_name}'.split(' '))


if __name__ == '__main__':
    import plac; plac.call(main)
