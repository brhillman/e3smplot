#!/usr/bin/env python3
import xarray
import numpy
from glob import glob
from ..e3sm_utils import get_data, get_coords, get_area_weights, area_average, calculate_zonal_mean
from matplotlib import pyplot
import cftime
from ..utils import apply_map


def open_dataset(files, **kwargs):
    # Open dataset as a dask array
    ds = xarray.open_mfdataset(sorted(files), combine='by_coords', drop_variables='P3_output_dim', **kwargs) 

    # Rename coordinate variables
    if 'latitude' in ds.dims: ds = ds.rename({'latitude': 'lat'})
    if 'longitude' in ds.dims: ds = ds.rename({'longitude': 'lon'})

    if 'msshf' in ds.variables.keys():
        ds['msshf'].values = -ds['msshf'].values

    # Convert cftime coordinate
    if isinstance(ds.time.values[0], cftime._cftime.DatetimeNoLeap):
        try:
            ds['time'] = ds.indexes['time'].to_datetimeindex()
        except:
            print('Could not convert times to datetimeindex. But proceeding anyway...')

    # Return dataset
    return ds


def plot_zonal_mean(lat, data, **kwargs):
    ax = pyplot.gca()
    pl = ax.plot(lat, data, **kwargs)
    ax.set_xlabel(f'{data.long_name} ({data.units})')
    ax.set_ylabel('Latitude')
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


def to_latlon(data, x, y):

    if len(data.shape) == 1:
        # Data is unstructured, we need to remap
        
    else:
        return data, x, y


def main(test_files, cntl_files, vname, fig_name, map_file=None, test_name='Model', cntl_name='Obs', **kwargs):

    # Load datasets
    datasets = [open_dataset(f) for f in (test_files, cntl_files)]

    # Subset for overlapping time periods
    t1 = max([ds.time[0].values for ds in datasets])
    t2 = min([ds.time[-1].values for ds in datasets])
    datasets = [ds.sel(time=slice(str(t1), str(t2))) for ds in datasets]
    #t1, t2 = datasets[0].time[0].values, datasets[0].time[-1].values
    #datasets = [ds.sel(time=slice(str(t1), str(t2))) for ds in datasets]
    #t1, t2 = datasets[1].time[0].values, datasets[1].time[-1].values
    #datasets = [ds.sel(time=slice(str(t1), str(t2))) for ds in datasets]

    # Select data
    data_arrays = [get_data(ds, vname) for ds in datasets]
    coords = [get_coords(ds) for ds in datasets]
    weights = [get_area_weights(ds) for ds in datasets]

    # Compute time averages
    means = [da.mean(dim='time', keep_attrs=True) for da in data_arrays]
    weights = [wgt.mean(dim='time', keep_attrs=True) if 'time' in wgt.dims else wgt for wgt in weights]
    lats = [c[1].mean(dim='time', keep_attrs=True) if 'time' in c[1].dims else c[1] for c in coords]

    # Remap data to a lat/lon grid to compute zonal means. 
    #means, lats = [to_latlon(m, coords) for (m, coords) in zip(means, lats)] 
    #weights, lats = [to_latlon(m, coords) for (m, coords) in zip(weights, lats)] 

    # Compute *zonal* average. Note that this is tricky for unstructured data.
    # Our approach is to compute means over latitude bins, and return a new
    # coordinate for latitude of the zonal mean. The function defined in
    # e3sm_utils takes the data, a set of latitude weights, and the input
    # latitudes, and optionally the number of new latitude bands within which to
    # average the data, and returns the binned/averaged data and the new
    # latitudes.
    if map_file is not None:
        print('Apply map...')
        means[0] = apply_map(means[0], map_file, template=means[1])
        weights[0] = apply_map(weights[0], map_file, template=means[1])
        lats[0] = lats[1]
        weights, *__ = zip(*[xarray.broadcast(w, d) for (w, d) in zip(weights, means)])
        means = [zonal_mean(d, weights=w) for (d, w) in zip(means, weights)]
    else:
        print('Try using our slow zonal mean routine...')
        means, lats = zip(*[calculate_zonal_mean(d, w, y) for (d, w, y) in zip(means, weights, lats)])

    # Make line plots of area averages
    figure, ax = pyplot.subplots(1, 1)
    plots = [plot_zonal_mean(y, m, label=l, **kwargs) for (y, m, l) in zip(lats, means, (test_name, cntl_name))]
    pyplot.legend()
    figure.savefig(fig_name, bbox_inches='tight')


if __name__ == '__main__':
    import plac; plac.call(main)
