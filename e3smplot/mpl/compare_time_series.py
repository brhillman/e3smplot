#!/usr/bin/env python3
import xarray
import numpy
from glob import glob
from ..e3sm_utils import get_data, get_coords, get_area_weights, area_average
from matplotlib import pyplot


def open_dataset(files, **kwargs):
    # Open dataset as a dask array
    ds = xarray.open_mfdataset(sorted(files), combine='by_coords', drop_variables='P3_output_dim', **kwargs) 

    # Rename coordinate variables
    if 'latitude' in ds.dims: ds = ds.rename({'latitude': 'lat'})
    if 'longitude' in ds.dims: ds = ds.rename({'longitude': 'lon'})

    if 'msshf' in ds.variables.keys():
        ds['msshf'].values = -ds['msshf'].values

    # Return dataset
    return ds


def plot_time_series(data, **kwargs):
    ax = pyplot.gca()
    pl = ax.plot(data['time'].values, data, **kwargs)
    return pl


def subset_time(datasets):

    # Subset based on date range of first dataset
    t1, t2 = datasets[0].time[0].values, datasets[0].time[-1].values
    for i in range(1, len(datasets)):
        datasets[i] = datasets[i].sel(time=slice(str(t1), str(t2)))
    return datasets


def main(test_files, cntl_files, vname, fig_name, test_name='Model', cntl_name='Obs', **kwargs):

    # Load datasets
    datasets = [open_dataset(f) for f in (test_files, cntl_files)]


    # Subset for overlapping time periods
    t1, t2 = datasets[0].time[0].values, datasets[0].time[-1].values
    datasets = [ds.sel(time=slice(str(t1), str(t2))) for ds in datasets]

    # Select data
    data_arrays = [get_data(ds, vname) for ds in datasets]
    coords = [get_coords(ds) for ds in datasets]
    weights = [get_area_weights(ds) for ds in datasets]

    # Compute area averages
    means = [area_average(da, wgt, dims=[d for d in da.dims if d != 'time']) for (da, wgt) in zip(data_arrays, weights)]

    # Try to convert time to datetime64
    means[0]['time'] = means[0].indexes['time'].to_datetimeindex()

    for m in means:
        print(m.time[0], m.time[-1])

    # Make line plots of area averages
    figure, ax = pyplot.subplots(1, 1)
    plots = [plot_time_series(m, label=l, **kwargs) for (m, l) in zip(means, (test_name, cntl_name))]
    pyplot.legend()
    figure.savefig(fig_name, bbox_inches='tight')


if __name__ == '__main__':
    import plac; plac.call(main)
