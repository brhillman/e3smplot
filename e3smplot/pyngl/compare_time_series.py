#!/usr/bin/env python3
import xarray
import numpy
from glob import glob
from ..e3sm_utils import get_data, get_coords, get_area_weights, area_average
import ngl

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


def plot_time_series(time, data, **kwargs):

    # Setup a plot resources object
    res = ngl.Resources()

    # Set resource values by keyword
    for key, val in kwargs.items():
        setattr(res, key, val)

    # Make plot
    pl = ngl.xy(wks, time, data, res)
    return pl


def main(test_files, cntl_files, vname, fig_name, test_name='Model', cntl_name='Obs', **kwargs):

    # Load datasets
    datasets = [open_dataset(f) for f in (test_files, cntl_files)]
    data_arrays = [get_data(ds, vname) for ds in datasets]
    coords = [get_coords(ds) for ds in datasets]
    weights = [get_area_weights(ds) for ds in datasets]

    # Compute area averages
    means = [area_average(da, wgt, dims=[d for d in da.dims if d != 'time']) for (da, wgt) in zip(data_arrays, weights)]

    # Make line plots of area averages
    times = [da['time'] for da in data_arrays]
    labels = (test_name, cntl_name)
    plot_format = fig_name.split('.')[-1]
    wks = ngl.open_wks(plot_format, fig_name)


if __name__ == '__main__':
    import plac; plac.call(main)
