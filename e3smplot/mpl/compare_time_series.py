#!/usr/bin/env python3
import xarray
import numpy
from glob import glob
from ..e3sm_utils import get_data, get_coords, get_area_weights, area_average
from ..e3sm_utils import open_dataset, mask_all_zero
from matplotlib import pyplot


def plot_time_series(data, **kwargs):
    ax = pyplot.gca()
    pl = ax.plot(data['time'].values, data, **kwargs)
    ax.set_ylabel(f'{data.long_name} ({data.units})')
    ax.set_xlabel('Time')
    pyplot.xticks(rotation=45)
    return pl


def main(files, vname, fig_name, labels, time_offsets=None, **kwargs):

    # Load datasets
    if time_offsets is None: time_offsets = [None for x in range(len(files))]
    datasets = [open_dataset(f, time_offset=dt) for f, dt in zip((files), time_offsets)]

    # Subset for overlapping time periods
    t1 = max(ds.time[0].values for ds in datasets)
    t2 = min(ds.time[-1].values for ds in datasets)
    datasets = [ds.sel(time=slice(str(t1), str(t2))) for ds in datasets]

    # Select data
    data_arrays = [mask_all_zero(get_data(ds, vname)) for ds in datasets]
    coords = [get_coords(ds) for ds in datasets]
    weights = [get_area_weights(ds) for ds in datasets]

    # Compute area averages
    means = [area_average(da, wgt, dims=[d for d in da.dims if d != 'time']) for (da, wgt) in zip(data_arrays, weights)]

    # Make line plots of area averages
    figure, ax = pyplot.subplots(1, 1)
    plots = [plot_time_series(m, label=l, **kwargs) for (m, l) in zip(means, labels)]
    pyplot.legend()
    figure.savefig(fig_name, bbox_inches='tight')

    # Clean up
    pyplot.close()
    for ds in datasets: ds.close()


if __name__ == '__main__':
    import plac; plac.call(main)
