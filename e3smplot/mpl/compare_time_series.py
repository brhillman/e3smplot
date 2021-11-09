#!/usr/bin/env python3
import xarray
import numpy
import sys
from glob import glob
from ..e3sm_utils import get_data, get_coords, get_area_weights, area_average
from ..e3sm_utils import open_dataset, mask_all_zero
from matplotlib import pyplot


def plot_time_series(data, **kwargs):
    ax = pyplot.gca()
    if len(data.shape) == 1:
        pl = ax.plot(data['time'].values, data, **kwargs)
        ax.set_ylabel(f'{data.long_name} ({data.units})')
        ax.set_xlabel('Time')
    else:
        pl = ax.contourf(data['time'].values, data.lev, data.transpose(), **kwargs)
        ax.set_ylabel(f'Vertical level')
        ax.set_xlabel('Time')
        cb = pyplot.colorbar(
            pl, orientation='horizontal',
            label=ax.set_ylabel(f'{data.long_name} ({data.units})')
        )
    pyplot.xticks(rotation=45)
    return pl


def main(files, vname, fig_name, labels, operation='average', time_offsets=None,
        dims=None, t1=None, t2=None, verbose=False, dpi=400, **kwargs):

    # Load datasets
    if verbose: print('Load datasets'); sys.stdout.flush()
    if time_offsets is None: time_offsets = [None for x in range(len(files))]
    datasets = [open_dataset(*f, time_offset=dt, chunks={'time': 1}) for f, dt in zip((files), time_offsets)]

    # Subset for overlapping time periods
    if verbose: print('Subset time'); sys.stdout.flush()
    if t1 is None: t1 = max([ds.time[0].values for ds in datasets])
    if t2 is None: t2 = min([ds.time[-1].values for ds in datasets])
    if verbose: print('Comparing period {} to {}'.format(str(t1), str(t2))); sys.stdout.flush()
    datasets = [ds.sel(time=slice(str(t1), str(t2))) for ds in datasets]

    # Select data
    if verbose: print('Select data'); sys.stdout.flush()
    data_arrays = [mask_all_zero(get_data(ds, vname)) for ds in datasets]
    coords = [get_coords(ds) for ds in datasets]
    weights = [get_area_weights(ds) for ds in datasets]

    # Compute area averages
    if verbose: print('Compute global statistics'); sys.stdout.flush()
    means = [area_average(da, wgt, dims=[d for d in da.dims if d != 'time']) 
             for (da, wgt) in zip(data_arrays, weights)]
    maxes = [da.max(dim=[d for d in da.dims if d != 'time']).values for da in data_arrays]
    mins  = [da.min(dim=[d for d in da.dims if d != 'time']).values for da in data_arrays]

    # Make line plots of area averages
    if verbose: print('Make plots'); sys.stdout.flush()
    figure, ax = pyplot.subplots(1, 1)
    if 'ymin' in kwargs.keys(): 
        ymin = kwargs.pop('ymin')
    else:
        ymin = None
    plots_avg = [plot_time_series(m, label=l, **kwargs) for (m, l) in zip(means, labels)]
    pyplot.legend()
    colors = [p[0].get_color() for p in plots_avg]
    #plots_max = [plot_time_series(m, color=c, linestyle='dashed', **kwargs) for (m, c) in zip(maxes, colors)]
    #plots_min = [plot_time_series(m, color=c, linestyle='dashed', **kwargs) for (m, c) in zip(mins, colors)]
    if ymin is not None: ax.set_ylim([ymin, ax.get_ylim()[1]])
    figure.savefig(fig_name, bbox_inches='tight', dpi=dpi)

    # Clean up
    if verbose: print('Clean up'); sys.stdout.flush()
    pyplot.close()
    for ds in datasets: ds.close()


if __name__ == '__main__':
    import plac; plac.call(main)
