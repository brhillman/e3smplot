#!/usr/bin/env python3
import xarray
import numpy
import sys
from glob import glob
from e3smplot.e3sm_utils import get_data, get_coords, get_area_weights, area_average
from e3smplot.e3sm_utils import open_dataset, mask_all_zero
from matplotlib import pyplot

def plot_profile(d, *args, **kwargs):
    if 'ax' in kwargs:
        ax = kwargs.pop('ax')
    else:
        ax = pyplot.gca()
    if 'lev' in d.dims:
        z = d['lev']
    elif 'ilev' in d.dims:
        z = d['ilev']
    elif 'level' in d.dims:
        z = d['level']
    elif 'plev' in d.dims:
        z = d['plev']
    else:
        raise RuntimeError(f'We cannot find a valid level coordinate in {d.dims}')
    pl = ax.plot(d, z, *args, **kwargs)
    if 'hybrid level at' in z.long_name:
        ax.set_ylabel(f'Reference level ({z.units})')
    else:
        ax.set_ylabel(f'{z.long_name} ({z.units})')
    # Reverse direction of axis for pressure coordinates
    ax.set_ylim(max(ax.get_ylim()), min(ax.get_ylim()))
    return pl


def compare_profiles(data_arrays, labels, **kwargs):
    figure = pyplot.figure()
    ax = figure.add_subplot(111)
    ax_diff = ax.twiny()
    linestyles = ('solid', 'dashed', 'dotted')
    for ic, (d, l) in enumerate(zip(data_arrays, labels)):
        pl = plot_profile(d, label=l, ax=ax, **kwargs) #ax.plot(d, range(len(d)), label=l, **kwargs)
        if ic == 0:
            d_cntl = d.copy(deep=True)
        else:
            d_diff = d - d_cntl
            pl_diff = plot_profile(d_diff, label=f'{l} - {labels[0]}', ax=ax_diff, color='0.5', linestyle=linestyles[ic], **kwargs)
    ax.set_xlabel(f'{d.long_name} ({d.units})')
    ax.legend(loc='upper right')
    ax_diff.set_xlabel('Difference')
    ax_diff.legend(loc='lower left')

    # Set axis colors for differences to match line color used
    ax_diff.spines['top'].set_color('0.5')
    ax_diff.xaxis.label.set_color('0.5')
    ax_diff.tick_params(axis='x', colors='0.5')

    return figure 
       

def main(test_files, cntl_files, vname, fig_name, labels, operation='average', time_offsets=None,
        dims=None, t1=None, t2=None, verbose=False, dpi=400, **kwargs):

    # Load datasets
    if verbose: print('Load datasets'); sys.stdout.flush()
    if time_offsets is None: time_offsets = [None for x in range(len(files))]
    datasets = [open_dataset(f, time_offset=dt, chunks={'time': 1}) for f, dt in zip((test_files, cntl_files), time_offsets)]

    # Subset for overlapping time periods
    if verbose: print('Subset time'); sys.stdout.flush()
    if t1 is None: t1 = max([ds.time[0].values for ds in datasets])
    if t2 is None: t2 = min([ds.time[-1].values for ds in datasets])
    if verbose: print('Comparing period {} to {}'.format(str(t1), str(t2))); sys.stdout.flush()
    datasets = [ds.sel(time=slice(str(t1), str(t2))) for ds in datasets]

    # Compute time average
    if verbose: print('Compute time average'); sys.stdout.flush()
    datasets = [ds.mean(dim='time', keep_attrs=True) for ds in datasets]

    # Select data
    if verbose: print('Select data'); sys.stdout.flush()
    data_arrays = [mask_all_zero(get_data(ds, vname)) for ds in datasets]
    coords = [get_coords(ds) for ds in datasets]
    weights = [get_area_weights(ds) for ds in datasets]

    # Compute area averages
    if verbose: print('Compute global statistics'); sys.stdout.flush()
    exclude_dims = ('time', 'lev', 'level', 'plev')
    means = [area_average(da, wgt, dims=[d for d in da.dims if d not in exclude_dims]) 
             for (da, wgt) in zip(data_arrays, weights)]
    maxes = [da.max(dim=[d for d in da.dims if d not in exclude_dims]).values for da in data_arrays]
    mins  = [da.min(dim=[d for d in da.dims if d not in exclude_dims]).values for da in data_arrays]

    # Make plots of area averages
    if verbose: print('Make plots'); sys.stdout.flush()
    figure, ax = pyplot.subplots(1, 1)
    plots_avg = [plot_profile(m, label=l, **kwargs) for (m, l) in zip(means, labels)]
    pyplot.legend()
    colors = [p[0].get_color() for p in plots_avg]
    #plots_max = [plot_time_series(m, color=c, linestyle='dashed', **kwargs) for (m, c) in zip(maxes, colors)]
    #plots_min = [plot_time_series(m, color=c, linestyle='dashed', **kwargs) for (m, c) in zip(mins, colors)]
    ax.set_ylim(ax.get_ylim()[::-1])
    figure.savefig(fig_name, bbox_inches='tight', dpi=dpi)

    # Clean up
    if verbose: print('Clean up'); sys.stdout.flush()
    pyplot.close()
    for ds in datasets: ds.close()


if __name__ == '__main__':
    import plac; plac.call(main)
