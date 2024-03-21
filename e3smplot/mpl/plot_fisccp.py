#!/usr/bin/env python3

from plot_utils import open_files, get_data, area_average
from matplotlib import pyplot

xbnds = (0, 0.3, 1.3, 3.6, 9.4, 23, 60, 100000)
ybnds = (100000, 80000, 68000, 56000, 44000, 31000, 18000, 0)

def plot_fisccp(data, ax=None, add_colorbar=True, **kwargs):
    #
    # Select axes
    #
    if ax is None: ax = pyplot.gca()
    #
    # Make pcolor plot
    #
    x = range(len(xbnds)-1)
    y = range(len(ybnds)-1)
    pl = ax.pcolormesh(x, y, data, **kwargs)
    ax.set_xticklabels([f'{d}' for d in xbnds])
    ax.set_yticklabels([f'{d}' for d in ybnds])
    #
    # Add colorbar to plot
    #
    if add_colorbar == True:
        cb = pyplot.colorbar(
            pl, ax=ax, orientation='horizontal',
            label=f'{data.long_name} ({data.units})',
            shrink=0.8
        )
    #
    # Return
    #
    return pl, cb

    
def main(varname, outputfile, *inputfiles, **kwargs):

    #
    # Read data
    #
    ds = open_files(*inputfiles).mean(dim='time', keep_attrs=True)
    data = get_data(ds, varname)
    area = get_data(ds, 'area')
    #x = get_data(ds, coord_names[varname][0])
    #y = get_data(ds, coord_names[varname][1])
    #
    # Reduce data
    #
    data_mean = area_average(data, area, dims='ncol')
    #
    # Plot
    #
    figure = pyplot.figure()
    ax = figure.add_subplot(111)
    pl, cb = plot_fisccp(data_mean.transpose('cosp_prs', 'cosp_tau'), **kwargs)
    figure.savefig(outputfile, bbox_inches='tight')


if __name__ == '__main__':
    import plac; plac.call(main)
