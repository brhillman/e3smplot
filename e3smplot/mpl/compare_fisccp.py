#!/usr/bin/env python3

from e3smplot.plot_utils import get_data, area_average
from matplotlib import pyplot
from pylab import cm

xbnds = (0, 0.3, 1.3, 3.6, 9.4, 23, 60, 999)
ybnds = (1000, 800, 680, 560, 440, 310, 180, 0)

def plot_fisccp(data, ax=None, add_colorbar=True, **kwargs):
    #
    # Select axes
    #
    if ax is None: ax = pyplot.gca()
    #
    # Make pcolor plot
    #
    x = range(len(xbnds))
    y = range(len(ybnds))
    pl = ax.pcolormesh(x, y, data, **kwargs)
    xlabels = [f'{d}' for d in xbnds]
    ylabels = [f'{d}' for d in ybnds]
    ax.set_xticks(x); ax.set_xticklabels(xlabels)
    ax.set_yticks(y); ax.set_yticklabels(ylabels)
    ax.set_xlabel('Cloud optical depth')
    ax.set_ylabel('Cloud top pressure (hPa)')
    #
    # Add colorbar to plot
    #
    if add_colorbar == True:
        cb = pyplot.colorbar(
            pl, ax=ax, orientation='horizontal',
            label=f'Frequency ({data.units})',
            #shrink=0.8
        )
    #
    # Return
    #
    return pl, cb

import xarray 
def main(varname, outputfile, testfiles, cntlfiles, casenames=('test', 'cntl'), **kwargs):

    # Setup figure
    figure, axes = pyplot.subplots(1, 3, sharex=True, sharey=True)

    # Read data
    datasets = [xarray.open_mfdataset(inputfiles).mean(dim='time', keep_attrs=True) for inputfiles in (testfiles, cntlfiles)]
    dataarrays = [get_data(ds, varname) for ds in datasets]
    areaarrays = [ds['area'] for ds in datasets]

    # Compute spatial means
    datameans = [area_average(d, a, dims='ncol') for d, a
                 in zip(dataarrays, areaarrays)]

    # Get data range
    vmin = min([d.min().values for d in datameans])
    vmax = max([d.max().values for d in datameans])

    # Compute diff
    datameans.append(datameans[0] - datameans[1])
    datameans[2].attrs = datameans[0].attrs

    for icase, (case, data) in enumerate(zip(casenames, datameans)):

        # Plot
        ax = figure.add_axes(axes[icase])
        cmap = cm.get_cmap('viridis', 11)
        pl, cb = plot_fisccp(data.transpose('cosp_prs', 'cosp_tau'), vmin=vmin,
                             vmax=vmax, cmap=cmap, **kwargs)

        ax.set_title(case)

    # Make difference plot
    dmax = abs(datameans[2]).max().values
    dmin = -dmax
    ax = figure.add_axes(axes[-1])
    cmap = cm.get_cmap('bwr', 11)
    pl, cb = plot_fisccp(datameans[2].transpose('cosp_prs', 'cosp_tau'), vmin=dmin,
                         vmax=dmax, cmap=cmap, **kwargs)
    ax.set_title(f'{casenames[0]} - {casenames[1]}') 
    axes[1].set_ylabel('')
    axes[2].set_ylabel('')

    # Save figure
    figure.savefig(outputfile, bbox_inches='tight')


if __name__ == '__main__':
    import plac; plac.call(main)
