#!/usr/bin/env python3

from plot_utils import open_files, get_data, area_average
from matplotlib import pyplot

#
# Dictionary of x and y axis names
#
coord_names = {
    'CFAD_DBZE94_CS': ('cosp_dbze', 'cosp_ht'),
    'CFAD_SR532_CAL': ('cosp_sr', 'cosp_ht'),
    'CLD_MISR': ('cosp_tau', 'cosp_htmisr'),
    'FISCCP1': ('cosp_tau', 'cosp_prs'),
    'CLMODIS': ('cosp_tau_modis', 'cosp_prs'),
}

def plot_cfad(xbnds, ybnds, data, ax=None, add_colorbar=True, **kwargs):
    #
    # Select axes
    #
    if ax is None: ax = pyplot.gca()
    #
    # Make pcolor plot
    #
    pl = ax.pcolormesh(xbnds, ybnds, data, **kwargs)
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
    ds = open_files(*inputfiles)
    data = get_data(ds, varname)
    area = get_data(ds, 'area')
    x = get_data(ds, coord_names[varname][0])
    y = get_data(ds, coord_names[varname][1])
    #
    # Reduce data
    #
    data = data.mean(dim='time', keep_attrs=True)
    data_mean = area_average(data, area, dims='ncol')
    #
    # Plot
    #
    figure = pyplot.figure()
    ax = figure.add_subplot(111)
    pl, cb = plot_cfad(x, y, data_mean.transpose(*coord_names[varname][::-1]), **kwargs)
    figure.savefig(outputfile, bbox_inches='tight')


if __name__ == '__main__':
    import plac; plac.call(main)
