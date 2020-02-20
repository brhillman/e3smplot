#!/usr/bin/env python3

import plac
import ngl
import numpy
import xarray

def plot_unstructured(
        xv, yv, data, plot_format='x11',
        plot_name='dummy', **kwargs):

    # Setup the canvas
    wks = ngl.open_wks(plot_format, plot_name)

    # Set up annoying plot resources
    res = ngl.Resources()
    res.cnFillOn = True
    res.cnLinesOn = False
    res.cnFillPalette = 'MPL_viridis'
    res.cnFillMode = 'CellFill'
    res.sfXCellBounds = xv
    res.sfYCellBounds = yv

    # Tweak plot appearance
    res.mpGridAndLimbOn = False
    res.mpPerimOn = False

    # Additional options passed via kwargs
    for key, val in kwargs.items():
        setattr(res, key, val)

    # Make the plot
    plot = ngl.contour_map(wks, data, res)
    ngl.end()

    return plot


def main(datafile, gridfile, varname):
    # Read data
    ds_data = xarray.open_dataset(datafile)
    ds_grid = xarray.open_dataset(gridfile).rename({'grid_size': 'ncol'})
    data = ds_data[varname]
    lon_edges = ds_grid['grid_corner_lon']
    lat_edges = ds_grid['grid_corner_lat']

    # Make plot
    plot = plot_unstructured(
        lon_edges.values, lat_edges.values, data.values,
        plot_format='png', plot_name=varname, 
        mpGeophysicalLineColor='white',
        lbOrientation='horizontal', 
        lbTitleString='%s (%s)'%(data.long_name, data.units)
    )


if __name__ == '__main__':
    plac.call(main)
