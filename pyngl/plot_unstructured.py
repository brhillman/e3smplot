#!/usr/bin/env python3

import plac
import ngl
import numpy
import xarray

def plot_unstructured(
        x, y, data, plot_format='x11',
        plot_name='dummy', **kwargs):

    # Setup the canvas
    wks = ngl.open_wks(plot_format, plot_name)

    # Set up annoying plot resources
    res = ngl.Resources()
    res.cnFillOn = True
    res.cnLinesOn = False
    res.cnFillPalette = 'MPL_viridis'

    # If passed 2d coordinate arrays assume they represent cell vertices, 
    # otherwise assume cell centers
    if len(x.shape) == 2:
        res.cnFillMode = 'CellFill'
        res.sfXCellBounds = x
        res.sfYCellBounds = y
    else:
        res.sfXArray = x
        res.sfYArray = y

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


def main(datafile, varname, gridfile=None):

    # Read data
    ds_data = xarray.open_dataset(datafile)
    data = ds_data[varname]
    if gridfile is not None:
        ds_grid = xarray.open_dataset(gridfile).rename({'grid_size': 'ncol'})
        lon = ds_grid['grid_corner_lon']
        lat = ds_grid['grid_corner_lat']
    else:
        lon = ds_data['lon']
        lat = ds_data['lat']

    # Make sure we don't have time or level dimensions
    if 'time' in data.dims:
        data = data.isel(time=0).squeeze()
    if 'lev' in data.dims:
        data = data.isel(lev=-1).squeeze()

    # Make plot
    plot = plot_unstructured(
        lon.values, lat.values, data.values,
        plot_format='png', plot_name=varname, 
        mpGeophysicalLineColor='white',
        lbOrientation='horizontal', 
        lbTitleString='%s (%s)'%(data.long_name, data.units)
    )


if __name__ == '__main__':
    plac.call(main)
