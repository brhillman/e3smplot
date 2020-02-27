#!/usr/bin/env python3

import plac
import ngl
import numpy
import xarray

def plot_unstructured(wks, xv, yv, data, **kwargs):

    # Set up annoying plot resources
    res = ngl.Resources()
    res.cnFillOn = True
    res.cnLinesOn = False
    res.cnFillPalette = 'MPL_viridis'
    res.cnFillMode = 'RasterFill'
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

    return plot


def main(datafile, gridfile, varname):

    # Read data
    ds_data = xarray.open_dataset(datafile)
    ds_grid = xarray.open_dataset(gridfile)
    data = ds_data[varname]

    if 'time' in data.dims:
        data = data.isel(time=0).squeeze()
    if 'lev' in data.dims:
        data = data.isel(lev=-1).squeeze()

    if 'lon' in ds_grid and 'lat' in ds_grid:
        x = ds_grid['lon']
        y = ds_grid['lat']
    elif 'grid_corner_lon' in ds_grid and 'grid_corner_lat' in ds_grid:
        x = ds_grid['grid_corner_lon'].rename({'grid_size': 'ncol'})
        y = ds_grid['grid_corner_lat'].rename({'grid_size': 'ncol'})

    # Setup the canvas
    plot_format='png'
    plot_name='./' + varname
    wks = ngl.open_wks(plot_format, plot_name)

    # Make plot
    plot = plot_unstructured(
        wks, x.values, y.values, data.values,
        mpGeophysicalLineColor='white',
        lbOrientation='horizontal', 
        lbTitleString='%s (%s)'%(data.long_name, data.units),
        cnFillMode='RasterFill',
        cnLineLabelsOn=False, cnLinesOn=False,
    )

    # Close things
    ngl.end()


if __name__ == '__main__':
    plac.call(main)
