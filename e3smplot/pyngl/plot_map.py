#!/usr/bin/env python3

import plac
import ngl
import numpy
import xarray

def plot_map(wks, x, y, data, **kwargs):

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
        res.cnFillMode = 'RasterFill'
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

    return plot


def main(datafile, varname, plotname, gridfile=None, **kwargs):

    # Read data
    ds_data = xarray.open_dataset(datafile)
    data = ds_data[varname]
    if gridfile is not None:
        ds_grid = xarray.open_dataset(gridfile).rename({'grid_size': 'ncol'})
        if 'lon' in ds_grid and 'lat' in ds_grid:
            x = ds_grid['lon']
            y = ds_grid['lat']
        elif 'grid_corner_lon' in ds_grid and 'grid_corner_lat' in ds_grid:
            x = ds_grid['grid_corner_lon']
            y = ds_grid['grid_corner_lat']
        else:
            raise RuntimeError('No valid coordinates in grid file.')
    else:
        x = ds_data['lon']
        y = ds_data['lat']

    # Make sure we don't have time or level dimensions
    if 'time' in data.dims:
        data = data.isel(time=0).squeeze()
    if 'lev' in data.dims:
        data = data.isel(lev=-1).squeeze()

    # Setup the canvas
    plot_format='png'
    wks = ngl.open_wks(plot_format, plotname)

    # Make plot
    plot = plot_map(
        wks, x.values, y.values, data.values,
        mpGeophysicalLineColor='white',
        lbOrientation='horizontal', 
        lbTitleString='%s (%s)'%(data.long_name, data.units),
        cnFillMode='RasterFill',
        cnLineLabelsOn=False, cnLinesOn=False, **kwargs
    )

    # Close things
    ngl.end()


if __name__ == '__main__':
    plac.call(main)
