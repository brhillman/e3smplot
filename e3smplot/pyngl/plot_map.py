#!/usr/bin/env python3

import ngl
import xarray
import os
from e3smplot.utils import nice_cntr_levels
from e3smplot.e3sm_utils import get_data


def plot_map(wks, x, y, data, ngl_function=ngl.contour_map, **kwargs):

    # Set up default plot resources as dict to allow overriding via kwargs
    resDict = dict(
        cnFillOn=True, cnLinesOn=False, cnFillPalette='MPL_viridis',
        mpGridAndLimbOn=False, mpPerimOn=False,
    )

    # Add cyclic point if we need to
    if len(data.shape) == 2 and x.max() < 360:
        data, x = ngl.add_cyclic(data, x)

    # If passed 2d coordinate arrays assume they represent cell vertices, 
    # otherwise assume cell centers
    if len(x.shape) == 2:
        resDict['cnFillMode'] = 'CellFill'
        resDict['sfXCellBounds'] = x
        resDict['sfYCellBounds'] = y
    else:
        resDict['cnFillMode'] = 'RasterFill'
        resDict['sfXArray'] = x
        resDict['sfYArray'] = y

    # Additional options passed via kwargs
    for key, val in kwargs.items(): resDict[key] = val

    # Set up plot resources for call to ngl functions
    res = ngl.Resources()
    for key, val in resDict.items(): setattr(res, key, val)

    # Make the plot
    plot = ngl_function(wks, data, res)

    return plot


def main(varname, plotname, *datafiles, gridfile=None,
         time_index=None, time_value=None,
         vmin=None, vmax=None, ncontours=13, **kwargs):

    # Read data
    ds_data = xarray.open_mfdataset(datafiles)
    data = get_data(ds_data, varname)
    if gridfile is not None:
        ds_grid = xarray.open_dataset(gridfile).rename({'grid_size': 'ncol'})
        if 'grid_corner_lon' in ds_grid and 'grid_corner_lat' in ds_grid:
            x = ds_grid['grid_corner_lon']
            y = ds_grid['grid_corner_lat']
        elif 'lon' in ds_grid and 'lat' in ds_grid:
            x = ds_grid['lon']
            y = ds_grid['lat']
        else:
            raise RuntimeError('No valid coordinates in grid file.')
    else:
        x = ds_data['lon']
        y = ds_data['lat']

    # Make sure we don't have time or level dimensions
    if 'time' in data.dims:
        if time_value is not None:
            data = data.sel(time=time_value)
        elif time_index is not None:
            data = data.isel(time=int(time_index))
        else:
            data = data.mean(dim='time', keep_attrs=True).squeeze()
    if 'lev' in data.dims:
        data = data.isel(lev=-1).squeeze()
    if 'time' in x.dims: x = x.isel(time=0)
    if 'time' in y.dims: y = y.isel(time=0)

    # Setup the canvas
    wks = ngl.open_wks(
        os.path.splitext(plotname)[1][1:],
        os.path.splitext(plotname)[0]
    )

    # Get contour levels; the explicit type casting deals with problems calling
    # this standalone code using subprocess.run() with string arguments, where
    # all kwargs are going to be interpreted as strings
    if vmin is None: vmin = data.min().values
    if vmax is None: vmax = data.max().values
    if float(vmin) < 0 and float(vmax) > 0:
        *__, clevels = nice_cntr_levels(float(vmin), float(vmax), returnLevels=True, max_steps=ncontours, aboutZero=True)
    else:
        *__, clevels = nice_cntr_levels(float(vmin), float(vmax), returnLevels=True, max_steps=ncontours)
    kwargs['cnLevels'] = clevels #get_contour_levels(data)
    kwargs['cnLevelSelectionMode'] = 'ExplicitLevels'

    # Make plot
    if 'lbTitleString' not in kwargs.keys():
        kwargs['lbTitleString'] = f'{data.long_name} ({data.units})'
    plot = plot_map(
        wks, x.values, y.values, data.values,
        mpGeophysicalLineColor='white',
        lbOrientation='horizontal', 
        cnLineLabelsOn=False, cnLinesOn=False, **kwargs
    )

    # Close things
    ngl.end()


if __name__ == '__main__':
    import plac; plac.call(main)
