#!/usr/bin/env python3

import plac
import ngl
import numpy
import xarray
import os
from e3smplot.utils import nice_cntr_levels
from e3smplot.e3sm_utils import get_data, get_area_weights, area_average
import dask

def plot_map(wks, x, y, data, **kwargs):

    # Set up annoying plot resources
    res = ngl.Resources()
    res.cnFillOn = True
    res.cnLinesOn = False
    res.cnFillPalette = 'MPL_viridis'

    # Add cyclic point if we need to
    if len(data.shape) == 2 and x.max() < 360:
        data, x = ngl.add_cyclic(data, x)

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


def plot_map_fig(x, y, data, plotname, **kwargs):
    # Setup the canvas
    wks = ngl.open_wks(
        os.path.splitext(plotname)[1][1:],
        os.path.splitext(plotname)[0]
    )
    # Make the plot
    return plot_map(wks, x, y, data, **kwargs)


def main(varname, plotname, *datafiles, gridfile=None,
         time_index=None, lev_index=0, functions=None,
         vmin=None, vmax=None, **kwargs):

    # Read data
    ds_data = xarray.open_mfdataset(
        sorted(datafiles), chunks={'time': 1},
        drop_variables=('P3_input_dim', 'P3_output_dim'),
    )
    data = get_data(ds_data, varname)
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
        if time_index is None:
            print("Doing time average...")
            data = data.mean(dim='time', keep_attrs=True).squeeze()
        else:
            print(f"Selecting time index {time_index}...")
            data = data.isel(time=int(time_index))
    if 'lev' in data.dims:
        print(f"Selecting level index {lev_index}...")
        data = data.isel(lev=lev_index).squeeze()
    if 'time' in x.dims: x = x.isel(time=0)
    if 'time' in y.dims: y = y.isel(time=0)

    # Setup the canvas
    wks = ngl.open_wks(
        os.path.splitext(plotname)[1][1:],
        os.path.splitext(plotname)[0]
    )

    if functions is not None: data = eval(functions)

    # Get contour levels; the explicit type casting deals with problems calling
    # this standalone code using subprocess.run() with string arguments, where
    # all kwargs are going to be interpreted as strings
    if vmin is not None and vmax is not None:
        if float(vmin) < 0 and float(vmax) > 0:
            *__, clevels = nice_cntr_levels(float(vmin), float(vmax), returnLevels=True, max_steps=13, aboutZero=True)
        else:
            *__, clevels = nice_cntr_levels(float(vmin), float(vmax), returnLevels=True, max_steps=13)
        kwargs['cnLevels'] = clevels #get_contour_levels(data)
        kwargs['cnLevelSelectionMode'] = 'ExplicitLevels'

    if 'cnLevels' in kwargs.keys():
        kwargs['cnLevels'] = kwargs['cnLevels'].split(',')

    if 'mpMinLonF' in kwargs.keys():
        # subset data for min/max
        mx = (x > float(kwargs['mpMinLonF'])) & (x < float(kwargs['mpMaxLonF']))
        my = (y > float(kwargs['mpMinLatF'])) & (y < float(kwargs['mpMaxLatF']))
        dsub = data.where(mx).where(my)
        print(f'data.min = {dsub.min().values}; data.max = {dsub.max().values}')
        kwargs['tiMainString'] = f'min = {dsub.min().values:.02f}; max = {dsub.max().values:.02f}'

    # Get title string
    if 'tiMainString' not in kwargs.keys():
        dmin = data.min().values
        dmax = data.max().values
        wgts = get_area_weights(ds_data)
        davg = area_average(data, wgts).values
        kwargs['tiMainString'] = f'min = {dmin:.02f}; max = {dmax:.02f}; mean = {davg:.02f}'

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
    ngl.destroy(wks)

    # Trim whitespace
    os.system(f"convert -trim {plotname} {plotname}")

    # Fix permissions
    os.system(f'chmod a+r {plotname}')


if __name__ == '__main__':
    plac.call(main)
