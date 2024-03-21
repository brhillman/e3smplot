#!/usr/bin/env python3

import ngl
import xarray
import os
from e3smplot.utils import nice_cntr_levels
from e3smplot.e3sm_utils import get_data

def get_cloud_cmap():

    # Start with MPL_Gray colormap

def plot_map_overlay(ds, image_file, **kwargs):

    # Open workstation

    # Plot image
    pl_img = plot_map_image(wks, imgfile, nglDraw=False, nglFrame=False, **kwargs)

    # Plot data; need to plot as a regular contour plot so we can overlay with
    # background image above
    cmap = get_cloud_cmap()
    pl_data = plot_map(
        wks, x.values, y.values, data.values, ngl_function=ngl.contour,
        nglDraw=False, nglFrame=False,
        **kwargs
    )

    ngl.overlay(pl_img, pl_data)
    ngl.draw(wks)
    ngl.frame(wks)




def main(varname, plotname, imgfile, *datafiles, gridfile=None,
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

    #pl = plot_over_image(wks, x, y, data, image_file, **kwargs)

    # Setup the canvas
    wks = ngl.open_wks(
        os.path.splitext(plotname)[1][1:],
        os.path.splitext(plotname)[0]
    )

    # Plot background image, but do not advance frame or draw because we will
    # need to overlay the data on top
    pl_img = plot_map_image(wks, imgfile, nglDraw=False, nglFrame=False, **kwargs)

    # Get contour levels; the explicit type casting deals with problems calling
    # this standalone code using subprocess.run() with string arguments, where
    # all kwargs are going to be interpreted as strings
    if vmin is None: vmin = data.min().values
    if vmax is None: vmax = data.max().values
    aboutZero = float(vmin) < 0 and float(vmax) > 0
    *__, clevels = nice_cntr_levels(float(vmin), float(vmax), returnLevels=True,
                                    max_steps=ncontours, aboutZero=aboutZero)
    kwargs['cnLevels'] = clevels
    kwargs['cnLevelSelectionMode'] = 'ExplicitLevels'

    # Make plot; need to plot as a regular contour plot so we can overlay with
    # background image above
    pl_data = plot_map(
        wks, x.values, y.values, data.values, ngl_function=ngl.contour,
        cnLineLabelsOn=False, cnLinesOn=False,
        nglDraw=False, nglFrame=False,
        **kwargs
    )
    ngl.overlay(pl_img, pl_data)
    ngl.draw(wks)
    ngl.frame(wks)

    # Close things
    ngl.end()


if __name__ == '__main__':
    import plac; plac.call(main)
