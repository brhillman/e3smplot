#!/usr/bin/env python3

import xarray, ngl, os
from e3smplot.e3sm_utils import get_data, get_area_weights, area_average
from e3smplot.pyngl.plot_map import plot_map
from e3smplot.utils import nice_cntr_levels

def compare_maps(wks, xcoords, ycoords, data_arrays, mapfile=None, labels=None, **kwargs):

    # remap or whatever

    # Append difference matrix
    data_arrays.append(data_arrays[0] - data_arrays[1])
    data_arrays[-1].attrs = data_arrays[0].attrs
    xcoords.append(xcoords[0])
    ycoords.append(ycoords[0])

    if 'clevels' not in kwargs.keys():
        vmin = min(da.min().values for da in data_arrays[:-1])
        vmax = max(da.max().values for da in data_arrays[:-1])
        *_, clevels = nice_cntr_levels(vmin, vmax, outside=True, max_steps=15, cint=None,
                         returnLevels=True, aboutZero=False)
        kwargs['clevels'] = clevels
    else:
        kwargs['clevels'] = kwargs['clevels'].split(',')
    if 'dlevels' not in kwargs.keys():
        vmax = abs(data_arrays[-1]).max().values
        vmin = -vmax
        *_, dlevels = nice_cntr_levels(vmin, vmax, outside=True, max_steps=15, cint=None,
                         returnLevels=True, aboutZero=True)
        kwargs['dlevels'] = dlevels
    else:
        kwargs['dlevels'] = kwargs['dlevels'].split(',')

    # Make the plot
    kwargs.update(dict(nglDraw=False, nglFrame=False))
    kws = [{k: v for k, v in kwargs.items()} for dummy in data_arrays]
    kws[-1]['cnFillPalette'] = 'MPL_coolwarm'
    if 'clevels' in kwargs.keys():
        del kws[0]['clevels']
        del kws[1]['clevels']
        del kws[2]['clevels']
        kws[0]['cnLevels'] = clevels
        kws[1]['cnLevels'] = clevels
        kws[0]['cnLevelSelectionMode'] = 'ExplicitLevels'
        kws[1]['cnLevelSelectionMode'] = 'ExplicitLevels'
    if 'dlevels' in kwargs.keys():
        del kws[0]['dlevels']
        del kws[1]['dlevels']
        del kws[2]['dlevels']
        kws[2]['cnLevels'] = dlevels
        kws[2]['cnLevelSelectionMode'] = 'ExplicitLevels'

    if labels is not None:
        [kw.update(dict(tiMainString=l)) for kw, l in zip(kws, labels)]
    print(kws)
    plots = [plot_map(wks, x.values, y.values, d.values, **kw) for (x, y, d, kw)
             in zip(xcoords, ycoords, data_arrays, kws)]
    # Panel the plots together
    ngl.panel(wks, plots, [3, 1])

    return plots


def get_coords(ds):
    if 'lon' in ds and 'lat' in ds:
        x = ds['lon']
        y = ds['lat']
    elif 'grid_corner_lon' in ds and 'grid_corner_lat' in ds:
        x = ds['grid_corner_lon']
        y = ds['grid_corner_lat']
    elif 'grid_center_lon' in ds and 'grid_center_lat' in ds:
        x = ds['grid_center_lon']
        y = ds['grid_center_lat']
    else:
        raise RuntimeError('No valid coordinates in grid file.')
    return x, y


def main(vname, plotname, testfiles, cntlfiles, gridfiles=(None, None),
         mapfile=None, time_indices=(None, None), **kwargs):

    # Open datasets
    datasets = [xarray.open_mfdataset(f) for f in (testfiles, cntlfiles)]

    # Select time index or average
    datasets = [ds.isel(time=i) if i is not None else
                ds.mean(dim='time', keep_attrs=True) for (ds, i) in
                zip(datasets, time_indices)]

    # Open grid datasets
    grid_datasets = [ds if grid_file is None else
                     xarray.open_dataset(grid_file).rename({'grid_size': 'ncol'})
                     for ds, grid_file in zip(datasets, gridfiles)]

    # Select data
    data_arrays = [get_data(ds, vname) for ds in datasets]
    area_weights = [get_area_weights(ds) for ds in datasets]
    xcoords, ycoords = map(list, zip(*[get_coords(ds) for ds in grid_datasets]))

    # Remap, if needed
    #if mapfile is not None:

    # Setup the canvas
    wks = ngl.open_wks(
        os.path.splitext(plotname)[1][1:],
        os.path.splitext(plotname)[0]
    )

    # Append title strings
    labels = [f'min = {da.min().values:.02f}; max = {da.max().values:.02f}; mean = {area_average(da, wgt).values:.02f}'
              for da, wgt in zip(data_arrays, area_weights)]

    # Create figure
    fig = compare_maps(
        wks, xcoords, ycoords, data_arrays,
        mapfile=mapfile, labels=labels,
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
    import plac; plac.call(main)
