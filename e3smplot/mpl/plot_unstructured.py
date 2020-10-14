#!/usr/bin/env python3

from matplotlib import pyplot
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection, PolyCollection
from cartopy import crs
import numpy
import xarray


def plot_unstructured(
        xv, yv, data, antialiased=False, vmin=None, vmax=None, **kwargs):
    """
    Plot unstructured data. xv and yv specify the x and y coordinates
    of the vertices of each cell and should have shape [ni,nv] where ni
    is the size of the grid (number of cells) and nv is the number of
    vertices for each cell. Data contains the values of the data you
    want to plot, and should have dimension [ni,]. The intent with this
    function is that you will need to read in an auxillary SCRIP-format
    file that describes the grid (unless this information is present in
    the same file that contains the data; unlikely) as well as the file
    that contains the data, and grab the cell corner information from
    the SCRIP file and the data from the data file. This function will
    then plot the data on the native grid by looping over each cell and
    drawing patches for each. Note that this will probably be really
    slow for big grids! Seems to work alright up to about ne120 or so,
    but really some more clever techniques should probably be used here
    (parallelism?).
    
    NOTE: To avoid artifacts due to antialiasing, you should probably pass
    antialiaseds=False to **kwargs.
    """

    # Create array of cell vertices, indexed [npoints, ncorners, 2]
    corners = numpy.stack([xv, yv], axis=2)
    
    # Go back and fix corners where they wrap; we shouldn't have to do
    # this with cartopy, but it seems we do...
    for i in range(corners.shape[0]):
        if any(corners[i,:,0] < -90) and any(corners[i,:,0] > 90):
            corners[i,:,0] = numpy.where(corners[i,:,0] < -90, corners[i,:,0] + 360, corners[i,:,0])
        if any(corners[i,:,1] < -45) and any(corners[i,:,1] > 45):
            corners[i,:,1] = numpy.where(corners[i,:,1] < -45, corners[i,:,1] + 90, corners[i,:,1])

    # Create a PatchCollection from our aggregated list of PathPatches
    p = PolyCollection(corners, array=data, antialiaseds=antialiased, **kwargs)

    # Set scale, mimicking vmin and vmax plot kwargs
    if vmin is not None and vmax is not None:
        p.set_clim([vmin, vmax])

    # Add the collection to the axes
    ax = pyplot.gca()
    ax.add_collection(p)

    # Set sane axes limits
    ax.set_xlim([xv.min(), xv.max()])
    ax.set_ylim([yv.min(), yv.max()])
    
    # Return collection of patches
    return p


def fix_lon(lon):
    return numpy.where(lon > 180, lon - 360, lon)
 

def main(datafile, gridfile, varname, plotfile=None):
    # Read data
    ds_data = xarray.open_dataset(datafile)
    ds_grid = xarray.open_dataset(gridfile).rename({'grid_size': 'ncol'})
    data = ds_data[varname]
    lon_edges = ds_grid['grid_corner_lon']
    lat_edges = ds_grid['grid_corner_lat']

    # Reduce data if we need to
    if 'time' in data.dims: data = data.isel(time=0).squeeze()
    if 'lev' in data.dims: data = data.isel(lev=-1).squeeze()

    # Make plot
    figure, ax = pyplot.subplots(
        1, 1,
        #subplot_kw=dict(projection=crs.Orthographic(central_longitude=-100))
        subplot_kw=dict(projection=crs.PlateCarree())
    )
    pl = plot_unstructured(
        fix_lon(lon_edges.values), lat_edges.values, data.values,
        transform=crs.PlateCarree()
    )
    ax.set_global()
    ax.coastlines(color='white', linewidth=0.5)

    # Add colorbar to plot
    cb = pyplot.colorbar(
        pl, orientation='horizontal', 
        label='%s (%s)'%(data.long_name, data.units)
    )

    # Save plot
    if plotfile is None:
        plotfile = '%s.png'%varname

    figure.savefig(plotfile, bbox_inches='tight')


if __name__ == '__main__':
    import plac; plac.call(main)
