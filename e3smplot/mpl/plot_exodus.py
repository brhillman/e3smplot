#!/usr/bin/env python3

def plot_exodus(ds, ax=None, **plot_kwargs):
    import numpy
    from matplotlib import pyplot
    from cartopy import crs, feature

    # set up the axes and map
    if ax is None: 
        ax = pyplot.subplot(projection=crs.PlateCarree())
        ax.set_global()

        # add ocean and land fills
        ax.add_feature(feature.OCEAN, zorder=0)
        ax.add_feature(feature.LAND, zorder=0)

    # cartesian coordinates (radians?)
    x = ds['coord'][0,:].squeeze()
    y = ds['coord'][1,:].squeeze()
    z = ds['coord'][2,:].squeeze()

    # convert to longitude and latitude
    lon = numpy.arctan2(y, x) * 180.0 / numpy.pi
    lat = numpy.arcsin(z) * 180.0 / numpy.pi

    print('Number of vertices: %i'%lon.size)
    print('Longitude range: %i to %i'%(lon.min().values, lon.max().values))
    print('Latitude range: %i to %i'%(lat.min().values, lat.max().values))

    # corner indices
    corner_indices = ds['connect1']

    if 'color' not in plot_kwargs.keys():
        plot_kwargs['color'] = 'black'
        
    for i in range(corner_indices.shape[0]):
        # get element corners
        lon_corners = lon[corner_indices[i,:] - 1]
        lat_corners = lat[corner_indices[i,:] - 1]

        # map corners to projection
        xx, yy = lon_corners.values, lat_corners.values
        xx = tuple(xx) + (xx[0],)
        yy = tuple(yy) + (yy[0],)

        # draw element boundaries as great circles
        for j in range(len(xx)-1):
            pl = ax.plot(
                [xx[j], xx[j+1]], [yy[j], yy[j+1]], 
                transform=crs.Geodetic(), **plot_kwargs
            )

    return pl, ax


def draw_element_patch(ds, i):
    from matplotlib import pyplot
    ax = pyplot.gca()
    
    # cartesian coordinates (radians?)
    x = ds['coord'][0,:].squeeze()
    y = ds['coord'][1,:].squeeze()
    z = ds['coord'][2,:].squeeze()

    # convert to longitude and latitude
    lon = numpy.arctan2(y, x) * 180.0 / numpy.pi
    lat = numpy.arcsin(z) * 180.0 / numpy.pi

    # get corners
    lon_corners = lon[corner_indices[i,:] - 1]
    lat_corners = lat[corner_indices[i,:] - 1]
    
    xx, yy = lon_corners.values, lat_corners.values
    xx = tuple(xx) + (xx[0],)
    yy = tuple(yy) + (yy[0],)

    # draw element as a patch from path
    vertices = [(x, y) for (x, y) in zip(xx, yy)]
    
    # path instructions
    pathCodes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]
    
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    path = Path(vertices, pathCodes)
    patch = PathPatch(path, **kwargs)
    ax.add_patch(patch, transform=crs.Geodetic(), **kwargs)
    

def main(inputfile, outputfile, **plot_kwargs):
    """
    Plot exodus-format grid files generated by SQuadGen.
    """

    import matplotlib; matplotlib.use('Agg')
    from matplotlib import pyplot
    from xarray import open_dataset

    with open_dataset(inputfile) as ds:
        figure = pyplot.figure()
        ax = plot_exodus(ds, **plot_kwargs)
        figure.savefig(outputfile, bbox_inches='tight')


if __name__ == '__main__':
    import plac; plac.call(main)