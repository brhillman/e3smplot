#!/usr/bin/env python3

def plot_centers(dataset, ax=None, **kwargs):
    from matplotlib import pyplot
    from cartopy import crs
    if ax is None:
        ax = pyplot.subplot(projection=crs.PlateCarree())

    x_centers = dataset['grid_center_lon']
    y_centers = dataset['grid_center_lat']
    for ielement in range(dataset.dims['grid_size']):
        # plot centers
        pl = ax.plot(
            x_centers[ielement], y_centers[ielement],
            transform=crs.PlateCarree(), 
            **kwargs
        )

    return ax


def plot_corners(dataset, ax=None, **kwargs):
    from matplotlib import pyplot
    from cartopy import crs
    if ax is None:
        ax = pyplot.subplot(projection=crs.PlateCarree())

    for ielement in range(dataset.dims['grid_size']):
        # grid corners
        x_corners = dataset['grid_corner_lon'][ielement, :]
        y_corners = dataset['grid_corner_lat'][ielement, :]

        # close polygons by appending first entry to the end
        from xarray import concat
        x_corners = concat([x_corners, x_corners[0]], dim='grid_corners')
        y_corners = concat([y_corners, y_corners[0]], dim='grid_corners')

        # draw corners and polygons
        pl = ax.plot(
            x_corners, y_corners, transform=crs.Geodetic(), 
            **kwargs
        )

    return ax


def main(inputfile, outputfile):
    from xarray import open_dataset
    import matplotlib; matplotlib.use('Agg')
    from matplotlib import pyplot
    with open_dataset(inputfile) as ds:
        figure = pyplot.figure()
        ax = plot_centers(dataset, marker='.', color='blue', linestyle='none')
        ax = plot_corners(
            dataset, marker='.', color='red', 
            linestyle='solid', linewidth=0.1
        )

    figure.savefig(outputfile, bbox_inches='tight', dpi=400)


if __name__ == '__main__':
    import plac; plac.call(main)
