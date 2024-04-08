#!/usr/bin/env python3
import plac, numpy
from matplotlib import pyplot
from matplotlib.tri import Triangulation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cartopy import crs
from cartopy.util import add_cyclic_point
import xarray
from time import perf_counter
from scipy.interpolate import griddata
from e3smplot.e3sm_utils import get_data, get_area_weights, area_average


def fix_longitudes(lon):
    return ((lon + 180) % 360) - 180


def plot_map(lon, lat, data, axes=None, plot_method='triangulation', nlon=360, nlat=180, title=None, **kwargs):

    # Get current axes
    if axes is None:
        axes = pyplot.gca()

    # Draw coastlines on map
    axes.coastlines(color='white', linewidth=0.5)

    # Fix longitudes so we are always working in -180 to 180 space
    lon = fix_longitudes(lon)

    # Plot data
    if len(data.shape) == 1:
        if plot_method == 'triangulation':
            # Calculate triangulation
            triangulation = Triangulation(lon.values, lat.values)

            # Plot using triangulation
            pl = axes.tripcolor(
                triangulation, data,
                transform=crs.PlateCarree(), **kwargs
            )
        elif plot_method == 'regrid':
            xi = numpy.linspace(-180, 180, int(nlon))
            yi = numpy.linspace(-90, 90, int(nlat))
            data_regridded = griddata((lon, lat), data, (xi[None,:], yi[:,None]), method='nearest')
            pl = axes.pcolormesh(xi, yi, data_regridded, transform=crs.PlateCarree(), **kwargs)
        else:
            raise ValueError('method %s not known'%method)
    #elif ('lon' in data.dims) and ('lat' in data.dims):
    elif len(data.shape) == 2:
        # Need to add a cyclic point?
        #_data, _lon = add_cyclic_point(data.transpose('lat', 'lon').values, coord=lon.values)

        # Figure out dim indices for lon, lat to transpose
        transpose_dims = data.dims[data.shape.index(len(lat))], data.dims[data.shape.index(len(lon))]
        pl = axes.pcolormesh(
            lon, lat, data.transpose(*transpose_dims), 
            transform=crs.PlateCarree(), **kwargs
        )
    else:
        raise ValueError('Dimensions invalid.')

    # Label plot
    if 'time' in data.dims and title is None:
        axes.set_title('time = %04i-%02i-%02i %02i:%02i:%02i'%(
            data['time.year'], data['time.month'], data['time.day'],
            data['time.hour'], data['time.minute'], data['time.second']
        ))
    else:
        axes.set_title(title)

    # Add a colorbar
    divider = make_axes_locatable(axes)
    cax = divider.append_axes("bottom", size="5%", pad=0.1, axes_class=pyplot.Axes)
    cb = pyplot.colorbar(
        pl, cax=cax, orientation='horizontal',
        label='%s (%s)'%(data.long_name, data.units),
    )

    # Return plot and colorbar
    return pl, cb


def main(varname, outputfile, *inputfiles, **kwargs):

    print(f'Plot {varname} to {outputfile}...')

    # Open dataset
    ds = xarray.open_mfdataset(inputfiles, data_vars='minimal', coords='minimal', compat='override')

    # Read selected data from file
    data = get_data(ds, varname)
    lon  = get_data(ds, 'lon')
    lat  = get_data(ds, 'lat')

    # Reduce data; TODO: apply vertical reduction?
    data = data.mean(dim='time', keep_attrs=True)

    # Compute area average
    weights = get_area_weights(ds)
    data_mean = area_average(data, weights)
    plot_title = f'mean = {data_mean.values:.2f}; min = {data.min().values:.2f}; max = {data.max().values:.2f}'

    # Make figure
    figure, axes = pyplot.subplots(1, 1, subplot_kw=dict(projection=crs.PlateCarree()))
    pl, cb = plot_map(lon, lat, data, title=plot_title, **kwargs)
    figure.savefig(outputfile, bbox_inches='tight')

    # Clean up
    pyplot.close(figure)
    ds.close()


if __name__ == '__main__':
    import plac; plac.call(main)
