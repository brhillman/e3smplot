#!/usr/bin/env python3
import plac, numpy, xarray
from matplotlib import pyplot
from matplotlib.tri import Triangulation
from cartopy import crs
from cartopy.util import add_cyclic_point
from time import perf_counter
from scipy.interpolate import griddata
from e3smplot.e3sm_utils import infer_grid_file, infer_grid_coords, get_data


def open_files(*inputfiles):
    return open_mfdataset(
        sorted(inputfiles), combine='by_coords', 
        drop_variables=['P3_input_dim', 'P3_output_dim'],
        chunks={'time': 1}
    )


def fix_longitudes(lon):
    #return lon.assign_coords(lon=((lon + 180) % 360) - 180) #numpy.where(lon > 180, lon - 360, lon))
    return ((lon + 180) % 360) - 180 #numpy.where(lon > 180, lon - 360, lon))


def plot_map(lon, lat, data, axes=None,
        plot_method='regrid', nlon=360, nlat=180, draw_colorbar=True,
        cb_kwargs={'pad': 0.02, 'shrink': 0.8, 'orientation': 'horizontal'}, **kwargs):

    # Get current axes
    if axes is None:
        axes = pyplot.gca()

    # Draw coastlines on map, only works on GeoAxes
    try:
        axes.coastlines(color='white')
        axes.set_global()
    except:
        pass

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
                **kwargs
            )
        elif plot_method == 'regrid':
            xi = numpy.linspace(-180, 180, int(nlon))
            yi = numpy.linspace(-90, 90, int(nlat))
            data_regridded = griddata((lon, lat), data, (xi[None,:], yi[:,None]), method='nearest')
            pl = axes.pcolormesh(xi, yi, data_regridded, **kwargs)
        else:
            raise ValueError('method %s not known; please choose one of triangulation or regrid'%method)
    elif len(data.shape) == 2:
        # Need to add a cyclic point
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
    if 'time' in data.dims:
        axes.set_title('time = %04i-%02i-%02i %02i:%02i:%02i'%(
            data['time.year'], data['time.month'], data['time.day'],
            data['time.hour'], data['time.minute'], data['time.second']
        ))

    # Add a colorbar
    if draw_colorbar:
        cb = pyplot.colorbar(
            pl, ax=axes,
            label=f'{data.long_name} ({data.units})',
            **cb_kwargs,
        )
    else:
        cb = None
    # Return plot and colorbar
    return pl, cb



def plot_maps(coords, data_arrays, labels, figshape=None, figsize=None, subplot_kw=None, **kwargs):
    if figshape is None: figshape = (len(data_arrays), 1)
    figure, axes = pyplot.subplots(*figshape, figsize=figsize, subplot_kw=subplot_kw)

    # Get mins/maxes over all datasets
    #vmin = min([d.min().values for d in data_arrays])
    #vmax = max([d.max().values for d in data_arrays])

    for i, (c, d, l) in enumerate(zip(coords, data_arrays, labels)):
        ax = figure.add_axes(axes.ravel()[i])
        pl, cb = plot_map(*c, d, **kwargs)
        ax.set_title(l)

    return figure


def open_dataset(f, **kwargs):
    return xarray.open_mfdataset(f, data_vars='minimal', coords='minimal', compat='override')


def main(varnames, outputfile, *inputfiles, **kwargs):

    print(f'Plot {varnames} to {outputfile}...')
    #
    # Open dataset and compute time average if needed
    #
    datasets = [open_dataset(f) for f in inputfiles]
    #
    # Read selected data from file
    #
    data_arrays = [get_data(ds, varname) for varname, ds in zip(varnames, datasets)]
    lons, lats = zip(*[(get_data(ds, 'longitude'), get_data(ds, 'latitude')) for ds in datasets])
    #
    # Reduce data
    #
    data_arrays = [da.mean(dim='time', keep_attrs=True) for da in data_arrays]
    #
    #
    # Make figure
    #
    figure, axes = pyplot.subplots(len(inputfiles), 1, figsize=(8, 11), subplot_kw=dict(projection=crs.PlateCarree(central_longitude=180)))
    for icase, (lon, lat, data) in enumerate(zip(lons, lats, data_arrays)):
        ax = figure.add_axes(axes[icase])
        pl, cb = plot_map(lon, lat, data, transform=crs.PlateCarree(), **kwargs)
    #
    # Save figure
    #
    figure.savefig(outputfile, bbox_inches='tight')


if __name__ == '__main__':
    import plac; plac.call(main)
