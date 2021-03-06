#!/usr/bin/env python3
import plac, numpy
from matplotlib import pyplot
from matplotlib.tri import Triangulation
from cartopy import crs
from cartopy.util import add_cyclic_point
from xarray import open_mfdataset
from time import perf_counter
from scipy.interpolate import griddata


def open_files(*inputfiles):
    return open_mfdataset(
        sorted(inputfiles), combine='by_coords', 
        drop_variables=['P3_input_dim', 'P3_output_dim'],
        chunks={'time': 1}
    )


def fix_longitudes(lon):
    lon.values = numpy.where(lon > 180, lon - 360, lon)
    return lon


def get_data(ds, variable_name):
    if variable_name in ds.variables.keys():
        data = ds[variable_name]
    elif variable_name == 'PRECT':
        precc = get_data(ds, 'PRECC')
        precl = get_data(ds, 'PRECL')
        data = precc + precl
        data.attrs = precc.attrs
        data.attrs['long_name'] = 'Total precipitation rate'
    else:
        raise NameError('%s not found in dataset'%variable_name)

    # Adjust units
    if variable_name in ('PRECC', 'PRECL', 'PRECT'):
        if data.attrs['units'].lower() == 'm/s':
            attrs = data.attrs
            data = 60 * 60 * 24 * 1e3 * data
            data.attrs = attrs
            data.attrs['units'] = 'mm/day'
    return data


def plot_map(lon, lat, data, axes=None, plot_method='triangulation', nlon=360, nlat=180, **kwargs):

    # Get current axes
    if axes is None:
        axes = pyplot.gca()

    # Draw coastlines on map
    axes.coastlines(color='white')

    # Fix longitudes so we are always working in -180 to 180 space
    lon = fix_longitudes(lon)

    # Plot data
    if 'ncol' in data.dims:
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
    elif ('lon' in data.dims) and ('lat' in data.dims):
        # Need to add a cyclic point
        #_data, _lon = add_cyclic_point(data.transpose('lat', 'lon').values, coord=lon.values)
        pl = axes.pcolormesh(
            lon, lat, data.transpose('lat', 'lon'),
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
    cb = pyplot.colorbar(
        pl, ax=axes, orientation='horizontal',
        label='%s (%s)'%(data.long_name, data.units),
        shrink=0.8, pad=0.02
    )

    # Return plot and colorbar
    return pl, cb


def main(varname, outputfile, *inputfiles, **kwargs):

    print(f'Plot {varname} to {outputfile}...')
    #
    # Open dataset
    #
    ds = open_files(*inputfiles)
    #
    # Read selected data from file
    #
    data = get_data(ds, varname)
    lon = get_data(ds, 'lon')
    lat = get_data(ds, 'lat')
    #
    # Reduce data
    #
    data = data.mean(dim='time', keep_attrs=True)
    # TODO: apply vertical reduction here
    #
    # Make figure
    #
    figure, axes = pyplot.subplots(1, 1, subplot_kw=dict(projection=crs.PlateCarree()))
    pl, cb = plot_map(lon, lat, data, **kwargs)
    #
    # Save figure
    #
    figure.savefig(outputfile, bbox_inches='tight')


if __name__ == '__main__':
    import plac; plac.call(main)
