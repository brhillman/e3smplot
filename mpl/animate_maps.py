#!/usr/bin/env python3

import plac, os, imageio, sys, numpy, functools
from matplotlib import pyplot
from matplotlib.tri import Triangulation
from cartopy import crs
from cartopy.util import add_cyclic_point
from xarray import open_mfdataset
from time import perf_counter
from scipy.interpolate import griddata


def remove_dims(ds):
    return ds.drop_vars(['P3_input_dim', 'P3_output_dim'], errors='ignore')


def open_files(*inputfiles):
    print('Found %i files'%len(inputfiles))
    return open_mfdataset(
        sorted(inputfiles), combine='by_coords', 
        preprocess=remove_dims,
        chunks={'time': 1}
    )


def memoize(func):
    cache = dict()    
    @functools.wraps(func)
    def memoized_func(*args, **kwargs):
        key = str(args) + str(kwargs)
        if key not in cache:
            cache[key] = func(*args, **kwargs)
        return cache[key]
    return memoized_func


@memoize
def get_triangulation(lon, lat):
    return Triangulation(lon, lat)


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

def plot_frame(dataset, variable_name, frame_name,
               lat_min=None, lat_max=None, lon_min=None, lon_max=None, 
               nlon=360, nlat=180,
               plot_method='regrid',
               **kwargs):

    if 'projection' in kwargs.keys():
        projection=kwargs['projection']
        kwargs.pop('projection')
    else:
        projection=crs.PlateCarree()

    # Select data
    data = get_data(dataset, variable_name).squeeze()
    lon = get_data(dataset, 'lon').squeeze()
    lat = get_data(dataset, 'lat').squeeze()

    # Reduce data
    if 'lev' in data.dims:
        data = data.isel(lev=-1)

    # Fix coordinates
    lon = fix_longitudes(lon)

    # Open figure
    figure = pyplot.figure(figsize=(10, 8))
    axes = pyplot.axes(projection=projection)
    axes.coastlines()
    axes.set_global()

    # Set extent
    if all(v is not None for v in [lat_min, lat_max, lon_min, lon_max]):
        axes.set_extent([float(lon_min), float(lon_max), float(lat_min), float(lat_max)])
        if 'ncol' in data.dims:
            criteria = ((lon >= float(lon_min)) & 
                        (lon <= float(lon_max)) &
                        (lat >= float(lat_min)) &
                        (lat <= float(lat_max)))
                       
            data = data.where(criteria).dropna('ncol')
            lon = lon.where(criteria).dropna('ncol')
            lat = lat.where(criteria).dropna('ncol')

    # Plot data
    if 'ncol' in data.dims:
        if plot_method == 'triangulation':
            # Calculate triangulation
            triangulation = get_triangulation(lon.values, lat.values)

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
    else:
        # Need to add a cyclic point
        #_data, _lon = add_cyclic_point(data.transpose('lat', 'lon').values, coord=lon.values)
        pl = axes.pcolormesh(
            lon, lat, data.transpose('lat', 'lon'),
            transform=crs.PlateCarree(), **kwargs
        )

    # Label plot
    axes.set_title('time = %04i-%02i-%02i %02i:%02i:%02i'%(
        data['time.year'], data['time.month'], data['time.day'],
        data['time.hour'], data['time.minute'], data['time.second']
    ))

    # Add a colorbar
    cb = pyplot.colorbar(
        pl, orientation='horizontal',
        label='%s (%s)'%(data.long_name, data.units),
        shrink=0.8, pad=0.02
    )

    # Save figure
    figure.savefig(frame_name, dpi=100)
    pyplot.close()


# Get a time-varying longitude to be used to rotate map center to mimic earth
# rotation in animations
def rotate_longitude(itime, samples_per_day, start_lon=360):
    central_longitude = (start_lon - itime * 360 / samples_per_day) % 360
    if central_longitude > 180: central_longitude = central_longitude - 360
    return central_longitude
   

def plot_frames(
        dataset, variable_name, 
        rotate=False, samples_per_day=48, 
        **kwargs):

    # Get data range so that all frames will be consistent
    if 'vmin' not in kwargs.keys(): kwargs['vmin'] = get_data(dataset, variable_name).min().values
    if 'vmax' not in kwargs.keys(): kwargs['vmax'] = get_data(dataset, variable_name).max().values

    # Loop over time series and make a plot for each
    frames = []
    print('Looping over %i time indices'%len(dataset.time)); sys.stdout.flush()
    for i in range(len(dataset.time)):
        frame_name = 'tmp_frames/%s.%i.png'%(variable_name, i)
        os.makedirs(os.path.dirname(frame_name), exist_ok=True)
        if rotate:
            central_longitude = rotate_longitude(i, samples_per_day)
            plot_frame(
                dataset.isel(time=i), variable_name, frame_name, 
                projection=crs.Orthographic(central_longitude=central_longitude, central_latitude=20),
                **kwargs
            )
        else:
            plot_frame(
                dataset.isel(time=i), variable_name, frame_name, 
                **kwargs
            )
        frames.append(frame_name)
        update_progress(i+1, len(dataset.time))

    # Return list of frames
    return frames


def animate_frames(outputfile, frames, **kwargs):
    print('Stitching %i frames together...'%(len(frames)), end='')
    sys.stdout.flush()
    images = [imageio.imread(frame) for frame in frames]
    imageio.mimsave(outputfile, images, **kwargs)
    print('Done.'); sys.stdout.flush()


def remove_frames(frames):
    for frame in frames: os.remove(frame)


# update_progress() : Displays or updates a console progress bar
def update_progress(iteration, num_iterations, bar_length=10):

    # Get progress as a fraction and compute size of "block" filled to visually
    # represent fraction completed.
    progress = 1.0 * iteration / num_iterations
    block = int(round(bar_length * progress))

    # Get status; if status < 1 there will be no newline, so we need to add one
    # explicitly on the last iteration.
    if progress >= 1:
        status = "\r\n"
    else:
        status = ""

    # Display appropriate text to build status bar
    text = "\rPercent: [%s] %i of %i (%.2f%%)%s"%( 
        "#"*block + "-"*(bar_length-block), 
        iteration, num_iterations,
        progress*100, status
    )
    sys.stdout.write(text)
    sys.stdout.flush()


def main(outputfile, variable_name, *inputfiles, **kwargs):

    # Open files
    dataset = open_files(*inputfiles)

    # Pull out keyward args
    animate_kw = {}
    for key in ('time_per_frame',):
        if key in kwargs.keys():
            animate_kw[key] = kwargs[key]
            kwargs.pop(key)

    # Plot frames
    if 'rotate' in kwargs.keys() and 'samples_per_day' in kwargs.keys():
        rotate=kwargs['rotate']
        samples_per_day=float(kwargs['samples_per_day'])
        kwargs.pop('rotate')
        kwargs.pop('samples_per_day')
        frames = plot_frames(
            dataset, variable_name, rotate=rotate, samples_per_day=samples_per_day, **kwargs
        )
    else:
        frames = plot_frames(dataset, variable_name, **kwargs)

    # Stitch together frames into single animation
    animate_frames(outputfile, frames, **animate_kw)

    # Clean up
    remove_frames(frames)


if __name__ == '__main__':
    plac.call(main)
