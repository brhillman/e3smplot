#!/usr/bin/env python3

import plac, os, imageio, sys, numpy
from matplotlib import pyplot
from matplotlib.tri import Triangulation
from cartopy import crs
from xarray import open_mfdataset
from time import perf_counter

def open_files(*inputfiles):

    print('Found %i files'%len(inputfiles))
    return open_mfdataset(sorted(inputfiles), chunks={'time': 1})


def memoize(func):
    cache = dict()    
    def memoized_func(*args):
        if args in cache:
            return cache[args]
        result = func(*args)
        cache[args] = result
        return result    
    
    return memoized_func


@memoize
def get_triangulation(lon, lat):
    return Triangulation(lon, lat)


@memoize
def fix_longitudes(lon):
    lon.values = numpy.where(lon > 180, lon - 360, lon)
    return lon


def plot_frame(dataset, variable_name, frame_name, 
               lat_min=None, lat_max=None, lon_min=None, lon_max=None, 
               projection=crs.PlateCarree(central_longitude=180),
               **kwargs):

    # Select data
    data = get_data(dataset, variable_name).squeeze()
    lon = get_data(dataset, 'lon').squeeze()
    lat = get_data(dataset, 'lat').squeeze()

    # Fix coordinates
    lon = fix_longitudes(lon)

    # Open figure
    figure, axes = pyplot.subplots(
        figsize=(10, 8),
        subplot_kw=dict(projection=projection)
    )

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
        # Calculate triangulation
        triangulation = get_triangulation(lon, lat)

        # Plot using triangulation
        pl = axes.tripcolor(
            triangulation, data,
            transform=crs.PlateCarree(), **kwargs
        )
    else:
        pl = axes.pcolormesh(
            lon, lat, data.transpose('lat', 'lon'),
            transform=crs.PlateCarree(), **kwargs
        )

    # Label plot
    axes.set_title('time = %s'%(data['time'].values))
    axes.coastlines()

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
def rotate_longitude(itime, samples_per_day, start_lon=0):
    central_longitude = (start_lon + itime * 360 / samples_per_day) % 360
    if central_longitude > 180: central_longitude = central_longitude - 360
    return central_longitude
   

def plot_frames(dataset, variable_name, **kwargs):

    # Get data range so that all frames will be consistent
    if 'vmin' not in kwargs.keys(): kwargs['vmin'] = get_data(dataset, variable_name).min().values
    if 'vmax' not in kwargs.keys(): kwargs['vmax'] = get_data(dataset, variable_name).max().values

    # Loop over time series and make a plot for each
    frames = []
    print('Looping over %i time indices'%len(dataset.time)); sys.stdout.flush()
    for i in range(len(dataset.time)):
        frame_name = '%s.%i.png'%(variable_name, i)
        plot_frame(dataset.isel(time=i), variable_name, frame_name, **kwargs)
        frames.append(frame_name)
        update_progress(i+1, len(dataset.time))

    # Return list of frames
    return frames


def animate_frames(outputfile, frames):
    print('Stitching %i frames together...'%(len(frames)), end='')
    sys.stdout.flush()
    images = [imageio.imread(frame) for frame in frames]
    imageio.mimsave(outputfile, images)
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

    rotate = False

    # Open files
    dataset = open_files(*inputfiles)

    # Plot frames
    frames = plot_frames(dataset, variable_name, **kwargs)

    # Stitch together frames into single animation
    animate_frames(outputfile, frames)

    # Clean up
    remove_frames(frames)


if __name__ == '__main__':
    plac.call(main)
