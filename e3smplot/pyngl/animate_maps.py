#!/usr/bin/env python3

import plac, os, imageio, sys, numpy, functools
from xarray import open_mfdataset
from time import perf_counter
from e3smplot.e3sm_utils import get_data
from e3smplot.pyngl.plot_map import plot_map
from e3smplot.utils import nice_cntr_levels
import ngl


def remove_dims(ds):
    return ds.drop_vars(['P3_input_dim', 'P3_output_dim'], errors='ignore')


def open_files(*inputfiles):
    print('Found %i files'%len(inputfiles))
    return open_mfdataset(
        sorted(inputfiles), combine='by_coords', 
        drop_variables=('P3_input_dim', 'P3_output_dim'),
        chunks={'time': 1}
    )


def plot_frame(lon, lat, data, frame_name,
               lat_min=None, lat_max=None, lon_min=None, lon_max=None, 
               **kwargs):

    # Open figure
    wks = ngl.open_wks('png', os.path.splitext(frame_name)[0])

    # Plot data
    pl = plot_map(
        wks, lon, lat, data,
        #tiMainString=f'{str(data.time.values)}',
        lbOrientation='horizontal',
        #lbTitleString='%s (%s)'%(data.long_name, data.units),
        cnFillMode='RasterFill',
        **kwargs
    )

    # Save figure
    ngl.draw(pl)
    ngl.destroy(wks)

    return frame_name


# Get a time-varying longitude to be used to rotate map center to mimic earth
# rotation in animations
def rotate_longitude(itime, samples_per_day, start_lon=360):
    central_longitude = (start_lon - itime * 360 / samples_per_day) % 360
    if central_longitude > 180: central_longitude = central_longitude - 360
    return central_longitude
   
def get_contour_levels(data_arrays, percentile=2, **kwargs):
    # Try to get robust contour intervals
    if False:
        cmin = min([numpy.nanpercentile(da.values, percentile) for da in data_arrays])
        cmax = max([numpy.nanpercentile(da.values, 100-percentile) for da in data_arrays])
        if 'aboutZero' not in kwargs.keys(): kwargs['aboutZero'] = (cmin < 0 and cmax > 0)
        *__, clevels = nice_cntr_levels(cmin, cmax, returnLevels=True, max_steps=13, **kwargs)
    # Fall back to just using the min and max to set the limits
    else:
        cmin = min([da.min().values for da in data_arrays])
        cmax = max([da.max().values for da in data_arrays])
        if 'aboutZero' not in kwargs.keys(): kwargs['aboutZero'] = (cmin < 0 and cmax > 0)
        *__, clevels = nice_cntr_levels(cmin, cmax, returnLevels=True, max_steps=13, **kwargs)
    return clevels

def plot_frames(
        dataset, variable_name, 
        **kwargs):

    lon = get_data(dataset, 'lon')
    lat = get_data(dataset, 'lat')
    data = get_data(dataset, variable_name) 
    data = data.where(
        (lon > float(kwargs['mpMinLonF'])) & (lon <= float(kwargs['mpMaxLonF'])) &
        (lat > float(kwargs['mpMinLatF'])) & (lat <= float(kwargs['mpMaxLatF']))
    )
    print('Min of data:', data.min().values)

    # Get data range so that all frames will be consistent
    print('Get mins/maxes over entire dataset...'); sys.stdout.flush()
    kwargs['cnLevels'] = get_contour_levels(data)
    kwargs['cnLevelSelectionMode'] = 'ExplicitLevels'

    # Loop over time series and make a plot for each
    os.makedirs('tmp_frames', exist_ok=True)
    print('Looping over %i time indices'%len(dataset.time)); sys.stdout.flush()
    frames = []
    for i in range(len(dataset.time)):
        frame_name = f'tmp_frames/{variable_name}.{i}.png'
        pl = plot_frame(lon.isel(time=i).values, lat.isel(time=i).values, data.isel(time=i).values, frame_name, **kwargs)
        frames.append(frame_name)
        update_progress(i, len(dataset.time), bar_length=10)

    # Return list of frames
    return frames


def animate_frames(outputfile, frames, **kwargs):
    print('Stitching %i frames together...'%(len(frames)), end='')
    sys.stdout.flush()
    print(frames)
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
    import xarray
    import gc
    # Open files
    #dataset = open_files(sorted(inputfiles))
    with xarray.open_mfdataset(
        sorted(inputfiles), drop_variables=('P3_input_dim',
        'P3_output_dim'), chunks={'time': 1}
    ) as ds:

        # Get data
        data = get_data(ds, variable_name)
        lon = get_data(ds, 'lon')
        lat = get_data(ds, 'lat')

        # Subset
        data = data.where(
            (lon > float(kwargs['mpMinLonF'])) & (lon <= float(kwargs['mpMaxLonF'])) &
            (lat > float(kwargs['mpMinLatF'])) & (lat <= float(kwargs['mpMaxLatF']))
        )

        # Get mins and maxes
        minval = data.min().values
        maxval = data.max().values
        print(f'Data range: {minval} to {maxval}')

        # Get contour levels
        *__, clevels = nice_cntr_levels(minval, maxval, returnLevels=True, max_steps=13)
        kwargs['cnLevels'] = clevels #get_contour_levels(data)
        kwargs['cnLevelSelectionMode'] = 'ExplicitLevels'
        print('clevels: ', kwargs['cnLevels'])

        ntime = len(ds.time)

    # Loop
    print('Looping over %i time indices'%len(ds.time)); sys.stdout.flush()
    frames = []
    for i in range(ntime):
        gc.collect()
        with xarray.open_mfdataset(
            sorted(inputfiles), drop_variables=('P3_input_dim',
            'P3_output_dim'), chunks={'time': 1}
        ) as ds:
            data = get_data(ds, variable_name).isel(time=i).values
            lon = get_data(ds, 'lon').isel(time=i).values
            lat = get_data(ds, 'lat').isel(time=i).values
            plot_name = f'tmp_frames/{variable_name}.{i}.png'
            frames.append(plot_frame(lon, lat, data, plot_name, **kwargs))
            update_progress(i, ntime, bar_length=10)
            del data
            del lon
            del lat

    # Pull out keyward args
    #animate_kw = {}
    #for key in ('time_per_frame',):
    #    if key in kwargs.keys():
    #        animate_kw[key] = kwargs.pop(key)

    # Plot frames
    #frames = plot_frames(dataset, variable_name, **kwargs)

    # Stitch together frames into single animation
    animate_frames(outputfile, frames, **animate_kw)

    # Clean up
    remove_frames(frames)


if __name__ == '__main__':
    plac.call(main)
