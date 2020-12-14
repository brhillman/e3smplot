#!/usr/bin/env python3

from xarray import open_mfdataset
from e3smplot.e3sm_utils import get_data
from e3smplot.utils import update_progress
import sys
import os
import imageio
from e3smplot.pyngl.plot_map import main as plot_map
import gc
import dask
import subprocess

def animate_frames(outputfile, frames, **kwargs):
    print('Stitching %i frames together...'%(len(frames)), end=''); sys.stdout.flush()
    images = [imageio.imread(frame) for frame in frames]
    imageio.mimsave(outputfile, images, **kwargs)
    print('Done.'); sys.stdout.flush()


def main(var_name, animation_name, *inputfiles, **kwargs):
    with dask.config.set(scheduler='single-threaded'):

        # Get some info about dataset by opening once
        with open_mfdataset(sorted(inputfiles), drop_variables=('P3_input_dim', 'P3_output_dim'), chunks={'time': 1}) as ds:
            # Find mins/maxes
            data = get_data(ds, var_name)
            lon = get_data(ds, 'lon')
            lat = get_data(ds, 'lat')
            if 'mpMinLonF' in kwargs.keys(): data = data.where(lon >  float(kwargs['mpMinLonF']))
            if 'mpMaxLonF' in kwargs.keys(): data = data.where(lon <= float(kwargs['mpMaxLonF']))
            if 'mpMinLatF' in kwargs.keys(): data = data.where(lat >  float(kwargs['mpMinLatF']))
            if 'mpMaxLatF' in kwargs.keys(): data = data.where(lat <= float(kwargs['mpMaxLatF']))
            if 'vmin' not in kwargs.keys(): kwargs['vmin'] = data.min().values
            if 'vmax' not in kwargs.keys(): kwargs['vmax'] = data.max().values
            #percentile = 2
            #cmin = min([numpy.nanpercentile(data.values, percentile) for da in data_arrays])
            #cmax = max([numpy.nanpercentile(data.values, 100-percentile) for da in data_arrays])
            # Find length of time dimension to loop over later
            time = ds['time'].copy(deep=True)
            ntime = len(ds.time)

        # Make a bunch of plots, save separately, then stitch together later
        frames = []
        print('Loop over time series...'); sys.stdout.flush()
        for i in range(ntime):
            frame_name = f'tmp_frames/{var_name}.{i}.png'
            kwargs['tiMainString'] = str(time.isel(time=i).values)
            # Try to run this as a subprocess to prevent an excruciatingly
            # frustating memory leak
            #plot_map(var_name, frame_name, *inputfiles, time_index=i, **kwargs)
            args = ["./e3smplot/pyngl/plot_map.py", var_name, frame_name,
                    *inputfiles, f"time_index={i}", 
                    *[f"{k}={v}" for k,v in kwargs.items()]]
            subprocess.run(args)
            # Trim frame
            subprocess.run(f'convert -trim {frame_name} {frame_name}'.split(' '))
            frames.append(frame_name)
            update_progress(i+1, ntime)

    print('Animate frames...'); sys.stdout.flush()
    animate_frames(animation_name, frames)

    print('Remove temporary files...'); sys.stdout.flush()
    for frame in frames: os.remove(frame)

if __name__ == '__main__':
    import plac; plac.call(main)
