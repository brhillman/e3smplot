#!/usr/bin/env python3
import plac, numpy
from matplotlib import pyplot
from matplotlib.tri import Triangulation
from cartopy import crs
from cartopy.util import add_cyclic_point
from time import perf_counter
from scipy.interpolate import griddata
from plot_maps import plot_map
from plot_utils import open_files, get_data, fix_longitudes, area_average


def main(varname, outputfile, testfile, cntlfile, **kwargs):
    #
    # Open datasets
    #
    datasets = [open_files(f) for f in (testfile, cntlfile)]
    #
    # Read selected data from file
    # TODO: set case names
    #
    data_arrays = [get_data(ds, varname) for ds in datasets]
    lons = [get_data(ds, 'lon') for ds in datasets]
    lats = [get_data(ds, 'lat') for ds in datasets]
    #
    # Area needed for weighted average calculation
    #
    area_arrays = [get_data(ds, 'area') for ds in datasets]
    #
    # Reduce data
    # TODO: apply vertical reduction here
    #
    data_arrays = [d.mean(dim='time', keep_attrs=True) for d in data_arrays]
    #
    # Compute differences
    #
    data_arrays.append(data_arrays[0] - data_arrays[1])
    area_arrays.append(area_arrays[0])
    lons.append(lons[0])
    lats.append(lats[0])
    data_arrays[-1].attrs = data_arrays[0].attrs
    data_arrays[-1].attrs['description'] = 'test minus control'
    #
    # Pop labels out of kwargs
    #
    if 'labels' in kwargs.keys():
        labels = kwargs.pop('labels').split(',')
    else:
        labels = None
    #
    # Make figure
    #
    figure, axes = pyplot.subplots(1, 3, figsize=(15, 5), subplot_kw=dict(projection=crs.PlateCarree(central_longitude=180)))
    cmaps = ['viridis', 'viridis', 'RdBu_r']
    vmin = min([d.min().values for d in data_arrays[:-1]])
    vmax = max([d.max().values for d in data_arrays[:-1]])
    vmins = [vmin, vmin, -abs(data_arrays[-1].max().values)]
    vmaxs = [vmax, vmax,  abs(data_arrays[-1].max().values)]
    p = [plot_map(lons[i], lats[i], data_arrays[i], axes=axes[i], cmap=cmaps[i], vmin=vmins[i], vmax=vmaxs[i], **kwargs) for i in range(len(data_arrays))]
    #
    # Annotate maps
    #
    for i in range(len(data_arrays)):
        data_mean = area_average(data_arrays[i], area_arrays[i])
        if i < len(datasets):
            if labels is not None:
                label = labels[i]
            else:
                label = datasets[i].attrs['case']
        else:
            label = 'Difference'
        axes[i].set_title(
            f'{label}\nmin = {data_arrays[i].min().values:.2f}; max = {data_arrays[i].max().values:.2f}; mean = {data_mean.values:.2f}'
        )
    #
    # Save figure
    #
    figure.savefig(outputfile, bbox_inches='tight')
    pyplot.close()


if __name__ == '__main__':
    import plac; plac.call(main)
