#!/usr/bin/env python3
import plac, numpy, warnings, xarray
from matplotlib import pyplot
from matplotlib.tri import Triangulation
from cartopy import crs
from cartopy.util import add_cyclic_point
from time import perf_counter
from scipy.interpolate import griddata
from .plot_maps import plot_map
from .plot_utils import area_average
from e3smplot.e3sm_utils import open_dataset, get_data, get_area_weights
from ..utils import apply_map, myprint


def main(varname, outputfile, testfiles, cntlfiles, t1=None, t2=None, maps=None, percentile=5, verbose=False, **kwargs):
    #
    # Open datasets if needed (may also pass datasets directly rather than filenames)
    #
    if verbose: myprint('Open datasets...')
    datasets = [f if isinstance(f, xarray.Dataset) else open_dataset(*f) for f in (testfiles, cntlfiles)]
    #
    # Subset data
    #
    if verbose: myprint('Subset consistent time periods...')
    if t1 is None: t1 = max([ds.time[0].values for ds in datasets])
    if t2 is None: t2 = min([ds.time[-1].values for ds in datasets])
    if verbose: myprint('Comparing period {} to {}'.format(str(t1), str(t2)))
    datasets = [ds.sel(time=slice(str(t1), str(t2))) for ds in datasets]
    #
    # Compute time average
    #
    if verbose: myprint('Compute time averages...')
    datasets = [ds.mean(dim='time', keep_attrs=True) for ds in datasets]
    #
    # Read selected data from file
    # TODO: set case names
    #
    if verbose: myprint('Get data...')
    data_arrays = [get_data(ds, varname) for ds in datasets]
    lons = [get_data(ds, 'lon') for ds in datasets]
    lats = [get_data(ds, 'lat') for ds in datasets]
    #
    # Area needed for weighted average calculation
    #
    if verbose: myprint('Get area weights...')
    area_arrays = [get_area_weights(ds) for ds in datasets]
    #
    # Remap if needed
    #
    if verbose: myprint('Remap to lat/lon grid if needed...')
    if maps is not None:
        map_datasets = [xarray.open_dataset(f) if f is not None else None for f in maps]
        area_arrays = [apply_map(m, f)[0] if f is not None else m for (m, f) in zip(area_arrays, map_datasets)]
        data_arrays = [apply_map(m, f)[0] if f is not None else m for (m, f) in zip(data_arrays, map_datasets)]

    # redefine lons and lats after remap
    # TODO: this is not unstructured grid-friendly
    lons = [d.lon for d in data_arrays]
    lats = [d.lat for d in data_arrays]
    #
    # Compute differences
    #
    if verbose: myprint('Compute differences...')
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
        labels = (*kwargs.pop('labels'), 'Difference')
    else:
        labels = ('test', 'cntl', 'diff')
    #
    # Make figure
    #
    if verbose: myprint('Make plots...')
    figure, axes = pyplot.subplots(1, 3, figsize=(15, 5), subplot_kw=dict(projection=crs.PlateCarree(central_longitude=180)))
    cmaps = ['viridis', 'viridis', 'RdBu_r']
    vmin = min([numpy.nanpercentile(da.values, percentile) for da in data_arrays[:-1]])
    vmax = max([numpy.nanpercentile(da.values, 100-percentile) for da in data_arrays[:-1]])
    vmins = [vmin, vmin, -abs(data_arrays[-1].max().values)]
    vmaxs = [vmax, vmax,  abs(data_arrays[-1].max().values)]
    plots = [plot_map(lons[i], lats[i], data_arrays[i], axes=axes[i], cmap=cmaps[i], vmin=vmins[i], vmax=vmaxs[i], **kwargs) for i in range(len(data_arrays))]
    #
    # Annotate maps
    #
    if verbose: myprint('Annotate plots...')
    for i in range(len(data_arrays)):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            data_mean = area_average(data_arrays[i], area_arrays[i])
        label = labels[i]
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
