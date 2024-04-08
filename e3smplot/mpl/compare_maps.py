#!/usr/bin/env python3
import plac, numpy, warnings, xarray
from matplotlib import pyplot
from matplotlib.tri import Triangulation
from cartopy import crs
from cartopy.util import add_cyclic_point
from time import perf_counter
from scipy.interpolate import griddata
from e3smplot.mpl.plot_maps import plot_map
from e3smplot.plot_utils import area_average
from e3smplot.e3sm_utils import get_data, get_area_weights
from e3smplot.utils import apply_map, myprint

def compare_maps(coords, data_arrays, labels, figsize=None, cb_kwargs=None, **kwargs):
    
    figure, axes = pyplot.subplots(len(data_arrays)+1, 1, figsize=figsize)
    vmin = min([d.min().values for d in data_arrays])
    vmax = max([d.max().values for d in data_arrays])
    for i, ((x, y), d, l) in enumerate(zip(coords, data_arrays, labels)):
        
        # Plot full fields
        ax = figure.add_axes(axes.ravel()[i])
        if i == 1:
            cb_label='%s (%s)'%(d.long_name, d.units)
            _cb_kwargs = {'label': cb_label, **cb_kwargs}
            pl, *__ = plot_map(x, y, d, vmin=vmin, vmax=vmax, cb_kwargs=_cb_kwargs, **kwargs)
        else:
            pl, *__ = plot_map(x, y, d, vmin=vmin, vmax=vmax, cb_kwargs=cb_kwargs, **kwargs)

        ax.set_title(l)
        
        # Plot diffs
        if i == 0:
            d_cntl = d.copy(deep=True)
        else:
            ax = figure.add_axes(axes.ravel()[-1])
            d_diff = d - d_cntl
            d_diff.attrs = d.attrs
            #d_diff.attrs['long_name'] = 'Difference'
            dmax = abs(d_diff).max().values
            pl = plot_map(x, y, d_diff, cmap='bwr', vmin=-dmax, vmax=dmax,
                    cb_kwargs=cb_kwargs, **kwargs)
            ax.set_title(f'Difference')
    #figure.suptitle(f'{d.long_name} ({d.units})', x=0, y=0, ha='left', va='bottom')

    return figure

def open_dataset(*inputfiles):
    return xarray.open_mfdataset(inputfiles, use_cftime=True, data_vars='minimal', coords='minimal', compat='override')

def main(varnames, outputfile, testfiles, cntlfiles,
         t_index=None, t1=None, t2=None, mapfiles=None, percentile=5, verbose=False, **kwargs):
    #
    # Open datasets if needed (may also pass datasets directly rather than filenames)
    #
    if verbose: myprint('Open datasets...')
    datasets = [f if isinstance(f, xarray.Dataset) else open_dataset(*f) for f in (testfiles, cntlfiles)]
    #
    # Compute time average
    #
    if t_index is None:
        if verbose: myprint('Compute time averages...')
        datasets = [ds.mean(dim='time', keep_attrs=True) for ds in datasets]
    else:
        if verbose: myprint(f'Select time index {t_index}...')
        datasets = [ds.isel(time=t_index) for ds in datasets]
    #
    # Read selected data from file
    # TODO: set case names
    #
    if verbose: myprint('Get data...')
    data_arrays = [get_data(ds, varname) for ds, varname in zip(datasets,
                                                                varnames)]
    lons = [get_data(ds, 'longitude') for ds in datasets]
    lats = [get_data(ds, 'latitude') for ds in datasets]
    #
    # Area needed for weighted average calculation
    #
    if verbose: myprint('Get area weights...')
    area_arrays = [get_area_weights(ds) for ds in datasets]
    #
    # Remap if needed
    #
    if verbose: myprint('Remap to lat/lon grid if needed...')
    if mapfiles is not None:
        map_datasets = [xarray.open_dataset(f) if f is not None else None for f in mapfiles]
        area_arrays, lons, lats = zip(*[apply_map(m, f) if f is not
                                   None else (m, x, y) for (m, x, y, f) in
                                   zip(area_arrays, lons, lats, map_datasets)])
        data_arrays, lons, lats = zip(*[apply_map(m, f) if f is not
                                   None else (m, x, y) for (m, x, y, f) in
                                   zip(data_arrays, lons, lats, map_datasets)])
    #
    # Compute differences
    #
    if verbose: myprint('Compute differences...')
    data_arrays = (*data_arrays, data_arrays[0].copy(deep=True))
    data_arrays[-1].data = data_arrays[0].data - data_arrays[1].data
    data_arrays[-1].attrs = data_arrays[0].attrs
    data_arrays[-1].attrs['description'] = 'test minus control'
    area_arrays = (*area_arrays, area_arrays[0])
    lons = (*lons, lons[0])
    lats = (*lats, lats[0])
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
    figure, axes = pyplot.subplots(len(data_arrays), 1, figsize=(8, 5*len(data_arrays)), subplot_kw=dict(projection=crs.PlateCarree(central_longitude=180)))
    cmaps = ['viridis', 'viridis', 'RdBu_r']
    vmin = min([numpy.nanpercentile(da.values, percentile) for da in data_arrays[:-1]])
    vmax = max([numpy.nanpercentile(da.values, 100-percentile) for da in data_arrays[:-1]])
    vmins = [vmin, vmin, -abs(data_arrays[-1].max().values)]
    vmaxs = [vmax, vmax,  abs(data_arrays[-1].max().values)]
    plots = [
        plot_map(lons[i], lats[i], data_arrays[i], axes=axes[i], cmap=cmaps[i], vmin=vmins[i], vmax=vmaxs[i], **kwargs)
        for i in range(len(data_arrays))
    ]
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
    pyplot.tight_layout()
    figure.savefig(outputfile, bbox_inches='tight')
    pyplot.close()


if __name__ == '__main__':
    import plac; plac.call(main)
