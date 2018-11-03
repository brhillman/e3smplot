#!/usr/bin/env python3
from matplotlib import pyplot
from cartopy import crs

def plot_field(dataset, field, plot_type='map', **kwargs):

    # Get data
    from e3sm_utils import get_data
    data = get_data(dataset, field)

    if plot_type == 'map':
        lon = get_data(dataset, 'lon')
        lat = get_data(dataset, 'lat')
        pl = plot_map(lon, lat, data, **kwargs)

    return pl
    

def plot_profile(data, levels, ax=None, **kwargs):
    
    # Get axes to plot in if not passed
    if ax is None: ax = pyplot.gca()
        
    # Make plot
    pl = ax.plot(data, levels, **kwargs)
    
    # Label plot
    ax.set_xlabel('%s (%s)'%(data.long_name, data.units))
    ax.set_ylabel('%s (%s)'%(levels.long_name, levels.units))
    
    return pl

def compare_profiles(data_arrays, level_arrays, labels=None, plot_differences=False, **kwargs):
    figure = pyplot.figure()
    ax = figure.add_subplot(111)
    
    for icase, (data, levels, label) in enumerate(zip(data_arrays, level_arrays, labels)):
        
        # Figure out levels to plot against
        #if levels is None:
        #    from xarray import DataArray
        #    levels = DataArray(range(len(data)), attrs={'long_name': 'Level index', 'units': 1})
        
        # Make plot for this case
        pl = plot_profile(data, levels, label=label, **kwargs)
        
        # If adding a difference line, do all that...
        if plot_differences is True:
            if icase == 0:
                data_cntl = data.copy(deep=True)
            else:
                ax_diff = ax.twiny()
                data_diff = data - data_cntl
                pl = ax_diff.plot(data_diff, levels, label='difference',
                                  color='0.5', linestyle='solid')

                # add a vertical line
                pl_diff = ax_diff.plot([0, 0], [levels[0], levels[-1]], color='0.5', linestyle='dashed')

                # set limits
                if ax_diff is not ax:
                    x2 = abs(data_diff).max()
                    x1 = -x2
                    ax_diff.set_xlim(x1, x2)

                    # fix labels for difference axes
                    ax_diff.set_xlabel('Difference', color='0.5')
                    ax_diff.tick_params('x', colors='0.5')
                    ax_diff.ticklabel_format(axis='x', style='sci', scilimits=(-3, 3))

    # Fix axes ticklabels
    ax.ticklabel_format(style='sci', axis='x', scilimits=(-2,2))    

    # Fix axes
    if levels.units in ('hPa', 'mb', 'Pa'):
        y1, y2 = ax.get_ylim()
        y1, y2 = max(y1, y2), min(y1, y2)
        ax.set_ylim(y1, y2)
        if plot_differences: ax_diff.set_ylim(y1, y2)
    
    return figure


# define a function to set a circular plot area
# We will make polar map plots using a sterographic projection. 
# We would like for these to be round, but we need to define a 
# function to fix the axes bounds when using the cartopy map 
# projection library to make this happen.
def cartopy_circular(ax):
    import matplotlib.path as mpath
    import matplotlib.pyplot as plt
    import numpy as np

    # Compute a circle in axes coordinates, which we can use as a boundary
    # for the map. We can pan/zoom as much as we like - the boundary will be
    # permanently circular.
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)
    
    



def compare_maps_diff(lon, lat, data1, data2, labels=('case 1', 'case 2'), **kwargs):
    
    from matplotlib import pyplot
    from cartopy import crs
    
    # setup figure
    figure, axes = pyplot.subplots(
        1, 3, figsize=(15, 3), 
        subplot_kw=dict(projection=crs.PlateCarree())
    )
    
    # Find data limits
    vmin = min(data1.min().values, data2.min().values)
    vmax = max(data1.max().values, data2.max().values)
    
    # Loop over cases and plot
    for icase, (data, case) in enumerate(zip((data1, data2), labels)):
        ax = figure.add_axes(axes[icase])
        pl, cb = plot_map(
            lon, lat, data, vmin=vmin, vmax=vmax, 
            transform=crs.PlateCarree(), 
            **kwargs
        )
        
        # Label plot
        ax.set_title('%s (min = %.3f, max = %.3f)'%(
                     case, data.min().values, data.max().values))
    
        # Plot differences
        if icase == 0:
            data_cntl = data.copy(deep=True)   
        else:
            # calculate differences
            data_diff = data - data_cntl
            data_diff.attrs = data.attrs
            
            # get data limits
            vmax = abs(data_diff).max().values
            vmin = -vmax
            
            # plot
            ax = figure.add_axes(axes[-1])
            pl, cb = plot_map(
                lon, lat, data_diff, 
                cmap='RdBu_r', vmin=vmin, vmax=vmax,
                transform=crs.PlateCarree(),
                **kwargs
            )
            ax.set_title('Difference (min = %.3f, max = %.3f)'%(
                data_diff.min().values, data_diff.max().values
            ))
            
    # Return figure object
    return figure


def plot_map_native(lon_corners, lat_corners, data, **kwargs):
    from matplotlib import pyplot, patches, path    
    import numpy
    
    # Fix longitude coordinates; we need longitudes to run from
    # -180 to 180, not from 0 to 360
    lon_corners.values = numpy.where(lon_corners > 180, lon_corners - 360, lon_corners)

    # Loop over GLL nodes and create paths that describe the boundaries
    # of each node; collect these into a list to plot all at once
    path_list = []
    for icol in range(lon_corners.shape[0]):
       
        # Get corners for this node
        x_corners = lon_corners[icol,:].values
        y_corners = lat_corners[icol,:].values
        
        # Repeat first vertex at end of array to close the path
        x_corners = numpy.append(x_corners, x_corners[0])
        y_corners = numpy.append(y_corners, y_corners[0])
            
        # Create paths connecting the corners and append to list
        vertices = numpy.column_stack([x_corners, y_corners])
        path_list.append(path.Path(vertices, closed=True))
        
    # Plot collection of patches
    from matplotlib.collections import PathCollection
    collection = PathCollection(path_list, transform=crs.Geodetic(), **kwargs)
    collection.set_array(data)
    
    ax = pyplot.gca()
    pl = ax.add_collection(collection)
    
    return pl


def plot_map(lon, lat, data, lon_corners=None, lat_corners=None, **kwargs):

    # Fix longitudes
    import numpy
    new_lon = numpy.where(lon > 180, lon - 360, lon)
    
    # Setup plot axes
    ax = pyplot.gca()
    ax.set_global()
    ax.coastlines()

    # Make plot
    if all([v is not None for v in (lon_corners, lat_corners)]):
        pl = plot_map_native(lon_corners, lat_corners, data, **kwargs)
    elif 'ncol' in data.dims:
        pl = ax.tripcolor(new_lon, lat, data, **kwargs)
    else:
        pl = ax.pcolormesh(new_lon, lat, data.transpose('lat', 'lon'), **kwargs)
    
    # Return plot handle
    return pl


def compare_maps(data_arrays, labels=None, 
                 ncols=None, nrows=None,
                 lat_bounds=None,
                 projection=crs.PlateCarree(), vmin=None, vmax=None, **kwargs):
    
    # TODO: allow for difference plots
    
    # Figure out size of figure
    if ncols is None: ncols = len(data_arrays)
    if nrows is None: nrows = 1
        
    # Find common mins and maxes
    if vmin is None: vmin = min([data.min().values for data in data_arrays])
    if vmax is None: vmax = max([data.max().values for data in data_arrays])
    
    # Open figure
    figure, axes = pyplot.subplots(nrows, ncols, subplot_kw=dict(projection=projection))
    
    # Loop and plot
    for icase, data in enumerate(data_arrays):
        
        # Make plot of this data array
        ax = figure.add_axes(axes.ravel()[icase])
        
        # fix plot area; not sure why all this is needed, bugs in cartopy?
        if lat_bounds is not None: 
            ax.set_extent([-180, 180, lat_bounds[0], lat_bounds[1]], crs.PlateCarree())
        if projection == crs.NorthPolarStereo(): cartopy_circular(ax)
        
        pl = plot_map(data.lon, data.lat, data, 
                      transform=crs.PlateCarree(), vmin=vmin, vmax=vmax)
        
        # Label plot
        if labels is not None:
            ax.set_title(labels[icase])
        elif 'case' in data.attrs:
            ax.set_title(data.case)
            
    # Common colorbar
    cb = pyplot.colorbar(pl, ax=axes.ravel().tolist(),
                         orientation='horizontal', shrink=0.8, pad=0.02, 
                         label='%s (%s)'%(data.long_name, data.units))

    return figure


def compare_maps_from_ds(datasets, field, plot_diffs=False, **kwargs):
    
    from e3sm_utils import get_data
    
    # Get datarrays
    data_arrays = []
    for ds in datasets:
        da = get_data(ds, field)
    
        # If lat and lon are not packaged with data_arrays as coordinate variables,
        # then package them up as attributes
        if not hasattr(da, 'lon'):
            da['lon'] = get_data(ds, 'lon')
        if not hasattr(da, 'lat'):
            da['lat'] = ds['lat'] #get_data(ds, 'lat')
                        
        # Append this DataArray to list
        data_arrays.append(da)
        
    # Now we can call our compare_maps function that operates on data_arrays        
    if plot_diffs:
        return compare_maps_diff(data_arrays[0].lon, data_arrays[0].lat, *data_arrays, **kwargs)
    else:
        return compare_maps(data_arrays, **kwargs)
    

    
def compare_maps_from_datasets(datasets, field, labels=None, projection=crs.PlateCarree(), 
                               ncols=None, nrows=None, vmin=None, vmax=None,
                               lat_bounds=None, lon_bounds=None, plot_differences=False,
                               **kwargs):
    
    # Figure out size of figure
    if ncols is None: ncols = len(datasets)
    if nrows is None: nrows = 1
        
    # Find common mins and maxes; first need data_arrays for all cases
    from e3sm_utils import get_data
    data_arrays = [get_data(dataset, field) for dataset in datasets]
    if vmin is None: vmin = min([data.min().values for data in data_arrays])
    if vmax is None: vmax = max([data.max().values for data in data_arrays])
    
    # Open figure
    figure, axes = pyplot.subplots(nrows, ncols, subplot_kw=dict(projection=projection))
    
    # Loop and plot
    for icase, dataset in enumerate(datasets):
        
        # Get data
        data = get_data(dataset, field)
        
        # Get coordinate variables
        lon = get_data(dataset, 'lon')
        lat = get_data(dataset, 'lat')
        
        # Make plot of this data array
        ax = figure.add_axes(axes.ravel()[icase])
        
        # fix plot area; not sure why all this is needed, bugs in cartopy?
        if lat_bounds is not None: 
            ax.set_extent([-180, 180, lat_bounds[0], lat_bounds[1]], crs.PlateCarree())
        if projection == crs.NorthPolarStereo(): cartopy_circular(ax)
        
        # Make plot
        pl = plot_map(lon, lat, data, 
                      transform=crs.PlateCarree(), vmin=vmin, vmax=vmax)
        
        # Label plot
        if labels is not None:
            ax.set_title(labels[icase])
        elif 'case' in data.attrs:
            ax.set_title(data.case)
            
    # Common colorbar
    cb = pyplot.colorbar(pl, ax=axes.ravel().tolist(),
                         orientation='horizontal', shrink=0.8, pad=0.02, 
                         label='%s (%s)'%(data.long_name, data.units))

    return figure
            
def compare_timeseries_2d(data_arrays, cases, **kwargs):
    from matplotlib import pyplot

    # Make sure data_arrays match up with cases
    assert(len(data_arrays) == len(cases))
    
    # Make sure we only have two cases so we can make differences
    assert(len(cases) == 2)
    
    # Open figure
    figure, axes = pyplot.subplots(len(cases) + 1, 1, sharex=True, sharey=True)
    
    # Plot each case, plus differences
    for icase, (case, data) in enumerate(zip(cases, data_arrays)):

        ax = figure.add_axes(axes[icase])
        pl = ax.pcolor(data.time, data.z, data.transpose('z', 'time'))
        cb = pyplot.colorbar(pl, orientation='vertical')

        if icase == 0:
            data_cntl = data.copy(deep=True)
        else:
            data_diff = data - data_cntl
            data_diff.attrs = data.attrs
            vmax = max(abs(data_diff.max().values), abs(data_diff.min().values))
            vmin = -vmax
            ax = figure.add_axes(axes[icase + 1])
            pl = ax.pcolor(data_diff.time, data_diff.z, data_diff.transpose('z', 'time'), cmap='RdBu_r', vmin=vmin, vmax=vmax)
            cb = pyplot.colorbar(pl, orientation='vertical')
            
    return figure


def compare_zonal_means(data_arrays, labels=None, **kwargs):
    
    figure = pyplot.figure()
    ax = figure.add_subplot(111)
    
    for icase, data in enumerate(data_arrays):
        
        # Make plot
        if labels is not None: label = labels[icase]
        pl = ax.plot(data.lat, data, label=label, **kwargs)
        
    # Label using last used data
    ax.set_xlabel('Latitude')
    ax.set_ylabel('%s (%s)'%(data.long_name, data.units))
    ax.legend()
    
    return figure

    
def compare_zonal_profiles(data_arrays, labels=None, nrows=None, ncols=None, **kwargs):
        
    if nrows is None: nrows = len(data_arrays)
    if ncols is None: ncols = 1 #len(data_arrays)
    figure, axes = pyplot.subplots(nrows, ncols, sharex=True, sharey=True)
    
    # Find min and max over all data arrays
    vmin = min([data.min().values for data in data_arrays])
    vmax = max([data.max().values for data in data_arrays])
    
    for icase, data in enumerate(data_arrays):
                
        # Make plot
        ax = figure.add_axes(axes.ravel()[icase])
        pl = ax.pcolormesh(
            data.lat, data.lev, data.transpose(data.lev.name, data.lat.name),
            vmin=vmin, vmax=vmax, cmap='viridis'
        )
        
        # Label this plot
        if labels is not None:
            ax.set_title(labels[icase])
        elif 'case' in data.attrs:
            ax.set_title(data.case)
    
        ax.set_xlabel('Latitude')
        ax.set_ylabel('Hybrid level')
        
    # Make a common label axes
    #ax_common = suplabels(xlabel='Latitude', ylabel=data.lev.name)
    
    pyplot.subplots_adjust(hspace=0.5)
    
    # common colorbar
    cb = pyplot.colorbar(
        pl, ax=axes.ravel().tolist(),
        orientation='horizontal', 
        label='%s (%s)'%(data.long_name, data.units),
        shrink=0.8
    )
    ax.set_ylim(ax.get_ylim()[::-1])
    
    return figure


def suplabels(xlabel=None, ylabel=None, xpad=None, ypad=None):
    # add a big axes, hide frame
    figure = pyplot.gcf()
    ax_common = figure.add_subplot(111, frameon=False)
    ax_common.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    if xlabel is not None: ax_common.set_xlabel(xlabel, labelpad=xpad)
    if ylabel is not None: ax_common.set_ylabel(ylabel, labelpad=ypad)
    return ax_common

