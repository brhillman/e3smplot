#!/usr/bin/env python3
from matplotlib import pyplot
from cartopy import crs

from .e3sm_utils import get_data

def plot_field(dataset, field, plot_type='map', **kwargs):

    # Get data
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
    
    



def compare_maps_diff(lon, lat, data1, data2, labels=('case 1', 'case 2'), 
                      ncols=3, nrows=1, vmin=None, vmax=None, 
                      cmap='viridis', cmap_diff='RdBu_r', **kwargs):
    
    from matplotlib import pyplot
    from cartopy import crs
    
    # setup figure
    figure, axes = pyplot.subplots(
        nrows, ncols, figsize=(15, 3), 
        subplot_kw=dict(projection=crs.PlateCarree())
    )
    
    # Find data limits
    if vmin is None: vmin = min(data1.min().values, data2.min().values)
    if vmax is None: vmax = max(data1.max().values, data2.max().values)
    
    # Loop over cases and plot
    for icase, (data, case) in enumerate(zip((data1, data2), labels)):
        ax = figure.add_axes(axes[icase])
        pl = plot_map(
            lon, lat, data, vmin=vmin, vmax=vmax, 
            transform=crs.PlateCarree(), cmap=cmap,
            **kwargs
        )
        
        # Label plot
        ax.set_title('%s (min = %.3f, max = %.3f)'%(
                     case, data.min().values, data.max().values))

        # Add a colorbar
        cb = pyplot.colorbar(pl, ax=ax, orientation='horizontal', 
                             label='%s (%s)'%(data.long_name, data.units))
    
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
            
            # Make sure cmap is not in kwargs
            if 'cmap' in kwargs.keys(): kwargs.pop('cmap')
            
            # plot
            ax = figure.add_axes(axes[-1])
            pl = plot_map(
                lon, lat, data_diff, 
                cmap=cmap_diff, vmin=vmin, vmax=vmax,
                transform=crs.PlateCarree(),
                **kwargs
            )
            ax.set_title('Difference (min = %.3f, max = %.3f)'%(
                data_diff.min().values, data_diff.max().values
            ))
            
            # Add a colorbar
            cb = pyplot.colorbar(pl, ax=ax, orientation='horizontal', 
                                 label='%s (%s)'%(data.long_name, data.units))
    
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
    from xarray import DataArray
    new_lon = DataArray(numpy.where(lon > 180, lon - 360, lon), dims=('ncol'))
    
    # Setup plot axes
    ax = pyplot.gca()
    ax.set_global()
    ax.coastlines()

    # Make plot
    if all([v is not None for v in (lon_corners, lat_corners)]):
        pl = plot_map_native(lon_corners, lat_corners, data, **kwargs)
    elif 'ncol' in data.dims:
        # Drop missing data
        new_lon = new_lon.where(data.squeeze().notnull()).dropna('ncol')
        lat = lat.where(data.squeeze().notnull()).dropna('ncol')
        data = data.squeeze().dropna('ncol')

        # Plot
        pl = ax.tripcolor(new_lon.squeeze(), lat.squeeze(), data.squeeze(), **kwargs)
    else:
        pl = ax.pcolormesh(new_lon.squeeze(), lat.squeeze(), numpy.ma.masked_invalid(data.squeeze().transpose('lat', 'lon')), **kwargs)
    
    # Return plot handle
    return pl


def compare_maps(data_arrays, labels=None, 
                 ncols=None, nrows=None,
                 lat_bounds=None,
                 projection=crs.PlateCarree(), 
                 vmin=None, vmax=None, label_minmax=False, **kwargs):
    
    # TODO: allow for difference plots
    
    # Figure out size of figure
    if ncols is None: ncols = len(data_arrays)
    if nrows is None: nrows = 1
        
    # Find common mins and maxes
    if vmin is None: vmin = min([data.min().values for data in data_arrays])
    if vmax is None: vmax = max([data.max().values for data in data_arrays])
    
    # Open figure
    figure, axes = pyplot.subplots(nrows, ncols, subplot_kw=dict(projection=projection), figsize=(15,5))
    
    # Loop and plot
    for icase, data in enumerate(data_arrays):
        
        # Make plot of this data array
        ax = figure.add_axes(axes.ravel()[icase])
        
        # fix plot area; not sure why all this is needed, bugs in cartopy?
        if lat_bounds is not None: 
            ax.set_extent([-180, 180, lat_bounds[0], lat_bounds[1]], crs.PlateCarree())
        if projection == crs.NorthPolarStereo(): cartopy_circular(ax)
        
        pl = plot_map(data.lon, data.lat, data, 
                      transform=crs.PlateCarree(), 
                      vmin=vmin, vmax=vmax, **kwargs)
        
        # Label plot
        if labels is None:
            if 'case' in data.attrs: label = data.case
        else:
            label = labels[icase]

        if label_minmax:
            label = '%s (%.2e, %.2e)'%(label, data.min().values, data.max().values)

        ax.set_title(label)
            
    # Common colorbar
    cb = pyplot.colorbar(pl, ax=axes.ravel().tolist(),
                         orientation='horizontal', shrink=0.8, pad=0.02, 
                         label='%s (%s)'%(data.long_name, data.units))

    return figure


def compare_maps_from_ds(datasets, field, plot_diffs=False, vmin=None, vmax=None, cmap='viridis', **kwargs):
    
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
        return compare_maps_diff(data_arrays[0].lon, data_arrays[0].lat, *data_arrays, 
                                 vmin=vmin, vmax=vmax, cmap=cmap, **kwargs)
    else:
        return compare_maps(data_arrays, vmin=vmin, vmax=vmax, cmap=cmap, **kwargs)
    

    
def compare_maps_from_datasets(datasets, field, labels=None, projection=crs.PlateCarree(), 
                               ncols=None, nrows=None, vmin=None, vmax=None,
                               lat_bounds=None, lon_bounds=None, plot_differences=False,
                               **kwargs):
    
    # Figure out size of figure
    if ncols is None: ncols = len(datasets)
    if nrows is None: nrows = 1
        
    # Find common mins and maxes; first need data_arrays for all cases
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
                      transform=crs.PlateCarree(), 
                      vmin=vmin, vmax=vmax, **kwargs)
        
        # Label plot
        if labels is not None:
            ax.set_title(labels[icase])
        elif 'case' in data.attrs:
            ax.set_title(data.case)
            
    # Common colorbar
    if not plot_diffs:
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


def calculate_zonal_mean(data, weights, old_lat, lat_edges=None):
    
    import numpy as numpy
    
    # Mask data weights
    weights = weights.where(data.notnull())
    
    # Calculate new latitudes
    if lat_edges is None:
        lat_edges = numpy.linspace(-90, 90, 31)
        
    # Calculating zonal mean for each latitude by binning data according to lat values
    lat_centers = numpy.zeros(len(lat_edges) - 1)
    if 'lev' in data.dims:
        data_zonal = numpy.zeros([len(lat_edges) - 1, len(data['lev'])])
    else:
        data_zonal = numpy.zeros(len(lat_edges) - 1)
    for ilat in range(len(lat_edges) - 1):
        
        # Find latitude bounds
        lat1 = lat_edges[ilat]
        lat2 = lat_edges[ilat+1]
        
        # Calculate latitude centers from bounds
        lat_centers[ilat] = (lat_edges[ilat+1] + lat_edges[ilat]) / 2.0

        # Calculate mean for this latitude band
        data_band = data.where(old_lat > lat1).where(old_lat <= lat2)
        weights_band = weights.where(old_lat > lat1).where(old_lat <= lat2)
        data_zonal[ilat, ...] = (weights_band * data_band).sum(dim='ncol') / weights_band.sum(dim='ncol')
        
    # Turn these into DataArrays
    from xarray import DataArray
    lat_centers = DataArray(
        lat_centers,
        dims=('lat',),
        attrs={'long_name': 'Latitude', 'units': 'Degrees north'}
    )
    if 'lev' in data.dims:
        dims = ('lat', 'lev')
        coords = {'lat': lat_centers, 'lev': data.lev}
    else:
        dims = ('lat')
        coords = {'lat': lat_centers}

    data_zonal = DataArray(
        data_zonal,
        dims=dims, coords=coords, attrs=data.attrs
    )
    return data_zonal, lat_centers


def compare_zonal_means(datasets, field, labels=None, **kwargs):

    figure = pyplot.figure()
    ax = figure.add_subplot(111)

    for icase, dataset in enumerate(datasets):

        # Get data
        data = get_data(dataset, field)
        area = get_data(dataset, 'area')
        lat = get_data(dataset, 'lat')
        
        # Get area weights for averaging
        from xarray import broadcast
        weights, *__ = broadcast(area, data)

        # Calculate zonal means
        data_zonal, lat_centers = calculate_zonal_mean(data, weights, lat)
        
        # Make plot
        if labels is not None: 
            label = labels[icase]
        else:
            if 'casename' in dataset.attrs.keys():
                label = dataset.attrs['casename']
            elif 'case' in dataset.attrs.keys():
                label = dataset.attrs['case']
            else:
                raise NameError('No valid label.')

        pl = ax.plot(lat_centers, data_zonal, label=label, **kwargs)
        
    # Label using last used data
    ax.set_xlabel('Latitude')
    ax.set_ylabel('%s (%s)'%(data.long_name, data.units))
    ax.legend()
    
    return figure

    
def compare_zonal_profiles(datasets, field, labels=None, 
                           nrows=None, ncols=None, 
                           common_colorbar=False,
                           plot_diffs=True, **kwargs):
        
    if nrows is None: 
        if plot_diffs is True:
            nrows = len(datasets) + 1
        else:
            nrows = len(datasets)
    if ncols is None: ncols = 1 #len(data_arrays)
    figure, axes = pyplot.subplots(nrows, ncols, sharex=True, sharey=True, figsize=(5*ncols, 3*nrows))

    # First, calculate zonal averages and build list of data_arrays to iterate over
    data_arrays = []
    latitudes = []
    for dataset in datasets:
        data = get_data(dataset, field)
        area = get_data(dataset, 'area')
        lat = get_data(dataset, 'lat')

        from xarray import broadcast
        weights, *__ = broadcast(area, data)
        data_zonal, lat_zonal = calculate_zonal_mean(data, weights, lat)
        data_arrays.append(data_zonal)
        latitudes.append(lat_zonal)

    # Find min and max over all data arrays
    vmin = min([data.min().values for data in data_arrays])
    vmax = max([data.max().values for data in data_arrays])
    
    # Loop and plot
    for icase, (data, lat) in enumerate(zip(data_arrays, latitudes)):
                
        # Make plot
        ax = figure.add_axes(axes.ravel()[icase])
        pl = ax.pcolormesh(
            lat, data.lev, data.transpose(data.lev.name, data.lat.name),
            vmin=vmin, vmax=vmax, **kwargs
        )
        cb = pyplot.colorbar(
            pl, ax=ax, orientation='horizontal',
            label='%s (%s)'%(data.long_name, data.units),
            pad=0.2,
        )
        
        # Label this plot
        if labels is not None:
            ax.set_title(labels[icase])
        elif 'case' in dataset.attrs:
            ax.set_title(dataset.case)
        elif 'casename' in dataset.attrs:
            ax.set_title(dataset.casename)
    
        ax.set_xlabel('Latitude')
        ax.set_ylabel('Hybrid level')

        # Maybe plot diffs
        if plot_diffs:
            if icase == 0:
                data_cntl = data.copy(deep=True)
            else:
                data_diff = data - data_cntl
                data_diff.attrs = data.attrs

                # Find min and max
                vmax = abs(data_diff).max().values
                vmin = -vmax

                # Plot diff
                ax = figure.add_axes(axes.ravel()[-1])
                pl = ax.pcolormesh(
                    lat, data_diff.lev, data_diff.transpose(data_diff.lev.name, data_diff.lat.name),
                    vmin=vmin, vmax=vmax, cmap='RdBu_r'
                )
                cb = pyplot.colorbar(
                    pl, ax=ax, orientation='horizontal',
                    label='%s (%s)'%(data_diff.long_name, data_diff.units),
                    pad=0.2,
                )

                # Label this plot
                dmin = data_diff.min().values
                dmax = data_diff.max().values
                davg = data_diff.mean().values
                ax.set_title('Difference (%0.2e, %0.2e, %0.2e)'%(dmin, dmax, davg))
                
                ax.set_xlabel('Latitude')
                ax.set_ylabel('Hybrid level')

    # Make a common label axes
    #ax_common = suplabels(xlabel='Latitude', ylabel=data.lev.name)

    #figure.subplots_adjust(hspace=0.1)
    
    if common_colorbar:
        cb = pyplot.colorbar(
            pl, ax=axes.ravel().tolist(),
            orientation='horizontal', 
            label='%s (%s)'%(data.long_name, data.units),
            shrink=0.8
        )
    ax.set_ylim(ax.get_ylim()[::-1])
    
    return figure


def compare_zonal_profiles_da(data_arrays, labels=None, nrows=None, ncols=None, **kwargs):
        
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
            vmin=vmin, vmax=vmax, **kwargs
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

