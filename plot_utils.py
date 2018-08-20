#!/usr/bin/env python3
from matplotlib import pyplot
from cartopy import crs

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


def plot_map(lon, lat, data, **kwargs):
    # Fix longitude
    import numpy
    new_lon = numpy.where(lon > 180, lon - 360, lon)
    
    # Make plot
    ax = pyplot.gca()
    ax.set_global()
    ax.coastlines()
    if 'ncol' in data.dims:
        pl = ax.tripcolor(new_lon, lat, data, **kwargs)
    else:
        pl = ax.pcolormesh(new_lon, lat, data.transpose('lat', 'lon'), **kwargs)
        
    # Add colorbar
    cb = pyplot.colorbar(
        pl, orientation='horizontal', shrink=0.8, pad=0.02, 
        label='%s (%s)'%(data.long_name, data.units)
    )
    
    return pl, cb


def compare_maps_diff(lon, lat, data1, data2, labels, **kwargs):
    
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


def plot_map(lon, lat, data, **kwargs):
    # Fix longitude
    import numpy
    new_lon = numpy.where(lon > 180, lon - 360, lon)
    
    # Make plot
    ax = pyplot.gca()
    ax.set_global()
    ax.coastlines()
    if 'ncol' in data.dims:
        pl = ax.tripcolor(new_lon, lat, data, **kwargs)
    else:
        pl = ax.pcolormesh(new_lon, lat, data.transpose('lat', 'lon'), **kwargs)
    
    return pl


def compare_maps(data_arrays, labels=None, 
                 ncols=None, nrows=None, 
                 projection=crs.PlateCarree(), **kwargs):
    
    # TODO: allow for difference plots
    
    # Figure out size of figure
    if ncols is None: ncols = len(data_arrays)
    if nrows is None: nrows = 1
        
    # Find common mins and maxes
    vmin = min([data.min().values for data in data_arrays])
    vmax = max([data.max().values for data in data_arrays])
    
    # Open figure
    figure, axes = pyplot.subplots(nrows, ncols, subplot_kw=dict(projection=projection))
    
    # Loop and plot
    for icase, data in enumerate(data_arrays):
        
        # Make plot of this data array
        ax = figure.add_axes(axes.ravel()[icase])
        pl = plot_map(data.lon, data.lat, data.transpose('lat', 'lon'), 
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