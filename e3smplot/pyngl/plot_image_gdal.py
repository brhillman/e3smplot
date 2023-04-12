#!/usr/bin/env python3
import ngl
import xarray as xr
import numpy as np
import warnings
import os

# Increase available buffer size for pyNGL workstations.
# This is a workaround for errors that look like: 
# fatal:ContourPlotDraw: Workspace reallocation would exceed maximum size
def increase_workspace_memory(value=10000000000):
    ws_id = ngl.get_workspace_id()
    rlist = ngl.Resources()
    rlist.wsMaximumSize = value
    ngl.set_values(ws_id,rlist)

# Plot a georeferenced RGB image
# Inspired by:
#  https://www.dkrz.de/up/services/analysis/visualization/sw/ncl/examples/source_code/dkrz-ncl-jpeg-images-as-background-map-overlayed-by-clouds-example
# NOTE: requires a georeferenced image converted to netcdf format, which can be
# done using gdal:
#
#     gdal_translate -if Int16 -of netcdf land_shallow_topo_2048.jpg land_shallow_topo_2048.nc
#
# where land_shallow_topo_2048.jpg is a blue marble jpg. See NCL documentation
# for more details:
#
#     https://www.ncl.ucar.edu/Applications/topo.shtml
#     https://www.ncl.ucar.edu/Applications/rgbacolor.shtml#ex9
#
def plot_image_gdal(
        wks, band1, band2, band3, 
        **kwargs):

    # Clip if over maximum
    band1[band1>255] = 255
    band2[band2>255] = 255
    band3[band3>255] = 255

    # Construct RGBA colormaps
    ramp = np.linspace(0.,1.,255)
    reds = np.zeros((255,4))
    greens = np.zeros((255,4))
    blues = np.zeros((255,4))
    reds[:,0]   = ramp
    greens[:,1] = ramp
    blues[:,2]  = ramp
    # Set alpha channel (a = 1: 100% opaque; a = 0: 100% transparent)
    reds[:,3]     = 1.0 #0.8     
    greens[:,3]   = 0.0 #0.44444
    blues[:,3]    = 0.0 #0.3077

    # Generate coordinate variables since they do not exist in the data
    y = ngl.fspan(-90,   90, band1.shape[0])
    x = ngl.fspan(-180, 180, band1.shape[1])

    # Set up plot resources
    mapres = ngl.Resources()
    mapres.nglDraw               = True
    mapres.nglFrame              = False
    mapres.cnFillOn              = True
    mapres.cnLinesOn             = False
    mapres.cnLineLabelsOn        = False
    mapres.cnInfoLabelOn         = False
    mapres.cnFillMode            = 'RasterFill'
    mapres.cnLevelSelectionMode  = "ExplicitLevels"
    mapres.cnLevels              = np.arange(0,254,1)
    mapres.cnFillBackgroundColor = [1., 1., 1., 1.]  #set it to black
    mapres.cnLinesOn             = False
    mapres.cnLineLabelsOn        = False
    mapres.lbLabelBarOn          = False
    mapres.sfXArray              = x
    mapres.sfYArray              = y
    # NOTE: pyngl will complain about these for the non-map plots,
    # but I want everything set ahead of time so that I can override
    # with kwargs passed into this function
    mapres.mpFillOn         = False
    mapres.mpGridAndLimbOn  = False
    # Set additional plot resource passed via kwargs
    for key, value in kwargs.items(): setattr(mapres, key, value)

    # The green and blue maps will be overlaid on top of the red map,
    # The red map serves as the base upon which the others are overlaid,
    # so we can make this one an actual map. NOTE: I am manually overlaying
    # these plots, because I want them all to use the same plot resources.
    mapres.cnFillColors = reds
    redMap = ngl.contour_map(wks,band1,mapres)

    mapres.cnFillColors = greens
    greenMap = ngl.contour_map(wks,band2,mapres)

    mapres.cnFillColors = blues
    blueMap  = ngl.contour_map(wks,band3,mapres)

    #map_plot = redMap
    #ngl.overlay(map_plot, redMap)
    #ngl.overlay(map_plot, greenMap)
    #ngl.overlay(map_plot, blueMap)

    # Return plot objects
    return redMap, greenMap, blueMap


def main(inputfile, outputfile, **kwargs):

    # output figure type and name
    fig_name = os.path.splitext(outputfile)[0]
    fig_type = os.path.splitext(outputfile)[1][1:]

    # create the plot workstation
    wks = ngl.open_wks(fig_type, fig_name)

    # Increase maximum buffer memory
    increase_workspace_memory(value=10000000000)

    # Read image data
    ds = xr.open_dataset(inputfile)
    band1 = ds['Band1'].values
    band2 = ds['Band2'].values
    band3 = ds['Band3'].values

    # Plot image
    # mpCenterLonF=150 for Pacific; 340 for Atlantic; 200 Pacific v2
    # mpCenterLatF=30 for Pacific; 45 for Atlantic; 40 Pacific v2
    pl = plot_image_gdal(wks, band1, band2, band3, **kwargs) #, mpCenterLonF=340, mpCenterLatF=45)

    # Finalize plot
    ngl.frame(wks)
    ngl.end()


if __name__ == '__main__':
    import plac; plac.call(main)
