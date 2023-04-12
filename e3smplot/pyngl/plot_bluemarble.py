#!/usr/bin/env python3
import plot_image_gdal
import os

# full path of input file
history_root = '/Users/bhillma/nersc_scratch/scream/regrid/fewer-p3-checks.FSCREAM-HR.ne1024np4_360x720cru_oRRS15to5.cori-knl_intel.1536x8x16.20200206-1553'
#map_file_name = f'{history_root}/BlueMarble_land_shallow_2011_180centered.nc'
#map_file_name = f'land_shallow_topo_2048.nc'
map_file_name = f'world.200401.3x5400x2700.nc'

# Make plot
plot_name = '/global/cfs/cdirs/e3sm/www/bhillma/test_graphics/bluemarble.png'
os.makedirs(f'{os.path.dirname(plot_name)}', exist_ok=True)
plot_image_gdal.main(
    map_file_name, plot_name, 
    mpProjection='Orthographic', mpCenterLatF=60, mpCenterLonF=-100
)
os.system(f'chmod a+rx {os.path.dirname(plot_name)}')
os.system(f'chmod a+r {plot_name}')
