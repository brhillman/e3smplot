#!/usr/bin/env python3
import plot_image_gdal

# full path of input file
history_root = '/Users/bhillma/nersc_scratch/scream/regrid/fewer-p3-checks.FSCREAM-HR.ne1024np4_360x720cru_oRRS15to5.cori-knl_intel.1536x8x16.20200206-1553'
map_file_name = f'{history_root}/BlueMarble_land_shallow_2011_180centered.nc'

# Make plot
plot_image_gdal.main(map_file_name, 'bluemarble_test.png')
