#!/usr/bin/env python3

from glob import glob
import sys, os
import xarray
import ngl
import numpy
from e3smplot.e3sm_utils import area_average, open_dataset, get_grid, get_mapping_file
from e3smplot.pyngl import compare_maps
from e3smplot.mpl import compare_time_series, compare_zonal_means
import os.path

#
# Stuff the user will probably want to change
#

# Where we keep our grids and map files. Make sure this is in scratch space,
# these things can get quite large
mapping_root = '/global/cscratch1/sd/bhillma/grids'

# A short, meaningful name with which to label plots and set output filenames.
# Does not need to match the case name of the model run
test_case_name = 'SCREAM-HR'

# Path were we can find netcdf files with output from the test/model case
test_data_path = '/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/master.ne1024pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.1536x8x16.DY2_Oct9.20201009-16/run'

# Where should we write the plot files?
graphics_root = './graphics'

# Flags to control what kind of plots to make. You probably want to make them
# all, but being able to disable for testing is useful. Some of them take quite
# a while to make.
do_contour_maps = True
do_time_series = False
do_zonal_means = True

# List of fields we want to make map plots for
variables = (
    'FSNTOA',
    #'FSNTOAC',
    #'FLNT',
    #'FLNTC',
    #'PRECT',
    #'CLDTOT',
    #'SHFLX',
    #'TREFHT',
    #'TMQ',
)

#
# Things the user may want to change, but probably not
#

# Dict specifying what obs we prefer to use for comparison
obs_sources = {
    'FSNTOA' : 'CERES-SYN',
    'FSNTOAC': 'CERES-SYN',
    'FLNT'   : 'CERES-SYN',
    'FLNTC'  : 'CERES-SYN',
    'PRECT'  : 'GPM',
    'CLDTOT' : 'CERES-SYN',
    'TMQ'    : 'ERA5',
    'SHFLX'  : 'ERA5',
    'TREFHT' : 'ERA5',
}

# Paths to obs datasets where we can find netcdf files. You may need to modify
# this for your system, but all users on the project should have access to these
# locations
obs_data_path = {
    'CERES-SYN': '/global/cfs/cdirs/e3sm/terai/Obs_datasets/CERES',
    'GPM': None, # GPM not supported yet
    'ERA5': '/global/cfs/cdirs/e3sm/bhillma/obs_datasets/ERA5',
}

# For some of the obs datasets, we have separate files for different categories
# of fields. Mostly this is just ERA5, but I could see needing to do something
# like this for things *other* than h1 files from EAM as well. This is here to
# help automate our searches, but he surefire way would still be to hard-code
# the glob lines below
obs_data_prefix = {
    'FSNTOA' : '',
    'FSNTOAC': '',
    'FLNT'   : '',
    'FLNTC'  : '',
    'PRECT'  : '',
    'CLDTOT' : '',
    'SHFLX': 'ERA5_surf',
    'TREFHT': 'ERA5_surf',
    'TMQ': 'ERA5_surf',
}

# Before we do anything, make sure our plot directory exists
os.makedirs(graphics_root, exist_ok=True)

if do_contour_maps:
    for vname in variables:

        print(f'Making contour maps for {vname}...')

        # Find files
        test_files = sorted(glob(f'{test_data_path}/*.h1.*.nc'))
        obs_files = sorted(glob(f'{obs_data_path[obs_sources[vname]]}/{obs_data_prefix[vname]}*.nc'))

        # Figure out what mapping files we need based on test and obs data
        map_file = get_mapping_file(test_files[0], obs_files[0], mapping_root)

        # Compare maps
        figname = f'{graphics_root}/{vname}_{test_case_name}_vs_{obs_sources[vname]}.png'
        compare_maps.main(test_files, obs_files, vname, figname,
                test_name=test_case_name, cntl_name=obs_sources[vname],
                map_file=map_file,
                lbAutoManage=False,
                lbTitleFontHeightF=0.02,
                lbLabelFontHeightF=0.02,
                )

# Make time series plots
if do_time_series:
    for vname in variables:

        # Find files
        print(f'Making time series for {vname}...')
        test_files = sorted(glob(f'{test_data_path}/*.h1.*.nc'))
        obs_files = sorted(glob(f'{obs_data_path[obs_sources[vname]]}/{obs_data_prefix[vname]}*.nc'))

        # Make time series plots
        figname = f'{graphics_root}/{vname}_{test_case_name}_vs_{obs_sources[vname]}_timeseries.png'
        compare_time_series.main(test_files, obs_files, vname, figname,
                test_name=test_case_name, cntl_name=obs_sources[vname],
                )

# Make zonal mean plots
if do_zonal_means:
    for vname in variables:

        # Find files
        print(f'Making zonal means for {vname}...')
        test_files = sorted(glob(f'{test_data_path}/*.h1.*.nc'))
        obs_files = sorted(glob(f'{obs_data_path[obs_sources[vname]]}/{obs_data_prefix[vname]}*.nc'))

        # Figure out what mapping files we need based on test and obs data
        map_file = get_mapping_file(test_files[0], obs_files[0], mapping_root)

        # Make time series plots
        figname = f'{graphics_root}/{vname}_{test_case_name}_vs_{obs_sources[vname]}_zonal.png'
        compare_zonal_means.main(test_files, obs_files, vname, figname,
                test_name=test_case_name, cntl_name=obs_sources[vname],
                map_file=map_file,
                )
