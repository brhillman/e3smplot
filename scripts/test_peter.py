#!/usr/bin/env python3

from glob import glob
import sys, os
import xarray
import ngl
import numpy
from e3smplot.e3sm_utils import area_average
from e3smplot.pyngl import compare_maps
from e3smplot.mpl import compare_time_series, compare_zonal_means

# List of fields we want to make map plots for
map_vars = (
    #'FSNTOA',
    #'FSNTOAC',
    #'FLNT',
    #'FLNTC',
    'PRECT',
    #'CLDTOT',
    #'SHFLX',
    #'TREFHT',
    #'TMQ',
)

time_series_vars = (
    #'FSNTOA',
    #'TREFHT',
)

zonal_mean_vars = (
    #'FSNTOA',
    #'TREFHT',
)

obs_sources = {
    'FSNTOA' : 'CERES-SYN',
    'FSNTOAC': 'CERES-SYN',
    'FLNT'   : 'CERES-SYN',
    'FLNTC'  : 'CERES-SYN',
    'PRECT'  : 'ERA5',
    'CLDTOT' : 'CERES-SYN',
    'TMQ'    : 'ERA5',
    'SHFLX'  : 'ERA5',
    'TREFHT' : 'ERA5',
}

# Offline mapping weights from model to obs sources
map_files = {
    'CERES-SYN': '/global/cscratch1/sd/bhillma/grids/map_ne256pg2_to_ceres1deg_nco.nc',
    'GPM': None,
    'ERA5': '/global/cscratch1/sd/bhillma/grids/map_ne256pg2_to_era5_nco.nc',
}

test_case_name = 'SCREAM-HR ac25'
#test_data_path = '/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/master.ne1024pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.1536x8x16.DY2_Oct9.20201009-16/run'
test_data_path = '/global/cscratch1/sd/petercal/E3SM_runs/autoconv_25perc.ne256pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.1024x16x8.DY2_Oct6.20201016-12/run'
obs_data_path = {
    'CERES-SYN': '/global/cfs/cdirs/e3sm/terai/Obs_datasets/CERES',
    'GPM': None,
    'ERA5': '/global/cfs/cdirs/e3sm/bhillma/obs_datasets/ERA5',
}
obs_data_prefix = {
    'FSNTOA' : '',
    'FSNTOAC': '',
    'FLNT'   : '',
    'FLNTC'  : '',
    'PRECT'  : 'ERA5_surf',
    'CLDTOT' : '',
    'SHFLX': 'ERA5_surf',
    'TREFHT': 'ERA5_surf',
    'TMQ': 'ERA5_surf',
}
graphics_root = './graphics'
os.makedirs(graphics_root, exist_ok=True)

for vname in map_vars:

    print(f'Processing {vname}...')

    # Find files
    test_files = sorted(glob(f'{test_data_path}/*.h1.*.nc'))
    obs_files = sorted(glob(f'{obs_data_path[obs_sources[vname]]}/{obs_data_prefix[vname]}*.nc'))

    # Compare maps
    figname = f'{graphics_root}/{vname}_{test_case_name}_vs_{obs_sources[vname]}.png'
    compare_maps.main(test_files, obs_files, vname, figname,
            test_name=test_case_name, cntl_name=obs_sources[vname],
            map_file=map_files[obs_sources[vname]],
            lbAutoManage=False,
            lbTitleFontHeightF=0.02,
            lbLabelFontHeightF=0.02,
            )

# Make time series plots
for vname in time_series_vars:

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
for vname in zonal_mean_vars:

    # Find files
    print(f'Making zonal means for {vname}...')
    test_files = sorted(glob(f'{test_data_path}/*.h1.*.nc'))
    obs_files = sorted(glob(f'{obs_data_path[obs_sources[vname]]}/{obs_data_prefix[vname]}*.nc'))

    # Make time series plots
    figname = f'{graphics_root}/{vname}_{test_case_name}_vs_{obs_sources[vname]}_zonal.png'
    compare_zonal_means.main(test_files, obs_files, vname, figname,
            test_name=test_case_name, cntl_name=obs_sources[vname],
            map_file=map_files[obs_sources[vname]],
            )
