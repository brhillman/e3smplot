#!/usr/bin/env python3

from glob import glob
import sys, os
import xarray
import ngl
import numpy
from e3smplot.e3sm_utils import area_average, open_dataset, get_mapping_file
from e3smplot.e3sm_utils import get_grid_name, can_retrieve_field, is_latlon
from e3smplot.pyngl import compare_maps
from e3smplot.mpl import compare_time_series, compare_zonal_means
import os.path
import datetime

#
# Stuff the user will probably want to change
#

# Where we keep our grids and map files. Make sure this is in scratch space,
# these things can get quite large
mapping_root = '/global/cscratch1/sd/bhillma/grids'

# A short, meaningful name with which to label plots and set output filenames.
# Does not need to match the case name of the model run, but does need to match
# one of the data_paths keys below!
test_case_name = 'ne256 rrtmgp'

# Set control case names. We will compare the test_case_name against each of
# these
cntl_case_names = ('ne256 rrtmg', 'CERES-SYN', 'ERA5')

# Path were we can find netcdf files with output from the test/model case
data_paths = {
    'ne1024'      : '/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/master.ne1024pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.1536x8x16.DY2_Oct9.20201009-16/run',
    'ne256 lamlow': '/global/cscratch1/sd/bogensch/E3SM_simulations/master.ne256pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.1024x16x8.DY2_Oct6.20201022-16.lamlow.001a/run',
    'ne256 rrtmg' : '/global/cscratch1/sd/bhillma/scream/cases/9bfb38267.ne256pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.192x16x8.rrtmg.20201021-0921/run',
    'ne256 rrtmgp': '/global/cscratch1/sd/bhillma/scream/cases/9bfb38267.ne256pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.192x16x8.rrtmgp.20201021-1722/run',
    'CERES-SYN': '/global/cfs/cdirs/e3sm/terai/Obs_datasets/CERES',
    'ERA5': '/global/cfs/cdirs/e3sm/bhillma/obs_datasets/ERA5',
    'GPM': None, # GPM not supported yet
}

# List of fields we want to make map plots for
variables = (
    'FSNTOA',
    'FSNTOAC',
    'FLNT',
    'FLNTC',
    'PRECT',
    'CLDTOT',
    'SHFLX',
    'TREFHT',
    'TMQ',
)

# Glob strings for searching for case files. For model cases, you may want to
# search for different history tapes (h0, h1, h2, etc). 
# For some of the obs datasets, we have separate files for different categories
# of fields.
glob_strings = {
    'ERA5': {'PRECT' : 'ERA5_surf*.nc',
             'SHFLX' : 'ERA5_surf*.nc',
             'TREFHT': 'ERA5_surf*.nc',
             'TMQ'   : 'ERA5_surf*.nc',},
    'CERES-SYN': {v: '*.nc' for v in variables},
    # For model-model comparisons, need to specify history tape number
    'ne1024 rrtmg': {v: '*.h1.*.nc' for v in variables},
    'ne256 rrtmg' : {v: '*.h1.*.nc' for v in variables},
    'ne256 rrtmgp': {v: '*.h1.*.nc' for v in variables},
    'ne256 lamlow': {v: '*.h1.*.nc' for v in variables},
}

# Time offsets; hack so we can compare F-cases.
# HOW TO USE THIS FEATURE: in the event that you want to compare a run for a
# different time period that what we have observations for, you can use this to
# manually adjust the time coordinate of the simulation. This is certainly a
# hack, but allows us to use my code that attempts to find an overlapping time
# period to compare with the correct observations. This is a list of time
# offsets for the "test" (run you want to evaluate) and "cntl" (run or obs you
# want to evaluate against) cases. It should be a datetime.timedelta object. As
# an example for usage, say you meant to start a run for the DYAMOND2 period
# (2020-01-20), but accidentally ran an F-case (0001-01-01). We want a time
# offset in that case of 365 * 2019 + 20 days. Any entries in this list that are
# None type will not add an offset (the default behavior).
#time_offsets = [datetime.timedelta(days=(365*2019+20)), None]
time_offsets = {c: None for c in [test_case_name, *cntl_case_names]}
time_offsets['ne256 lamlow'] = datetime.timedelta(days=(365*2019+20))

# Where should we write the plot files?
graphics_root = './graphics'

# Flags to control what kind of plots to make. You probably want to make them
# all, but being able to disable for testing is useful. Some of them take quite
# a while to make.
do_contour_maps = True
do_time_series = True
do_zonal_means = True

# Before we do anything, make sure our plot directory exists
os.makedirs(graphics_root, exist_ok=True)

for vname in variables:
    if do_contour_maps:

        # For map plots we only compare two at a time so we can look at differences
        for cntl_case_name in cntl_case_names:

            # Find files and make sure we can retrieve variable
            if vname not in glob_strings[test_case_name] or vname not in glob_strings[cntl_case_name]: continue
            test_files = sorted(glob(f'{data_paths[test_case_name]}/{glob_strings[test_case_name][vname]}'))
            cntl_files = sorted(glob(f'{data_paths[cntl_case_name]}/{glob_strings[cntl_case_name][vname]}'))
            if not can_retrieve_field(cntl_files[0], vname): continue


            # Figure out what mapping files we need based on test and obs data
            if get_grid_name(test_files[0]) != get_grid_name(cntl_files[0]):
                map_file = get_mapping_file(test_files[0], cntl_files[0], mapping_root)
            else:
                map_file = None

            # Compare maps
            print(f'Making {test_case_name} vs {cntl_case_name} contour maps for {vname}...')
            figname = f'{graphics_root}/{vname}_{test_case_name.replace(" ", "_")}_vs_{cntl_case_name.replace(" ", "_")}.png'
            compare_maps.main(test_files, cntl_files, vname, figname,
                    test_name=test_case_name, cntl_name=cntl_case_name,
                    map_file=map_file,
                    time_offsets=(time_offsets[test_case_name], time_offsets[cntl_case_name]),
                    lbAutoManage=False,
                    lbTitleFontHeightF=0.02,
                    lbLabelFontHeightF=0.02,
                    )

    # Make zonal mean plots
    if do_zonal_means:

        # Find files and make sure we can retrieve variable
        files, names = zip(*[
            (sorted(glob(f'{data_paths[n]}/{glob_strings[n][vname]}')), n) 
            for n in [test_case_name, *cntl_case_names] 
            if vname in glob_strings[n]
        ])
        # Only files we can find variable name in
        files, names = zip(*[(f,n) for (f,n) in zip(files, names) if  can_retrieve_field(f[0], vname)])

        # Figure out what mapping files we need based on test and obs data
        default_grid_file = sorted(glob(f'{data_paths["CERES-SYN"]}/*.nc'))[0]
        maps = [get_mapping_file(f[0], default_grid_file, mapping_root)
                if not is_latlon(f[0]) else None for f in files]

        # Make zonal mean plots
        print(f'Making zonal mean plots for {vname}...')
        figname = f'{graphics_root}/{vname}_{test_case_name.replace(" ", "_")}_vs_all_zonal.png'
        compare_zonal_means.main(files, names, vname, figname,
                maps=maps,
                time_offsets=[time_offsets[n] for n in names],
                )


    # Make time series plots
    if do_time_series:

        # Find files and make sure we can retrieve variable
        files, names = zip(*[
            (sorted(glob(f'{data_paths[n]}/{glob_strings[n][vname]}')), n) 
            for n in [test_case_name, *cntl_case_names] 
            if vname in glob_strings[n]
        ])
        # Only files we can find variable name in
        files, names = zip(*[(f,n) for (f,n) in zip(files, names) if  can_retrieve_field(f[0], vname)])

        # Make time series plots
        print(f'Making time series for {vname}...')
        figname = f'{graphics_root}/{vname}_{test_case_name.replace(" ", "_")}_vs_all_timeseries.png'
        compare_time_series.main(files, vname, figname, names, time_offsets=[time_offsets[n] for n in names])

