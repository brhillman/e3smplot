#!/usr/bin/env python3

from glob import glob
import sys, os
import xarray
import ngl
import numpy
from e3smplot.e3sm_utils import area_average, open_dataset, get_data, get_mapping_file
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
test_case_name = 'Run 2'

# Set control case names. We will compare the test_case_name against each of
# these
cntl_case_names = ('Run 1', 'CERES-SYN', 'ERA5') #'ne1024 Oct9', 'CERES-SYN', 'ERA5')

# Path were we can find netcdf files with output from the test/model case
data_paths = {
    # Raw files
    #'Run 1' : '/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/20201112.SCREAMv0dyamond2.F2010-SCREAM-HR-DYAMOND2.ne1024pg2_r0125_oRRS18to6v3.cori-knl.1536x8x16/run',
    #'Run 2': '/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/SCREAMv0.SCREAM-DY2.ne1024pg2.20201127/run',
    # Regridded files
    'Run 1' : '/global/cfs/cdirs/e3sm/terai/SCREAM/DYAMOND2/Output/20201112/regridded',
    'Run 2' : '/global/cfs/cdirs/e3sm/terai/SCREAM/DYAMOND2/Output/20201127/regridded',
    'ne1024 Oct9' : '/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/master.ne1024pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.1536x8x16.DY2_Oct9.20201009-16/run',
    'ne256 lamlow': '/global/cscratch1/sd/bogensch/E3SM_simulations/master.ne256pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.1024x16x8.DY2_Oct6.20201022-16.lamlow.001a/run',
    'ne256 rrtmg' : '/global/cscratch1/sd/bhillma/scream/cases/9bfb38267.ne256pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.192x16x8.rrtmg.20201021-0921/run',
    'ne256 rrtmgp': '/global/cscratch1/sd/bhillma/scream/cases/9bfb38267.ne256pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.192x16x8.rrtmgp.20201021-1722/run',
    'ne4 rrtmg': '/global/cscratch1/sd/bhillma/e3sm_scratch/cori-knl/SMS_Ld5.ne4_ne4.FC5AV1C-L.cori-knl_intel.add-optics-outputs/run',
    'ne4 rrtmgp': '/global/cscratch1/sd/bhillma/e3sm_scratch/cori-knl/SMS_Ld5.ne4_ne4.FC5AV1C-L.cori-knl_intel.cam-rrtmgp.add-optics-outputs/run',
    'ne1024 Nov05': '/global/cscratch1/sd/terai/E3SM_simulations/master.ne1024pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.1536x8x16.DY2_Nov05branch_SHOC_P3_AB_bugfix.20201105-16',
    'CERES-SYN': '/global/cfs/cdirs/e3sm/terai/Obs_datasets/CERES',
    'ERA5': '/global/cfs/cdirs/e3sm/bhillma/obs_datasets/ERA5',
    'GPM': '/global/cfs/cdirs/e3sm/terai/Obs_datasets/GPM',
}

# List of fields we want to make map plots for
variables = (
    # 3D vars
    # TODO: do something useful with these
    #'CLOUD', 'EMIS', 'TOT_CLD_VISTAU',
    #'TOT_ICLD_VISTAU', 'LIQ_ICLD_VISTAU', 'ICE_ICLD_VISTAU',
    # 
    # 2D vars
    'FSNTOA', 'FSNTOAC', 'FLNT', 'FLNTC',
    #'PRECT', 
    #'CLDTOT', 'CLDLOW', 'CLDMED', 'CLDHGH',
    'TMCLDLIQ', 'TMCLDICE',
    'SHFLX', 'TREFHT', 'TMQ',
)

# Glob strings for searching for case files. For model cases, you may want to
# search for different history tapes (h0, h1, h2, etc). 
# For some of the obs datasets, we have separate files for different categories
# of fields.
glob_strings = {
    'ERA5': {'PRECT' : 'ERA5_surf_2020*.nc',
             'SHFLX' : 'ERA5_surf_2020*.nc',
             'TREFHT': 'ERA5_surf_2020*.nc',
             'TMQ'   : 'ERA5_surf_2020*.nc',},
    'GPM': {'PRECT' : '*.nc',},
    'CERES-SYN': {v: '*.nc' for v in variables},
    # For model-model comparisons, need to specify history tape number
    'ne1024 Oct9': {v: '*.h1.*.nc' for v in variables},
    'ne1024 rrtmg': {v: '*.h1.*.nc' for v in variables},
    'ne256 rrtmg' : {v: '*.h1.*.nc' for v in variables},
    'ne256 rrtmgp': {v: '*.h1.*.nc' for v in variables},
    'ne256 lamlow': {v: '*.h1.*.nc' for v in variables},
    'ne4 rrtmg': {v: '*.h0.*.nc' for v in variables},
    'ne4 rrtmgp': {v: '*.h0.*.nc' for v in variables},
    # Map variables to files for SCREAM-HR production run.
    # TODO: make this more flexible
    'Run 1': {v: '*.eam.h[0-9]*.nc' for v in variables},
    'Run 2': {v: 'SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h[0-9]*.nc' for v in variables},
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
time_offsets['ne4 rrtmgp'] = datetime.timedelta(days=(365*2000))
time_offsets['ne4 rrtmg'] = datetime.timedelta(days=(365*2000))

# Where should we write the plot files?
tmp_name = test_case_name.replace(' ', '_')
graphics_root = f'./graphics/{tmp_name}'

# Flags to control what kind of plots to make. You probably want to make them
# all, but being able to disable for testing is useful. Some of them take quite
# a while to make.
do_contour_maps = False
do_time_series  = True
do_zonal_means  = False
do_anomalies    = False

# Before we do anything, make sure our plot directory exists
os.makedirs(graphics_root, exist_ok=True)

# Start with anomaly plots
if do_anomalies:

    # Find files
    files = [sorted(glob(f'{data_paths[c]}/{glob_strings[c]}')) for c in (test_case_name, *cntl_case_names)]
    print(files)

    exit()

for vname in variables:
    # Find files for which we can retrieve variable from
    if do_contour_maps:

        # For map plots we only compare two at a time so we can look at differences
        for cntl_case_name in cntl_case_names:

            # Find files and make sure we can retrieve variable
            if vname not in glob_strings[test_case_name] or vname not in glob_strings[cntl_case_name]: continue
            test_files = sorted(glob(f'{data_paths[test_case_name]}/{glob_strings[test_case_name][vname]}'))
            cntl_files = sorted(glob(f'{data_paths[cntl_case_name]}/{glob_strings[cntl_case_name][vname]}'))
            test_files = [f for f in test_files if can_retrieve_field(f, vname)]
            cntl_files = [f for f in cntl_files if can_retrieve_field(f, vname)]
            if len(test_files) == 0 or len(cntl_files) == 0:
                print('Could not find {} in test or cntl files'.format(vname))
                continue

            # Figure out what mapping files we need based on test and obs data
            if get_grid_name(test_files[0]) != get_grid_name(cntl_files[0]):
                map_file = get_mapping_file(test_files[0], cntl_files[0], mapping_root)
            else:
                map_file = None

            # Compare maps
            print(f'Making {test_case_name} vs {cntl_case_name} contour maps for {vname}...')
            figname = f'{graphics_root}/{vname}_{test_case_name.replace(" ", "_")}_vs_{cntl_case_name.replace(" ", "_")}_maps.png'
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
        #files, names = zip(*[(f,n) if can_retrieve_field(f, vname) else (None, None) 
        #                     for (f,n) in zip(files, names)])
        ffiles = []
        nnames = []
        for file_list, name in zip(files, names):
            files_with_var = [f for f in file_list if can_retrieve_field(f, vname)]
            if len(files_with_var) > 0:
                ffiles.append(files_with_var)
                nnames.append(name)
        files, names = ffiles, nnames
        if len(files) == 0: continue

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
        #files, names = zip(*[(f,n) for (f,n) in zip(files, names) if  can_retrieve_field(f, vname)])
        ffiles = []
        nnames = []
        for file_list, name in zip(files, names):
            files_with_var = [f for f in file_list if can_retrieve_field(f, vname)]
            if len(files_with_var) > 0:
                ffiles.append(files_with_var)
                nnames.append(name)
        files, names = ffiles, nnames
        if len(files) == 0: continue

        # Make time series plots
        print(f'Making time series for {vname}...')
        figname = f'{graphics_root}/time_series/{vname}_{test_case_name.replace(" ", "_")}_vs_all_timeseries.png'
        os.makedirs(os.path.dirname(figname), exist_ok=True)
        #if min([get_data(open_dataset(f), vname).min().values for f in files]) > 0: kwargs['ymin'] = 0
        compare_time_series.main(files, vname, figname, names, time_offsets=[time_offsets[n] for n in names])

