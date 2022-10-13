#!/usr/bin/env python3

from glob import glob
import sys, os
import xarray
import numpy
from e3smplot.e3sm_utils import area_average, open_dataset, get_data, get_mapping_file
from e3smplot.e3sm_utils import get_grid_name, can_retrieve_field, is_latlon
from e3smplot.mpl import compare_maps  # TODO: replace with mpl impl
from e3smplot.mpl import compare_time_series, compare_zonal_means
import os.path
import datetime

#
# Stuff the user will probably want to change
#

# Where we keep our grids and map files. Make sure this is in scratch space,
# these things can get quite large
mapping_root = '/gpfs/alpine/cli115/proj-shared/terai/maps'

# A short, meaningful name with which to label plots and set output filenames.
# Does not need to match the case name of the model run, but does need to match
# one of the data_paths keys below!
test_case_name = 'SCREAMv1.1_firstday'

# Set control case names. We will compare the test_case_name against each of these
cntl_case_names = ('CERES-SYN','ERA5',) #'ne1024 Oct9', 'CERES-SYN', 'ERA5')

# Period to compare
t1, t2 = ('2013-10-01 00:00:00', '2013-10-02 00:00:01')

# Specific time ranges for comparisons; this is OBS-dependent
time_ranges = {
    'CERES-EBAF': ('2013-10-01 00:00:00', '2013-10-02 00:00:01'),
    'CERES-SYN' : ('2013-10-01 00:00:00', '2013-10-02 00:00:01'),
    'ERA5'      : ('2013-10-01 00:00:00', '2013-10-02 00:00:01'),
    'GPM'       : ('2013-10-01 00:00:00', '2013-10-02 00:00:01'),
}

# Path were we can find netcdf files with output from the test/model case
data_paths = {
    # Raw files
    #'Run 1' : '/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/20201112.SCREAMv0dyamond2.F2010-SCREAM-HR-DYAMOND2.ne1024pg2_r0125_oRRS18to6v3.cori-knl.1536x8x16/run',
    #'Run 2': '/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/SCREAMv0.SCREAM-DY2.ne1024pg2.20201127/run',
    # Regridded files
    #'Run 1' : '/global/cfs/cdirs/e3sm/terai/SCREAM/DYAMOND2/Output/20201112/regridded',
    #'Run 2' : '/global/cfs/cdirs/e3sm/terai/SCREAM/DYAMOND2/Output/20201127/regridded',
    #'ne1024 Oct9' : '/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/master.ne1024pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.1536x8x16.DY2_Oct9.20201009-16/run',
    #'ne256 lamlow': '/global/cscratch1/sd/bogensch/E3SM_simulations/master.ne256pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.1024x16x8.DY2_Oct6.20201022-16.lamlow.001a/run',
    #'ne256 rrtmg' : '/global/cscratch1/sd/bhillma/scream/cases/9bfb38267.ne256pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.192x16x8.rrtmg.20201021-0921/run',
    #'ne256 rrtmgp': '/global/cscratch1/sd/bhillma/scream/cases/9bfb38267.ne256pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.192x16x8.rrtmgp.20201021-1722/run',
    #'ne4 rrtmg': '/global/cscratch1/sd/bhillma/e3sm_scratch/cori-knl/SMS_Ld5.ne4_ne4.FC5AV1C-L.cori-knl_intel.add-optics-outputs/run',
    #'ne4 rrtmgp': '/global/cscratch1/sd/bhillma/e3sm_scratch/cori-knl/SMS_Ld5.ne4_ne4.FC5AV1C-L.cori-knl_intel.cam-rrtmgp.add-optics-outputs/run',
    #'ne1024 Nov05': '/global/cscratch1/sd/terai/E3SM_simulations/master.ne1024pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.1536x8x16.DY2_Nov05branch_SHOC_P3_AB_bugfix.20201105-16',
    'SCREAMv1.1_firstday'       : '/gpfs/alpine/cli115/proj-shared/terai/SCREAMv1_output/ne1024pg2_ne1024pg2_7daySim_try2/changed_date',
    'CERES-SYN'      : '/gpfs/alpine/cli115/proj-shared/terai/Obs_datasets/CERES',
    #'CERES-EBAF'    : '/lus/theta-fs0/projects/ClimateEnergy_4/brhillman/ceres-ebaf',
    'ERA5'          : '/gpfs/alpine/cli115/proj-shared/terai/Obs_datasets/ERA5',
    #'GPM'           : '/lus/theta-fs0/projects/ClimateEnergy_4/brhillman/gpm/DYAMOND2_period',
}

# List of fields we want to make map plots for
variables = (
    # 3D vars
    # TODO: do something useful with these
    #'CLOUD', 'EMIS', 'TOT_CLD_VISTAU',
    #'TOT_ICLD_VISTAU', 'LIQ_ICLD_VISTAU', 'ICE_ICLD_VISTAU',
    # 
    # 2D vars
    'SW_flux_up@tom','SW_flux_dn@tom','LW_flux_up@tom','SW_clrsky_flux_up@tom','LW_clrsky_flux_up@tom',
    'VapWaterPath','T_2m','PRECT','ps','surf_evap','surf_sens_flux',
    #'surf_sens_flux',#'PRECT',
)

# Glob strings for searching for case files. For model cases, you may want to
# search for different history tapes (h0, h1, h2, etc). 
# For some of the obs datasets, we have separate files for different categories
# of fields.
glob_strings = {
    'ERA5': {#'PRECT'       : 'ERA5_surf_2020*.nc',
             'SHFLX'       : 'ERA5_surf_2020*.nc',
             'TREFHT'      : 'ERA5_surf_2020*.nc',
             'TMQ'         : 'ERA5_surf_2020*.nc',
             'PS'          : 'ERA5_surf_2020*.nc',
             'WINDSPD_10M' : 'ERA5_surf_2020*.nc',
             'LHFLX'       : 'ERA5_surf_2020*.nc',
             'T_2m'             : 'ERA5_humidity_*.nc',
             'ps'               : 'ERA5_humidity_*.nc',
             'surf_evap'        : 'ERA5_humidity_*.nc',
             'surf_sens_flux'   : 'ERA5_humidity_*.nc',
             'VapWaterPath'     : 'ERA5_humidity_*.nc',
             'PRECT'            : 'ERA5_humidity_*.nc',},
    'GPM'       : {'PRECT' : '*.nc'},
    'CERES-SYN' : {v: f'*.nc' for v in variables},
    'CERES-EBAF': {v: 'CERES_EBAF_Ed4.1_Subset_*.nc' for v in variables},
    # For model-model comparisons, need to specify history tape number
    'ne1024 Oct9' : {v: '*.h1.*.nc' for v in variables},
    'ne1024 rrtmg': {v: '*.h1.*.nc' for v in variables},
    'ne256 rrtmg' : {v: '*.h1.*.nc' for v in variables},
    'ne256 rrtmgp': {v: '*.h1.*.nc' for v in variables},
    'ne256 lamlow': {v: '*.h1.*.nc' for v in variables},
    'ne4 rrtmg'   : {v: '*.h0.*.nc' for v in variables},
    'ne4 rrtmgp'  : {v: '*.h0.*.nc' for v in variables},
    # Map variables to files for SCREAM-HR production run.
    # TODO: make this more flexible
    'Run 1'         : {v: '*.eam.h[0-9]*.nc' for v in variables},
    'Run 2'         : {v: 'SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h[0-9]*.nc' for v in variables},
    'SCREAMv0.1 DY2': {v: '*.eam.h[0-9]*.nc' for v in variables},
    'SCREAMv1.1_firstday'      : {v: 'rgr_output.scream.*.nc' for v in variables},
}

# Overwrite some glob strings
glob_strings['SCREAMv0.1 DY2']['PS'] = '*.eam.h4.*.nc'

# Time offsets; hack so we can compare F-cases.
# HOW TO USE THIS FEATURE: in the event that you want to compare a run for a
# different time period than what we have observations for, you can use this to
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
do_contour_maps = True
do_time_series  = False
do_zonal_means  = False
do_anomalies    = False

#
# Make anomaly plots
#
if do_anomalies:

    # Find files
    files = [sorted(glob(f'{data_paths[c]}/{glob_strings[c]}')) for c in (test_case_name, *cntl_case_names)]
    print(files)
    raise RuntimeError('Anomalies not implemented yet')

#
# Make contour maps
#
if do_contour_maps:
    for vname in variables:
        print(vname)
        # For map plots we only compare two at a time so we can look at differences
        for cntl_case_name in cntl_case_names:
            print(cntl_case_name)

            # Find files and make sure we can retrieve variable
            if vname not in glob_strings[test_case_name] or vname not in glob_strings[cntl_case_name]: continue
            test_files = sorted(glob(f'{data_paths[test_case_name]}/{glob_strings[test_case_name][vname]}'))
            cntl_files = sorted(glob(f'{data_paths[cntl_case_name]}/{glob_strings[cntl_case_name][vname]}'))
            test_files = [f for f in test_files if can_retrieve_field(f, vname)]
            cntl_files = [f for f in cntl_files if can_retrieve_field(f, vname)]
            if len(test_files) == 0 or len(cntl_files) == 0: continue

            # Figure out what mapping files we need based on test and obs data
            if get_grid_name(test_files[0]) != get_grid_name(cntl_files[0]):
                map_file = get_mapping_file(test_files[0], cntl_files[0], mapping_root, method='nearestdtos')
            else:
                map_file = None

            # Open datasets and subset?

            # Compare maps
            print(f'Making {test_case_name} vs {cntl_case_name} contour maps for {vname}...')
            figname = f'{graphics_root}/contour_maps/{vname}_{test_case_name.replace(" ", "_")}_vs_{cntl_case_name.replace(" ", "_")}_nearestdtos_maps.png'
            os.makedirs(os.path.dirname(figname), exist_ok=True)
            compare_maps.main(
                vname, figname, test_files, cntl_files, 
                maps=(map_file, None), labels=(test_case_name, cntl_case_name), 
                verbose=False, t1=t1, t2=t2,
            )

#
# Make zonal mean plots
#
if do_zonal_means:
    for vname in variables:

        # Find files, only those that contain variable we are interested in
        files = []
        names = []
        for n in (test_case_name, *cntl_case_names):
            if vname not in glob_strings[n]: continue
            files_with_var = [f for f in sorted(glob(f'{data_paths[n]}/{glob_strings[n][vname]}')) if can_retrieve_field(f, vname)]
            if len(files_with_var) > 0:
                files.append(files_with_var)
                names.append(n)
        if len(files) == 0: continue

        # Figure out what mapping files we need based on test and obs data
        default_grid_file = sorted(glob(f'{data_paths["CERES-SYN"]}/*.nc'))[0]
        maps = [get_mapping_file(f[0], default_grid_file, mapping_root)
                if not is_latlon(f[0]) else None for f in files]

        # Make zonal mean plots
        print(f'Making zonal mean plots for {vname}...')
        figname = f'{graphics_root}/{vname}_{test_case_name.replace(" ", "_")}_vs_all_zonal.png'
        os.makedirs(os.path.dirname(figname), exist_ok=True)
        compare_zonal_means.main(files, names, vname, figname,
                maps=maps,
                time_offsets=[time_offsets[n] for n in names],
                )

#
# Make time series plots
#
if do_time_series:
    for vname in variables:

        # Find files, only those that contain variable we are interested in
        files = []
        names = []
        for n in (test_case_name, *cntl_case_names):
            if vname not in glob_strings[n]: continue
            files_with_var = [f for f in sorted(glob(f'{data_paths[n]}/{glob_strings[n][vname]}')) if can_retrieve_field(f, vname)]
            if len(files_with_var) > 0:
                files.append(files_with_var)
                names.append(n)

        # Make time series plots
        print(f'Making time series for {vname}...')
        figname = f'{graphics_root}/time_series/{vname}_{test_case_name.replace(" ", "_")}_vs_all_timeseries.png'
        os.makedirs(os.path.dirname(figname), exist_ok=True)
        compare_time_series.main(files, vname, figname, names, time_offsets=[time_offsets[n] for n in names], verbose=True)
