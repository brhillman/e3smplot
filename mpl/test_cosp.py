#!/usr/bin/env python3
import os
import sys
from glob import glob
from compare_maps import main as compare_maps
from compare_cfads import main as compare_cfads

varnames = (
    #'CLDTOT'        , 'CLDLOW'        , 'CLDMED'        , 'CLDHGH'        ,
    #'CLDTOT_CAL'    , 'CLDLOW_CAL'    , 'CLDMED_CAL'    , 'CLDHGH_CAL'    ,
    #'CLDTOT_CAL_ICE', 'CLDLOW_CAL_ICE', 'CLDMED_CAL_ICE', 'CLDHGH_CAL_ICE',
    #'CLDTOT_CAL_LIQ', 'CLDLOW_CAL_LIQ', 'CLDMED_CAL_LIQ', 'CLDHGH_CAL_LIQ', 
    #'CLDTOT_CAL_UN' , 'CLDLOW_CAL_UN' , 'CLDMED_CAL_UN' , 'CLDHGH_CAL_UN' ,
    #'CLDTOT_CS', 'CLDTOT_CS2', 'CLDTOT_CALCS', 'CLDTOT_ISCCP',
    #'MEANCLDALB_ISCCP', 'MEANPTOP_ISCCP', 'MEANTAU_ISCCP', 
    #'MEANTBCLR_ISCCP', 'MEANTB_ISCCP',
    #'CLTMODIS', 'CLLMODIS', 'CLMMODIS', 'CLHMODIS',
    #'CLIMODIS', 'CLWMODIS',
    #'IWPMODIS', 'LWPMODIS', 'PCTMODIS', 'REFFCLIMODIS', 'REFFCLWMODIS',
    #'TAUILOGMODIS', 'TAUWLOGMODIS', 'TAUTLOGMODIS',
    #'TAUIMODIS', 'TAUTMODIS', 'TAUWMODIS', 
    # These need to be derived
    #'CLDTOT_MISR', 'CLDLOW_MISR', 'CLDMED_MISR', 'CLDHGH_MISR',
    'CLDTOT_ISCCP', 'CLDLOW_ISCCP', 'CLDMED_ISCCP', 'CLDHGH_ISCCP',
)

cases = ('RRTMGP', 'RRTMG')
files = {
    'RRTMGP': glob('/global/cscratch1/sd/bhillma/e3sm_scratch/cori-knl/SMS_Lm1_P64x2.ne4_ne4.FC5AV1C-L.cori-knl_intel.cam-rrtmgp.cosp/run/*.cam.h0.*.nc')[0],
    'RRTMG': glob('/global/cscratch1/sd/bhillma/e3sm_scratch/cori-knl/SMS_Lm1_P64x2.ne4_ne4.FC5AV1C-L.cori-knl_intel.cosp/run/*.cam.h0.*.nc')[0],
}

plot_root = './cosp_plots'
os.makedirs(plot_root, exist_ok=True)
if False:
    for v in varnames:
        print(f'Plot variable {v}...'); sys.stdout.flush()
        compare_maps(v, os.path.join(plot_root, f'{v}_{cases[0]}_vs_{cases[1]}.png'), files[cases[0]], files[cases[1]], plot_method='regrid', labels=f'{cases[0]}, {cases[1]}')

cfad_vars = (
    'CFAD_DBZE94_CS', 'CFAD_SR532_CAL',
)
for v in cfad_vars:
    print(f'Plot variable {v}...'); sys.stdout.flush()
    compare_cfads(v, os.path.join(plot_root, f'{v}_{cases[0]}_vs_{cases[1]}.png'), files[cases[0]], files[cases[1]], labels=f'{cases[0]}, {cases[1]}')

