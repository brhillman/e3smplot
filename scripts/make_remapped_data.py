#!/usr/bin/env python

from glob import glob
from subprocess import run, PIPE, STDOUT
import os

#native_dir = '/lus/theta-fs0/projects/ClimateEnergy_4/ndkeen/SCREAMv010.SCREAM-DY2.ne1024pg2.bigrid.oct27.n1024a16x8.dt100.ps.shoc1.MPASSI/run'
native_dir = '/lus/theta-fs0/projects/ClimateEnergy_4/wlin/20211104.V01.SCREAM-DY1.ne1024pg2.bigrid.n1024a16x8.dt100.ps.shoc1.MPASSI.invdist-1/run'
remapped_dir = '/lus/theta-fs0/projects/ClimateEnergy_4/brhillman/scream/SCREAMv0.1-DY1/remapped'
mapping_dir = '/lus/theta-fs0/projects/ClimateEnergy_4/brhillman/mapdata'
mapfile = f'{mapping_dir}/map_ne1024pg2_to_fv256x512_nco.20201201.nc'
grid_name = 'fv256x512'

os.makedirs(remapped_dir, exist_ok=True)

# Remap everything in inputdir
native_files = sorted(glob(f'{native_dir}/*.eam.h[0-9].*.nc'))
remapped_files = [f'{remapped_dir}/{os.path.basename(os.path.splitext(f)[0])}.{grid_name}.nc' for f in native_files]
for native_file, remapped_file in zip(native_files, remapped_files):
    if os.path.exists(remapped_file): continue
    print(f'{native_file} -> {remapped_file}')
    p = run(f'./remap.sh {mapfile} {native_file} {remapped_file}'.split(' '), capture_output=True, text=True)
    print(p.stdout)
    p.check_returncode()
