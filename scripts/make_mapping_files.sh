#!/bin/bash

export PATH=~zender/bin_cori:$PATH

# Output root
output_root=${SCRATCH}/grids
mkdir -p ${output_root}

# Model grid
model_name="ne1024pg2"
model_grid="${output_root}/ne1024pg2_scrip_fix.nc"
echo "Generate model_grid file ${model_grid}"
if [ ! -e ${model_grid} ]; then
    GenerateCSMesh --alt --res 1024 --file ${output_root}/ne1024.g
    GenerateVolumetricMesh --in ${output_root}/ne1024.g --out ${output_root}/ne1024pg2.g --np 2 --uniform
    ConvertExodusToSCRIP --in ${output_root}/ne1024pg2.g --out ${model_grid}
fi

# Obs grids
obs_name="ceres1deg"
obs_file="/global/cfs/cdirs/e3sm/terai/Obs_datasets/CERES/CERES_SYN1deg-3H_Terra-Aqua-MODIS_Ed4.1_Subset_20200101-20200229.nc"
obs_grid="${output_root}/ceres_scrip.nc"
echo "Generate obs_grid file ${obs_grid}"
if [ ! -e ${obs_grid} ]; then
    ncks --rgr infer --rgr scrip=${obs_grid} ${obs_file} ~/foo.nc
fi

# Generate mapping file
map_file="${output_root}/map_${model_name}_to_${obs_name}_nco.nc"
if [ ! -e ${map_file} ]; then
    echo "Generate mapping file ${map_file}"
    ncremap --fl_fmt=64bit_data --alg_typ=nco --src_grd=${model_grid} --dst_grd=${obs_grid} -m ${map_file}
fi


obs_name="era5"
obs_file="/global/cfs/cdirs/e3sm/bhillma/obs_datasets/ERA5/ERA5_surf_20200101_20200229.nc"
obs_grid="${output_root}/era5_scrip_fix.nc"
echo "Generate obs_grid file ${obs_grid}"
if [ ! -e ${obs_grid} ]; then
    ncks --rgr infer --rgr scrip=${obs_grid} ${obs_file} ~/foo.nc
fi

method="nco"
map_file="${output_root}/map_${model_name}_to_${obs_name}_${method}.nc"
if [ ! -e ${map_file} ]; then
    echo "Generate mapping file ${map_file}"
    ncremap --fl_fmt=64bit_data --alg_typ=${method} --src_grd=${model_grid} --dst_grd=${obs_grid} -m ${map_file}
fi
