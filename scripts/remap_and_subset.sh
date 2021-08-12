#!/bin/bash

module load python cmem
source activate e3smplot

export PATH=~zender/bin_cori:${PATH}

hpss_root="/home/t/terai/Production_runs/SCREAMv0.SCREAM-DY2.ne1024pg2.20201127"
output_root="/global/cscratch1/sd/bhillma/scream/dyamond2/remap/3072x6144"
map_file="/global/cfs/cdirs/e3sm/bhillma/maps/map_ne1024pg2_to_3072x6144_nco_con_20210616.nc"
do_remap=1
do_subset=0
do_remove_native=1
do_remove_derived=1
do_remove_remapped=0

# Hard-code these because I'm not sure how to get them from an hsi command
files=(
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-01-20-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-01-21-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-01-22-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-01-23-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-01-24-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-01-25-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-01-26-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-01-27-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-01-28-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-01-29-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-01-30-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-01-31-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-01-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-02-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-03-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-04-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-05-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-06-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-07-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-08-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-09-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-10-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-11-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-12-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-13-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-14-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-15-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-16-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-17-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-18-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-19-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-20-00000.nc
	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-21-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-22-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-23-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-24-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-25-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-26-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-27-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-02-28-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-03-01-00000.nc
 	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.2020-03-02-00000.nc
)

# Variables to pull
varnames=("PSL")
cmd="" #"echo" #"srun --nodes=1 --constraint=amd --time=01:00:00 --partition=bigmem"
# Loop over files
mkdir -p ${output_root}
for f in ${files[@]}; do
    echo "Working on ${f}..."

    # Grab a single file from HPSS to operate on
    tmp_file=${output_root}/${f}
    hpss_file=${hpss_root}/${f}

    # Work on one variable at a time
    for v in ${varnames[@]}; do

        derived_file=`echo ${tmp_file} | sed "s/eam.h4/eam.${v}/"`
        remapped_file=`echo ${derived_file} | sed "s/.nc/.remapped.nc/"`
        subset_file=`echo ${derived_file} | sed "s/.nc/.remapped_and_subset.nc/"`

        # If we have already created this file, move on
        #if [ -e ${subset_file} ]; then
        if [ -e ${remapped_file} ]; then
            echo "  ${remapped_file} already exists, skipping..."
            continue
        fi

        # Grab file from HPSS if we do not have it locally
        if [ ! -e ${tmp_file} ]; then
            echo "  --Grabbing ${f} from hpss..."
            $cmd /usr/common/mss/bin/hsi "get ${tmp_file} : ${hpss_file}" || exit 1
        fi

        # Extract only the data we want
        if [ ! -e ${derived_file} ]; then
            echo "  --Extracting ${v}..."
            $cmd ncks -O -v lat,lon,area,${v} ${tmp_file} ${derived_file} || exit 1
        fi

        # Remap from native grid to something smaller to work with
        if [ ${do_remap} -eq 1 ] && [ ! -e ${remapped_file} ] ; then
            echo "  --Remapping ${v}..."
            #srun --nodes 1 --constraint amd --time 00:20:00 --partition shared ./remap.sh ${map_file} ${derived_file} ${remapped_file} || exit 1
            echo $cmd ncks --map ${map_file} ${derived_file} ${remapped_file} || exit 1
            $cmd ncks --map ${map_file} ${derived_file} ${remapped_file} || exit 1
        fi

        # Subset for Arctic
        if [ ${do_subset} -eq 1 ] && [ ! -e ${subset_file} ]; then
            echo "  --Subsetting for arctic..."
            $cmd ncks -O -d lat,40.0,90.0 ${remapped_file} ${subset_file} || exit 1
        fi

        # Remove temporary files for this variable
        if [ -e ${derived_file} ] && [ ${do_remove_derived} -eq 1 ]; then rm -v ${derived_file}; fi
        if [ -e ${remapped_file} ] && [ ${do_remove_remapped} -eq 1 ]; then rm -v ${remapped_file}; fi

    done

    # Remove temporary file containing all variables
    if [ -e ${tmp_file} ] && [ ${do_remove_native} -eq 1 ]; then rm -v ${tmp_file}; fi

done
