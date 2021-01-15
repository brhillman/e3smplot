#!/bin/bash

hpss_root="/home/t/terai/Production_runs/SCREAMv0.SCREAM-DY2.ne1024pg2.20201127"
output_root="${CSCRATCH}/scream/dyamond2/post-processed"
map_file="${CSCRATCH}/scream/dyamond2/map_ne1024pg2_to_fv256x512_nco.20201201.nc"

# Hard-code these because I'm not sure how to get them from an hsi command
files=(
	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-01-20-00000.nc
	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-01-21-00000.nc
	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-01-22-00000.nc
	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-01-23-00000.nc
	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-01-24-00000.nc
	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-01-25-00000.nc
	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-01-26-00000.nc
	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-01-27-00000.nc
	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-01-28-00000.nc
	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-01-29-00000.nc
	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-01-30-00000.nc
	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-01-31-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-01-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-02-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-03-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-04-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-05-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-06-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-07-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-08-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-09-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-10-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-11-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-12-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-13-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-14-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-15-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-16-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-17-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-18-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-19-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-20-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-21-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-22-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-23-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-24-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-25-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-26-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-27-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-28-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-03-01-00000.nc
#	SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-03-02-00000.nc
)

# Variables to compute
varnames=("LIQ_CLD_MASK" "ICE_CLD_MASK" "TOT_CLD_MASK" "LIQ_CLD_AREA" "ICE_CLD_AREA" "TOT_CLD_AREA")

# Loop over files
mkdir -p ${output_root}
for f in ${files[@]}; do
    echo "Working on ${f}..."

    # Grab a single file from HPSS to operate on
    tmp_file=${output_root}/${f}
    hpss_file=${hpss_root}/${f}

    # Compute cloud masks and save to new file
    for v in ${varnames[@]}; do

        derived_file=`echo ${tmp_file} | sed "s/eam.h7/eam.${v}/"`
        remapped_file=`echo ${derived_file} | sed "s/.nc/.remapped.nc/"`

        # If we have already created this file, move on
        if [ -e ${remapped_file} ]; then
            continue
        fi

        # Grab file from HPSS if we do not have it locally
        if [ ! -e ${tmp_file} ]; then
            echo "  --Grabbing ${f} from hpss..."
            /usr/common/mss/bin/hsi "get ${tmp_file} : ${hpss_file}" || exit 1
        fi

        # Compute derived fields and save to temporary file
        if [ ! -e ${derived_file} ]; then
            echo "  --Deriving ${v}..."
            srun --nodes 1 --constraint haswell --time 00:30:00 ./compute_derived_field.py ${v} ${tmp_file} ${derived_file} || exit 1
        fi

        # Remap from native grid to something smaller to work with
        if [ ! -e ${remapped_file} ]; then
            echo "  --Remapping ${v}..."
            srun --nodes 1 --constraint haswell --time 00:30:00 ./remap.sh ${map_file} ${derived_file} ${remapped_file} || exit 1
        fi

        # Remove native grid derived file
        if [ $? -eq 0 ]; then
            if [ -e ${derived_file} ]; then rm ${derived_file}; fi
        fi
    done

    # Remove temporary copy of native grid file
    if [ -e ${tmp_file} ]; then rm -v ${tmp_file}; fi

done
