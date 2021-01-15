#!/bin/bash

function ncdmnsz { ncks --trd -m -M ${2} | grep -E -i ": ${1}, size =" | cut -f 7 -d ' ' | uniq ; }

mapfile=$1
inputfile=$2
outputfile=$3

cmd=""

# Loop over times; remap each separately and then stitch together
ntime=`ncdmnsz "time" $inputfile`
echo "Loop over $ntime time samples..."
outputfiles=""
for i in `seq $ntime`; do

    # Extract one time slice; unfortunately it does not look like we can combine
    # with the below remap command so split up first
    inputfile_tmp=`dirname ${outputfile}`/`basename ${inputfile}`.${i}
    $cmd ncks -d time,$(expr $i - 1) ${inputfile} ${inputfile_tmp}

    # Compute masks for single time sample
    outputfile_tmp=${outputfile}.${i}
    #$cmd ncks --map ${mapfile} ${inputfile_tmp} ${outputfile_tmp}
    $cmd ncap2 -s "LIQ_CLD_MASK = (CLDLIQ > 1e-5);" ${inputfile_tmp} ${outputfile_tmp}

    # Fix record dimension
    $cmd ncks -O --mk_rec_dmn time ${outputfile_tmp} ${outputfile_tmp}.rec && mv ${outputfile_tmp}.rec ${outputfile_tmp}

    # Append to list of output files to concatenate
    outputfiles="${outputfiles} ${outputfile_tmp}"

    # Clean up tmp files
    $cmd rm ${inputfile_tmp}

done

# Stitch together outputfiles now
echo "Stitch together outputfiles..."
$cmd ncrcat ${outputfiles} ${outputfile}
$cmd rm ${outputfiles}

echo "Done."
exit 0
