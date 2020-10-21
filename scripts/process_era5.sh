#!/bin/bash

# Location of original and reprocessed data
original_data_root="/global/cfs/cdirs/e3sm/terai/Obs_datasets/ERA5"
processed_data_root="/global/cfs/cdirs/e3sm/bhillma/obs_datasets/ERA5"
mkdir -p ${processed_data_root}

# Fix coordinate variables
for f in ${original_data_root}/*.nc; do
    echo "Reverse latitude order for ${f}..."
    ncpdq -a -latitude ${f} ${processed_data_root}/`basename ${f}`
done
echo "Done."
