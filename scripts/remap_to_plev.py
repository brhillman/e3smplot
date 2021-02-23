#!/usr/bin/env python

import numpy as np
import xarray
import os
from glob import glob
from e3smplot.utils import interp_along_axis, update_progress
from e3smplot.e3sm_utils import get_data

output_root = '/global/cfs/cdirs/e3sm/terai/SCREAM/DYAMOND2/Output/20201127/regridded'
postprocessed = '/global/cfs/cdirs/e3sm/bhillma/scream/dyamond2/post-processed'
os.makedirs(postprocessed, exist_ok=True)

# Repeat for h6, h7, h8
ps_series = 'h4'
for data_series in ('h6', 'h7', 'h8'):

    print(f'Process files in series {data_series}...')
    data_files = sorted(glob(f'{output_root}/*.{data_series}.*.nc'))
    ps_files = [f.replace(f'.{data_series}.', f'.{ps_series}.') for f in data_files]
    levels_file = '/global/cfs/cdirs/e3sm/bhillma/obs_datasets/ERA5/ERA5_TandQ_20200201_20200229.nc'

    ds_l = xarray.open_dataset(levels_file, chunks={'time': 1})
    levels = ds_l.level.rename({'level': 'plev'})

    exclude_variables = ('hyam', 'hybm', 'hyai', 'hybi', 'lev')

    for ifile, (f_d, f_p) in enumerate(zip(data_files, ps_files)):

        # Update progress bar
        update_progress(ifile+1, len(data_files))

        # Open file
        ds_v = xarray.open_dataset(f_d)
        ds_p = xarray.open_dataset(f_p)

        # Empty dict to hold remapped variables
        ds_out = xarray.Dataset({})

        # Loop over variables in file
        for v in ds_v.variables.keys():

            # If this variable does not have a level dimension, move on
            if 'lev' not in ds_v.variables[v].dims or v in exclude_variables: continue

            # Read data
            d = get_data(ds_v, v)
            p = get_data(ds_p, 'PMID')

            # Do vertical remap
            axis = 1
            shape_out = list(d.shape)
            shape_out[axis] = len(levels)
            dims_out = list(d.dims)
            dims_out[axis] = 'plev'
            coords_out = {c: (levels if c == 'plev' else d.coords[c]) for c in dims_out}
            d_on_levels = xarray.DataArray(
                np.empty(shape_out),
                name=v,
                dims=dims_out,
                coords=coords_out,
                attrs=d.attrs,
            )
            # Do linear interpolation in log(pressure)
            d_on_levels.data = interp_along_axis(
                    np.log(levels.data), np.log(p.transpose(*d.dims).data), d.data,
                    axis=axis
            )
            ds_out[v] = d_on_levels.astype(np.float32)

        # Save file
        f_out = f'{postprocessed}/{os.path.splitext(os.path.basename(f_d))[0]}.plev.nc'
        ds_out.to_netcdf(f_out)
        ds_v.close()
        ds_p.close()
