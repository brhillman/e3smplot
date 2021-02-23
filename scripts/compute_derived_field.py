#!/usr/bin/env python

import xarray
import numpy
from e3smplot.e3sm_utils import get_data

def main(vname, inputfile, outputfile, **kwargs):

    # Process kwargs
    for k, v in kwargs.items(): kwargs[k] = eval(v)

    # Open dataset
    ds = xarray.open_mfdataset(
        inputfile, drop_variables=('P3_input_dim', 'P3_output_dim'),
        **kwargs
    )

    # Compute cloud masks and save to disk
    da = get_data(ds, vname)
    ds_out = xarray.Dataset({vname: da})
    ds_out.to_netcdf(outputfile, encoding={vname: {'_FillValue': -9999}})

    # Clean up
    ds_out.close()

if __name__ == '__main__':
    import plac; plac.call(main)
