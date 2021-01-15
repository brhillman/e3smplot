#!/usr/bin/env python

import xarray
import numpy
import time 
import sys
from e3smplot.e3sm_utils import get_data

# Compute a cloud mask based on a threshold of liquid and ice water content
def get_liq_cld_mask(ds, threshold=1e-5):
    cldliq = get_data(ds, 'CLDLIQ')
    cld_mask = (cldliq > threshold) #(cldliq > threshold)
    cld_mask.attrs = {
        'name': 'liq_cld_mask',
        'long_name': 'Liquid cloud mask',
        'units': 'none',
        'description': f'CLDLIQ > {threshold}',
    }
    return cld_mask

def get_ice_cld_mask(ds, threshold=1e-5):
    cldice = get_data(ds, 'CLDICE')
    cld_mask = (cldice > threshold) #(cldice > threshold)
    cld_mask.attrs = {
        'name': 'ice_cld_mask',
        'long_name': 'Ice cloud mask',
        'units': 'none',
        'description': f'CLDICE > {threshold}',
    }
    return cld_mask

def get_tot_cld_mask(ds):
    liq_mask = get_liq_cld_mask(ds)
    ice_mask = get_ice_cld_mask(ds)
    cld_mask = (liq_mask + ice_mask) > 0 #((liq_mask > 0) | (ice_mask > 0))
    cld_mask.attrs = {
        'name': 'tot_cld_mask',
        'long_name': 'Cloud mask',
        'units': 'none',
        'description': f'{liq_mask.attrs["description"]} | {ice_mask.attrs["description"]}',
    }
    return cld_mask

# Vertically-projected cloud area (shadow)
def get_liq_cld_area(ds):
    # First get 3D cloud mask
    cld_mask = get_liq_cld_mask(ds)
    # Project down
    cld_area = (cld_mask > 0).any(dim='lev')
    cld_area.attrs = {
        'name': 'liq_cld_area',
        'long_name': 'Liquid cloud area mask',
        'units': 'none',
        'description': 'any(cld_mask > 0)',
    }
    return cld_area

# Vertically-projected cloud area (shadow)
def get_ice_cld_area(ds):
    # First get 3D cloud mask
    cld_mask = get_ice_cld_mask(ds)
    # Project down
    cld_area = (cld_mask > 0).any(dim='lev')
    cld_area.attrs = {
        'name': 'ice_cld_area',
        'long_name': 'Ice cloud area mask',
        'units': 'none',
        'description': 'any(cld_mask > 0)',
    }
    return cld_area

# Vertically-projected cloud area (shadow)
def get_tot_cld_area(ds):
    # First get 3D cloud mask
    cld_mask = get_tot_cld_mask(ds)
    # Project down
    cld_area = (cld_mask > 0).any(dim='lev')
    cld_area.attrs = {
        'name': 'tot_cld_area',
        'long_name': 'Cloud area mask',
        'units': 'none',
        'description': 'any(cld_mask > 0)',
    }
    return cld_area

def main(inputfile, outputfile, **kwargs):

    # Open dataset
    ds = xarray.open_mfdataset(
        inputfile, drop_variables=('P3_input_dim', 'P3_output_dim'),
        chunks={'time': 1, 'lev': 1}, engine='netcdf4',
    )

    # Compute cloud masks
    print('Compute cloud masks...', end=''); sys.stdout.flush()
    t1 = time.perf_counter()
    liq_cld_mask = get_liq_cld_mask(ds)#.astype(numpy.float32).to_netcdf(outputfile) #.close()
    ice_cld_mask = get_ice_cld_mask(ds)#.astype(numpy.float32).to_netcdf(outputfile).close()
    tot_cld_mask = get_tot_cld_mask(ds)#.astype(numpy.float32).to_netcdf(outputfile).close()
    liq_cld_area = get_liq_cld_area(ds)#.astype(numpy.float32).to_netcdf(outputfile).close()
    ice_cld_area = get_ice_cld_area(ds)#.astype(numpy.float32).to_netcdf(outputfile).close()
    tot_cld_area = get_tot_cld_area(ds)#.astype(numpy.float32).to_netcdf(outputfile).close()
    t2 = time.perf_counter()
    print('done; elapsed time = {}'.format(t2 - t1))

    # Create a new dataset with output fields
    print('Create dataset...', end=''); sys.stdout.flush()
    t1 = time.perf_counter()
    ds_out = xarray.Dataset({
        'liq_cld_mask': liq_cld_mask.astype(numpy.float32),
        'ice_cld_mask': ice_cld_mask.astype(numpy.float32),
        'tot_cld_mask': tot_cld_mask.astype(numpy.float32),
        'liq_cld_area': liq_cld_area.astype(numpy.float32),
        'ice_cld_area': ice_cld_area.astype(numpy.float32),
        'tot_cld_area': tot_cld_area.astype(numpy.float32),
    })
    t2 = time.perf_counter()
    print('done; elapsed time = {}'.format(t2 - t1))

    # Write the new dataset to file
    print('Create output netcdf...', end=''); sys.stdout.flush()
    t1 = time.perf_counter()
    ds_out.to_netcdf(outputfile)
    t2 = time.perf_counter()
    print('done; elapsed time = {}'.format(t2 - t1))

    # Clean up
    print('Close datasets...', end=''); sys.stdout.flush()
    t1 = time.perf_counter()
    ds.close()
    ds_out.close()
    t2 = time.perf_counter()
    print('done; elapsed time = {}'.format(t2 - t1))

if __name__ == '__main__':
    import plac; plac.call(main)
