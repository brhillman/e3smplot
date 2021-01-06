#!/usr/bin/env python

from e3smplot.e3sm_utils import get_data, open_dataset, can_retrieve_field

# Compute a cloud mask based on a threshold of liquid and ice water content
def get_liq_cloud_mask(ds, threshold=1e-5):
    cldliq = get_data(ds, 'CLDLIQ')
    cld_mask = cldliq.where(cldliq > threshold).notnull() #(cldliq > threshold)
    #cld_mask = cldliq.copy()
    cld_mask.attrs = {
        'long_name': 'Liquid cloud mask',
        'units': 'none',
        'description': f'CLDLIQ > {threshold}',
    }
    return cld_mask

def get_ice_cloud_mask(ds, threshold=1e-5):
    cldice = get_data(ds, 'CLDICE')
    cld_mask = cldice.where(cldice > threshold).notnull() #(cldice > threshold)
    cld_mask.attrs = {
        'long_name': 'Ice cloud mask',
        'units': 'none',
        'description': f'CLDICE > {threshold}',
    }
    return cld_mask

def get_cloud_mask(ds):
    liq_mask = get_liq_cloud_mask(ds)
    ice_mask = get_ice_cloud_mask(ds)
    cld_mask = (liq_mask * ice_mask).notnull() #((liq_mask > 0) | (ice_mask > 0))
    cld_mask.attrs = {
        'long_name': 'Cloud mask',
        'units': 'none',
        'description': f'{liq_mask.attrs["description"]} | {ice_mask.attrs["description"]}',
    }
    return cld_mask

from matplotlib import pyplot
from glob import glob
from e3smplot.e3sm_utils import open_dataset, get_data, get_area_weights, area_average, can_retrieve_field
# Compare time series of CLDTOT, LIQ_MASK, ICE_MASK, CLD_MASK
variable_names = ('LIQ_MASK',)#'CLDTOT', 'LIQ_MASK', 'ICE_MASK', 'CLD_MASK')
output_path = '/global/cfs/cdirs/e3sm/terai/SCREAM/DYAMOND2/Output/20201127'
all_files = glob(f'{output_path}/*.eam.h[0-9].*.nc')
figure, ax = pyplot.subplots(1, 1) # figsize=(10, 10))

for ivar, v in enumerate(variable_names):
    # find files
    print(v)
    these_files = [f for f in all_files if can_retrieve_field(f, v)]
    
    # Load files
    print('load dataset')
    ds = open_dataset(sorted(these_files), chunks={'time': 1, 'lev': 1})
    
    # Get data
    print('get data')
    data = get_data(ds, v)

    # Compute area averages
    print('compute averages')
    #w = get_area_weights(ds)
    #m = area_average(data, w, dims=[d for d in data.dims if d != 'time'])
    m = data.mean(dim=[d for d in data.dims if d != 'time'])

    # Convert units
    if False: #data.attrs['units'] == 'none':
        print('Convert units')
        m = 100.0 * m
        m.attrs['units'] = '%'
        
    # Plot timeseries
    print('print values')
    print(m.values)
    #print('plot')
    #pl = ax.plot(data.time, m, label=v)
    
    ds.close()
