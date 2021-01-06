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
    ds = open_dataset(sorted(these_files), chunks={'time': 1})
    
    # Get data
    print('get data')
    data = get_data(ds, v)

    # Compute area averages
    print('compute averages')
    w = get_area_weights(ds)
    m = area_average(data, w, dims=[d for d in data.dims if d != 'time'])
    
    # Convert units
    print('Convert units')
    if data.attrs['units'] == 'none':
        m = 100.0 * m
        m.attrs['units'] = '%'
        
    # Plot timeseries
    print('plot')
    pl = ax.plot(data.time, m, label=v)
    
    ds.close()
    
lh = pyplot.legend(pl)
