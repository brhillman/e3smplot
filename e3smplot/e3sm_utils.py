#!/usr/bin/env python3

import numpy
import xarray

def get_pressure(dataset, interfaces=False):
    if interfaces:
        a = dataset['hyai']
        b = dataset['hybi']
    else:
        a = dataset['hyam']
        b = dataset['hybm']
    p0 = dataset['P0']
    ps = dataset['PS']
    
    pressure = a * p0 + b * ps
    
    pressure = pressure * 1e-2
    pressure.attrs['units'] = 'hPa'
    pressure.attrs['long_name'] = 'Pressure'
    
    return pressure

def open_dataset(files, **kwargs):

    # Open dataset as a dask array
    ds = xarray.open_mfdataset(
        sorted(files), combine='by_coords',
        drop_variables=('P3_output_dim', 'P3_input_dim'), **kwargs
    ) 

    # Rename coordinate variables
    if 'latitude' in ds.dims: ds = ds.rename({'latitude': 'lat'})
    if 'longitude' in ds.dims: ds = ds.rename({'longitude': 'lon'})

    # Fix sensible heat flux
    # TODO: should not need to do this
    if 'msshf' in ds.variables.keys():
        ds['msshf'].values = -ds['msshf'].values

    # Convert cftime coordinate
    if isinstance(ds.time.values[0], cftime._cftime.DatetimeNoLeap):
        try:
            ds['time'] = ds.indexes['time'].to_datetimeindex()
        except:
            print('Could not convert times to datetimeindex. But proceeding anyway...')

    # Return dataset
    return ds


# Define function to read specialized data from E3SM files
def get_data(dataset, field):
    data = None
    if field in dataset.variables.keys():
        data = dataset[field]
    elif field == 'TSI':
        solar_insolation = get_data(dataset, 'SOLIN')
        cosine_solar_zenith = get_data(dataset, 'COSZRS')
        data = solar_insolation / cosine_solar_zenith
        data.attrs['long_name'] = 'Total solar irradiance'
        data.attrs['units'] = solar_insolation.attrs['units']
    elif field == 'PMID':
        data = get_pressure(dataset)
    elif field == 'PINT':
        data = get_pressure(dataset, interfaces=True)
    elif field == 'FNS':
        flux_up = get_data(dataset, 'FUS')
        flux_dn = get_data(dataset, 'FDS')
        data = flux_dn - flux_up
        data.attrs['long_name'] = 'Net SW flux'
        data.attrs['units'] = flux_up.attrs['units']
    elif field == 'FNSC':
        flux_up = get_data(dataset, 'FUSC')
        flux_dn = get_data(dataset, 'FDSC')
        data = flux_dn - flux_up
        data.attrs['long_name'] = 'Net clearky SW flux'
        data.attrs['units'] = flux_up.attrs['units']
    elif field == 'FNL':
        flux_up = get_data(dataset, 'FUL')
        flux_dn = get_data(dataset, 'FDL')
        data = flux_up - flux_dn
        data.attrs['long_name'] = 'Net LW flux'
        data.attrs['units'] = flux_up.attrs['units']
    elif field == 'FNLC':
        flux_up = get_data(dataset, 'FULC')
        flux_dn = get_data(dataset, 'FDLC')
        data = flux_up - flux_dn
        data.attrs['long_name'] = 'Net clearky LW flux'
        data.attrs['units'] = flux_up.attrs['units']
    elif field == 'CRESW':
        clear_sky_flux = get_data(dataset, 'FNSC')
        all_sky_flux = get_data(dataset, 'FNS')
        data = all_sky_flux - clear_sky_flux
        data.attrs['long_name'] = 'Net SW cloud radiative effect'
        data.attrs['units'] = all_sky_flux.attrs['units']
    elif field == 'CRELW':
        clear_sky_flux = get_data(dataset, 'FNLC')
        all_sky_flux = get_data(dataset, 'FNL')
        data = all_sky_flux - clear_sky_flux
        data.attrs['long_name'] = 'Net LW cloud radiative effect'
        data.attrs['units'] = all_sky_flux.attrs['units']
    elif field == 'FSDT':
        flux_dn = get_data(dataset, 'FDS')
        data = flux_dn.isel(ilev=0)
        data.attrs['long_name'] = 'Downward SW flux at top of model'
        data.attrs['units'] = flux_dn.attrs['units']
    elif field == 'FSUT':
        flux_up = get_data(dataset, 'FUS')
        data = flux_up.isel(ilev=0)
        data.attrs['long_name'] = 'Upward SW flux at top of model'
        data.attrs['units'] = flux_up.attrs['units']
    elif field == 'FLDT':
        flux_dn = get_data(dataset, 'FDL')
        data = flux_dn.isel(ilev=0)
        data.attrs['long_name'] = 'Downward LW flux at top of model'
        data.attrs['units'] = flux_dn.attrs['units']
    elif field == 'FLUT':
        flux_up = get_data(dataset, 'FUL')
        data = flux_up.isel(ilev=0)
        data.attrs['long_name'] = 'Upward LW flux at top of model'
        data.attrs['units'] = flux_up.attrs['units']
    elif field == 'FSDS':
        flux_dn = get_data(dataset, 'FDS')
        data = flux_dn.isel(ilev=-1)
        data.attrs['long_name'] = 'Downward SW flux at surface'
        data.attrs['units'] = flux_dn.attrs['units']
    elif field == 'FSUS':
        flux_up = get_data(dataset, 'FUS')
        data = flux_up.isel(ilev=-1)
        data.attrs['long_name'] = 'Upward SW flux at surface'
        data.attrs['units'] = flux_up.attrs['units']
    elif field == 'FLDS':
        flux_dn = get_data(dataset, 'FDL')
        data = flux_dn.isel(ilev=-1)
        data.attrs['long_name'] = 'Downward LW flux at surface'
        data.attrs['units'] = flux_dn.attrs['units']
    elif field == 'FLUS':
        flux_up = get_data(dataset, 'FUL')
        data = flux_up.isel(ilev=-1)
        data.attrs['long_name'] = 'Upward LW flux at surface'
        data.attrs['units'] = flux_up.attrs['units']
    elif field == 'FSDTC':
        flux_dn = get_data(dataset, 'FDSC')
        data = flux_dn.isel(ilev=0)
        data.attrs['long_name'] = 'Downward SW clearsky flux at top of model'
        data.attrs['units'] = flux_dn.attrs['units']
    elif field == 'FSUTC':
        flux_up = get_data(dataset, 'FUSC')
        data = flux_up.isel(ilev=0)
        data.attrs['long_name'] = 'Upward SW clearsky flux at top of model'
        data.attrs['units'] = flux_up.attrs['units']
    elif field == 'FLDTC':
        flux_dn = get_data(dataset, 'FDLC')
        data = flux_dn.isel(ilev=0)
        data.attrs['long_name'] = 'Downward LW clearsky flux at top of model'
        data.attrs['units'] = flux_dn.attrs['units']
    elif field == 'FLUTC':
        flux_up = get_data(dataset, 'FULC')
        data = flux_up.isel(ilev=0)
        data.attrs['long_name'] = 'Upward LW clearsky flux at top of model'
        data.attrs['units'] = flux_up.attrs['units']
    elif field == 'FSDSC':
        flux_dn = get_data(dataset, 'FDSC')
        data = flux_dn.isel(ilev=-1)
        data.attrs['long_name'] = 'Downward SW clearsky flux at surface'
        data.attrs['units'] = flux_dn.attrs['units']
    elif field == 'FSUSC':
        flux_up = get_data(dataset, 'FUSC')
        data = flux_up.isel(ilev=-1)
        data.attrs['long_name'] = 'Upward SW clearsky flux at surface'
        data.attrs['units'] = flux_up.attrs['units']
    elif field == 'FLDSC':
        flux_dn = get_data(dataset, 'FDLC')
        data = flux_dn.isel(ilev=-1)
        data.attrs['long_name'] = 'Downward LW clearsky flux at surface'
        data.attrs['units'] = flux_dn.attrs['units']
    elif field == 'FLUSC':
        flux_up = get_data(dataset, 'FULC')
        data = flux_up.isel(ilev=-1)
        data.attrs['long_name'] = 'Upward LW clearsky flux at surface'
        data.attrs['units'] = flux_up.attrs['units']
    elif field == 'RESTOM':
        net_sw_flux = get_data(dataset, 'FSNT')
        net_lw_flux = get_data(dataset, 'FLNT')
        data = net_sw_flux - net_lw_flux
        data.attrs['long_name'] = 'Top of model net radiative flux'
        data.attrs['units'] = net_sw_flux.attrs['units']
        data.attrs['formula'] = 'FSNT - FLNT'
    elif field == 'FSNTOA':
        if 'toa_sw_all_3h' in dataset.variables.keys() and 'toa_solar_all_3h' in dataset.variables.keys():
            flux_up = dataset['toa_sw_all_3h']
            flux_dn = dataset['toa_solar_all_3h']
            data = flux_dn - flux_up
            data.attrs['long_name'] = 'Net TOA SW flux'
            data.attrs['units'] = 'W/m2'
    elif field == 'FSNTOAC':
        if 'toa_sw_clr_3h' in dataset.variables.keys() and 'toa_solar_all_3h' in dataset.variables.keys():
            flux_up = dataset['toa_sw_clr_3h']
            flux_dn = dataset['toa_solar_all_3h']
            data = flux_dn - flux_up
            data.attrs['long_name'] = 'Net TOA SW clearsky flux'
            data.attrs['units'] = 'W/m2'
    elif field == 'FLNT':
        if 'toa_lw_all_3h' in dataset.variables.keys():
            data = dataset['toa_lw_all_3h']
            data.attrs['long_name'] = 'Longwave flux at TOA'
    elif field == 'FLNTC':
        if 'toa_lw_clr_3h' in dataset.variables.keys():
            data = dataset['toa_lw_clr_3h']
            data.attrs['long_name'] = 'Longwave clear-sky flux at TOA'
    elif field == 'CLDTOT':
        # CERES-SYN version of CLDTOT
        if 'cldarea_total_3h' in dataset.variables.keys():
            data = dataset['cldarea_total_3h']
            data.attrs['long_name'] = 'Cloud area fraction'
            data.attrs['units'] = '%'
    elif field == 'CLDLIQICE':
        cldliq = get_data(dataset, 'CLDLIQ')
        cldice = get_data(dataset, 'CLDICE')
        data = cldliq + cldice
        data.attrs = cldliq.attrs
        data.attrs['long_name'] = 'Combined liq+ice condensate'
    elif field == 'PRECT':
        precc = get_data(dataset, 'PRECC')
        precl = get_data(dataset, 'PRECL')
        data = precc + precl
        data.attrs = precc.attrs
        data.attrs['long_name'] = 'Total precipitation rate'
    elif field == 'TGCLDWP':
        clwp = get_data(dataset, 'TGCLDLWP')
        ciwp = get_data(dataset, 'TGCLDIWP')
        data = clwp + ciwp
        data.attrs = clwp.attrs
        data.attrs['long_name'] = 'Total gridbox cloud water path'
        
    # CALIPSO-simulated or CALIPSO-retrieved fields
    elif field == 'cltcalipso':
        data = get_data(dataset, 'CLDTOT_CAL')
    elif field == 'cltcalipso_liq':
        data = get_data(dataset, 'CLDTOT_CAL_LIQ')
    elif field == 'cltcalipso_ice':
        data = get_data(dataset, 'CLDTOT_CAL_ICE')
    elif field == 'clcalipso':
        data = get_data(dataset, 'CLD_CAL')
    elif field == 'clcalipso_liq':
        data = get_data(dataset, 'CLD_CAL_LIQ')
    elif field == 'clcalipso_ice':
        data = get_data(dataset, 'CLD_CAL_ICE')
        
    elif field == 'CRM_QCLD':
        crm_qc = get_data(dataset, 'CRM_QC')
        crm_qi = get_data(dataset, 'CRM_QI')
        data = crm_qc + crm_qi
        data.attrs = crm_qc.attrs
        data.attrs['long_name'] = 'Total CRM cloud water (liq + ice)'
    elif field == 'CRM_QPRC':
        crm_qpc = get_data(dataset, 'CRM_QPC')
        crm_qpi = get_data(dataset, 'CRM_QPI')
        data = crm_qpc + crm_qpi
        data.attrs = crm_qpc.attrs
        data.attrs['long_name'] = 'Total CRM precip water (liq + ice)'
    elif field == 'CRM_QTOT':
        crm_qcld = get_data(dataset, 'CRM_QCLD')
        crm_qprc = get_data(dataset, 'CRM_QPRC')
        crm_qv = get_data(dataset, 'CRM_QV')
        data = crm_qcld + crm_qprc + crm_qv
        data.attrs = crm_qv.attrs
        data.attrs['long_name'] = 'Total CRM water (cld + prec + qv)'
    elif field == 'TREFHT':
        if 't2m' in dataset.variables.keys():
            data = get_data(dataset, 't2m')
    elif field == 'SHFLX':
        if 'msshf' in dataset.variables.keys():
            data = get_data(dataset, 'msshf')
            #data.values = -data.values
    elif field == 'longitude':
        if 'lon' in dataset.variables.keys():
            data = get_data(dataset, 'lon')
        elif 'xc' in dataset.variables.keys():
            data = get_data(dataset, 'xc')
    elif field == 'latitude':
        if 'lat' in dataset.variables.keys():
            data = get_data(dataset, 'lat')
        elif 'yc' in dataset.variables.keys():
            data = get_data(dataset, 'yc')
    elif field == 'TMQ':
        if 'tcwv' in dataset.variables.keys():
            data = get_data(dataset, 'tcwv')

    
    # Check if we were able to find or derive the requested field
    if data is None:
        raise ValueError(f'{field} not found in dataset and no variables found to derive') 

    # Adjust units if necessary
    if field in ('TGCLDLWP', 'TGCLDIWP', 'TGCLDWP'):
        if data.units == 'kg/m2':
            attrs = data.attrs
            data = 1e3 * data
            data.attrs = attrs
            data.attrs['units'] = 'g/m2'
    elif field in ('CLDTOT', 'cltcalipso', 'cltcalipso_liq', 'cltcalipso_ice',
                   'clcalipso',  'clcalipso_liq',  'clcalipso_ice'):
        if data.units.lower() in ('1', 'fraction', 'none', '1 fraction'):
            attrs = data.attrs
            data = 100 * data
            data.attrs = attrs
            data.attrs['units'] = '%'
    elif field in ('PRECC', 'PRECL', 'PRECT'):
        if data.attrs['units'].lower() == 'm/s':
            attrs = data.attrs
            data = 60 * 60 * 24 * 1e3 * data
            data.attrs = attrs
            data.attrs['units'] = 'mm/day'
    elif field in ('QRS', 'QRSC', 'QRL', 'QRLC'):
        if data.attrs['units'].lower() == 'k/s':
            attrs = data.attrs
            data = 60 * 60 * 24 * data
            data.attrs = attrs
            data.attrs['units'] = 'K/day'
        
    return data


def get_coords(ds_data, ds_grid=None):
    if ds_grid is not None:
        if 'lon' in ds_grid and 'lat' in ds_grid:
            x = ds_grid['lon']
            y = ds_grid['lat']
        elif 'grid_corner_lon' in ds_grid and 'grid_corner_lat' in ds_grid:
            x = ds_grid['grid_corner_lon']
            y = ds_grid['grid_corner_lat']
        else:
            raise RuntimeError('No valid coordinates in grid file.')
    else:
        x = get_data(ds_data, 'longitude')
        y = get_data(ds_data, 'latitude')
    return x, y


def get_area_weights(ds):
    # Get weights; either use pre-computed or cosine(latitude) weights
    if 'area' in ds.variables.keys():
        wgt = ds.area
    else:
        wgt = numpy.cos(ds.lat * numpy.pi / 180.0)
    return wgt


def read_files(*files, fix_time=False, year_offset=None, **kwargs):
    # read files without decoding timestamps, because control dates are
    # probably outside pandas valid range
    from xarray import open_mfdataset
    ds = open_mfdataset(*files, decode_times=False, autoclose=True, **kwargs)

    # fix time_bnds
    if fix_time and 'time_bnds' in ds.variables.keys():
        ds['time_bnds'].attrs['units'] = ds['time'].attrs['units']

        # for monthly averages, the time is defined as the *end* of the
        # averaging interval...i.e., for the 0001-01 monthly mean, the time
        # coordinate has a value of 0001-02-01 00:00:00. This is inconvenient
        # because if we want to select all the January means to then do 
        # something like calculate a climatology, we would need to offset these
        # time values to get the month that the averaging is actually over.
        # As a work-around, we can redefine the time coordinate values to be the
        # average of the time_bnds variable.
        #ds['time_bnds'].attrs['units'] = ds['time'].attrs['units']
        bnd_dim = list(set(ds['time_bnds'].dims) - set(('time',)))
        attrs = ds['time'].attrs
        ds['time'] = ds['time_bnds'].astype('int64').mean(bnd_dim).astype('datetime64[ns]')
        ds['time'].attrs = attrs
        
    # manually decode the times
    from xarray import decode_cf
    if year_offset is not None:
        # Add an artificial offset to the time coordinate to handle cases where
        # the time units are out of bounds. This is needed for, e.g., F-compset
        # simulations that set the initial year as 0001.
        # as a hackish solution, we can just adjust the units of the time
        # axis...first, find the current base year
        time_units = ds['time'].attrs['units'] 
        base_year_idx = time_units.index('since ') + 6
        year = int(time_units[base_year_idx:base_year_idx+4])

        # add specified offset to base year
        new_year = year + year_offset
        new_units = time_units.replace(
            'since %04i'%year, 'since %04i'%new_year
        )

        # update the time attributes
        ds['time'].attrs['units'] = new_units
        ds['time_bnds'].attrs['units'] = new_units
        ds['time'].attrs['units_note'] = 'added a %i year offset to units to allow decoding with pandas'%year_offset
     
            
    # Finally, decode the dataset
    ds = decode_cf(ds)
        
    # update our dataset with the fixed up dataset read in by this method
    return ds


def area_average(data, weights, dims=None):
    
    '''Calculate area-weighted average over dims.'''
      
    if dims is None: dims = data.dims
        
    # Need to broadcast weights to make sure they have the
    # same size/shape as data. For example, data is (lat, lon)
    # but we passed weights with shape (lat,), or data is
    # (time, ncol) but we passed weights with shape (ncol,).
    from xarray import broadcast
    weights, *__ = broadcast(weights, data)
    
    # Mask weights consistent with data so we do not miscount
    # missing columns
    weights = weights.where(data.notnull())
    
    # Do the averaging        
    data_mean = (weights * data).sum(dim=dims) / weights.sum(dim=dims)
    
    # Copy over attributes, which we lose in the averaging
    # calculation
    data_mean.attrs = data.attrs
    
    # Return averaged data
    return data_mean


# Define a function to mask a list of data arrays consistent with one another
def mask_consistent(data_arrays):
    '''Mask a list of xarray.DataArrays so that they are consistent in terms
       of missing data'''
    for iarray in range(len(data_arrays)):
        for jarray in range(len(data_arrays)):
            if iarray == jarray: continue
            data_arrays[iarray] = data_arrays[iarray].where(data_arrays[jarray].notnull())
                                                            
    return data_arrays


def regrid_data(x1, y1, x2, y2, data):

    # Interpolatee to new grid
    from scipy.interpolate import griddata
    new_data = griddata((x1, y1), data, (x2.values[None,:], y2.values[:,None]), method='linear')

    # Turn these into DataArrays
    from xarray import DataArray
    new_data = DataArray(
        new_data, dims=('lat', 'lon'),
        coords={'lon': x2, 'lat': y2},
        attrs=data.attrs,
    )

    # Return DataArrays
    return new_data


def calculate_zonal_mean(data, weights, old_lat, avg_dims=None, lat_edges=None, nlat=180):
    
    import numpy as numpy

    # Mask data weights
    weights, *__ = xarray.broadcast(weights, data)
    weights = weights.where(data.notnull())
    
    # Calculate new latitudes
    if lat_edges is None:
        lat_edges = numpy.linspace(-90, 90, nlat+1)

    # Figure out what dimensions we are going to average over
    if avg_dims is None: avg_dims = data.dims
        
    # Calculating zonal mean for each latitude by binning data according to lat values
    # TODO: generalize this to arbitrary dimensions and ordering
    lat_centers = numpy.zeros(len(lat_edges) - 1)
    lat_centers2 = numpy.array([(y1 + y2) / 2.0 for (y1, y2) in zip(lat_edges[:-1], lat_edges[1:])])
    if 'lev' in data.dims:
        data_zonal = numpy.zeros([len(lat_edges) - 1, len(data['lev'])])
    else:
        data_zonal = numpy.zeros(len(lat_edges) - 1)

    # TODO: This will not work well for very large datasets; we need to figure
    # out a better way of doing this. May need to just remap the unstructured
    # grid first and compute a zonal mean in the usual way
    for ilat in range(len(lat_edges) - 1):
        
        # Find latitude bounds
        lat1 = lat_edges[ilat]
        lat2 = lat_edges[ilat+1]
        
        # Calculate latitude centers from bounds
        lat_centers[ilat] = (lat_edges[ilat+1] + lat_edges[ilat]) / 2.0

        # Calculate mean for this latitude band
        data_band = data.where(old_lat > lat1).where(old_lat <= lat2)
        weights_band = weights.where(old_lat > lat1).where(old_lat <= lat2)
        data_zonal[ilat, ...] = (weights_band * data_band).sum(dim=avg_dims) / weights_band.sum(dim=avg_dims)
        
    # Turn these into DataArrays
    # TODO: generalize to arbitrary dimensions and ordering
    from xarray import DataArray
    lat_centers = DataArray(
        lat_centers,
        dims=('lat',),
        attrs={'long_name': 'Latitude', 'units': 'Degrees north'}
    )
    if 'lev' in data.dims:
        dims = ('lat', 'lev')
        coords = {'lat': lat_centers, 'lev': data.lev}
    else:
        dims = ('lat')
        coords = {'lat': lat_centers}

    data_zonal = DataArray(
        data_zonal,
        dims=dims, coords=coords, attrs=data.attrs
    )
    return data_zonal, lat_centers



