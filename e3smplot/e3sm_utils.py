#!/usr/bin/env python3

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


# Define function to read specialized data from E3SM files
def get_data(dataset, field):
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
        
    else:
        raise NameError('%s not found in dataset.'%field)
    
    # Adjust units if necessary
    if field in ('TGCLDLWP', 'TGCLDIWP', 'TGCLDWP'):
        if data.units == 'kg/m2':
            attrs = data.attrs
            data = 1e3 * data
            data.attrs = attrs
            data.attrs['units'] = 'g/m2'
    elif field in ('cltcalipso', 'cltcalipso_liq', 'cltcalipso_ice',
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
        
    return data


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