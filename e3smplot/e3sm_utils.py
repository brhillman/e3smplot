#!/usr/bin/env python3

import numpy
import xarray, xarray.ufuncs
import os, os.path
import subprocess
import sys
import dask
import cftime
import warnings

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


def convert_time(ds):
    # Convert cftime coordinate
    if isinstance(ds.time.values[0], cftime._cftime.DatetimeNoLeap) or isinstance(ds.time.values[0], cftime._cftime.DatetimeJulian):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ds['time'] = ds.indexes['time'].to_datetimeindex()
        except:
            print('Could not convert times to datetimeindex. But proceeding anyway...')

    return ds


def open_dataset(*files, time_offset=None, **kwargs):
    # Open dataset as a dask array
    with dask.config.set(**{'array.slicing.split_large_chunks': True}):
        ds = xarray.open_mfdataset(
            sorted(files), combine='by_coords', 
            drop_variables=('P3_output_dim', 'P3_input_dim'), **kwargs
        )

    if time_offset is not None:
        print('Adding year offset...')
        ds['time'] = ds['time'] + time_offset

    # Rename coordinate variables
    if 'latitude' in ds.dims: ds = ds.rename({'latitude': 'lat'})
    if 'longitude' in ds.dims: ds = ds.rename({'longitude': 'lon'})

    if 'msshf' in ds.variables.keys():
        ds['msshf'].data = -ds['msshf'].data

    # Fix times so we can subset later
    ds = convert_time(ds)

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
    elif field == 'FSNT':
        try:
            data = get_data(dataset, 'FSNTOA')
            warnings.warn('Using FSNTOA instead of FSNT')
        except:
            pass
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
    elif field == 'CLDHGH':
        if 'cldarea_high_3h' in dataset.variables.keys():
            data = dataset['cldarea_high_3h']
            data.attrs['long_name'] = 'High-topped cloud area fraction'
            data.attrs['units'] = '%'
    elif field == 'CLDLIQICE':
        cldliq = get_data(dataset, 'CLDLIQ')
        cldice = get_data(dataset, 'CLDICE')
        data = cldliq + cldice
        data.attrs = cldliq.attrs
        data.attrs['long_name'] = 'Combined liq+ice condensate'
    elif field == 'PRECT':
        if 'PRECC' in dataset and 'PRECL' in dataset:
            precc = get_data(dataset, 'PRECC')
            precl = get_data(dataset, 'PRECL')
            data = precc + precl
            data.attrs = precc.attrs
            data.attrs['long_name'] = 'Total precipitation rate'
        elif 'mtpr' in dataset:
            data = get_data(dataset, 'mtpr')
            # Convert units from kg m^-2 s^-1 to mm/day
            density = 1.0e3 # kg / m^3
            data.data = data.data / density * 1e3 * 60 * 60 * 24
            data.name = 'PRECT'
            data.attrs['long_name'] = 'Total precipitation rate'
            data.attrs['units'] = 'mm/day'
        elif 'precipitationCal' in dataset:
            data = get_data(dataset, 'precipitationCal').copy()
            # Convert from mm/hr to mm/day
            data.data = data.data * 24.0
            data.name = 'PRECT'
            data.attrs['long_name'] = 'Total precipitation rate'
            data.attrs['units'] = 'mm/day'
    elif field == 'TGCLDWP':
        clwp = get_data(dataset, 'TGCLDLWP')
        ciwp = get_data(dataset, 'TGCLDIWP')
        data = clwp + ciwp
        data.attrs = clwp.attrs
        data.attrs['long_name'] = 'Total gridbox cloud water path'
    elif field == 'LIQ_CLD_MASK':
        data = get_liq_cld_mask(dataset)
    elif field == 'ICE_CLD_MASK':
        data = get_ice_cld_mask(dataset)
    elif field == 'TOT_CLD_MASK':
        data = get_tot_cld_mask(dataset)
    elif field == 'LIQ_CLD_AREA':
        data = get_liq_cld_area(dataset)
    elif field == 'ICE_CLD_AREA':
        data = get_ice_cld_area(dataset)
    elif field == 'TOT_CLD_AREA':
        data = get_tot_cld_area(dataset)

    # Derived MISR cloud types
    elif field == 'CLDTOT_MISR':
        tau1, tau2 = 0.3, numpy.inf
        cth1, cth2 = -numpy.inf, numpy.inf
        data = get_data(dataset, 'CLD_MISR').sel(cosp_tau=slice(tau1, tau2), cosp_htmisr=slice(cth1, cth2)).sum(dim=('cosp_tau', 'cosp_htmisr'), keep_attrs=True)
        data.attrs['long_name'] = 'MISR-simulated total cloud area'
    elif field == 'CLDLOW_MISR':
        tau1, tau2 = 0.3, numpy.inf
        cth1, cth2 = 0, 3000
        data = get_data(dataset, 'CLD_MISR').sel(cosp_tau=slice(tau1, tau2), cosp_htmisr=slice(cth1, cth2)).sum(dim=('cosp_tau', 'cosp_htmisr'), keep_attrs=True)
        data.attrs['long_name'] = 'MISR-simulated low-topped (< 3 km) cloud area'
    elif field == 'CLDMED_MISR':
        tau1, tau2 = 0.3, numpy.inf
        cth1, cth2 = 3000, 7000
        data = get_data(dataset, 'CLD_MISR').sel(cosp_tau=slice(tau1, tau2), cosp_htmisr=slice(cth1, cth2)).sum(dim=('cosp_tau', 'cosp_htmisr'), keep_attrs=True)
        data.attrs['long_name'] = 'MISR-simulated mid-topped (3 - 7 km) cloud area'
    elif field == 'CLDHGH_MISR':
        tau1, tau2 = 0.3, numpy.inf
        cth1, cth2 = 7000, numpy.inf
        data = get_data(dataset, 'CLD_MISR').sel(cosp_tau=slice(tau1, tau2), cosp_htmisr=slice(cth1, cth2)).sum(dim=('cosp_tau', 'cosp_htmisr'), keep_attrs=True)
        data.attrs['long_name'] = 'MISR-simulated high-topped (> 7 km) cloud area'
    
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
    elif field == 'T':
        if 't' in dataset.variables.keys():
            data = get_data(dataset, 't')
    elif field == 'Q':
        if 'q' in dataset.variables.keys():
            data = get_data(dataset, 'q')
    elif field == 'PS':
        if 'sp' in dataset.variables.keys():
            data = get_data(dataset, 'sp')
    elif field == 'WINDSPD_10M':
        if 'u10' in dataset.variables.keys() and 'v10' in dataset.variables.keys():
            zonalwind10m = get_data(dataset,'u10')
            meridwind10m = get_data(dataset,'v10')
            data = (zonalwind10m**2 + meridwind10m**2)**0.5
            data.attrs['long_name'] = '10m wind speed'
            data.attrs['units'] = 'm/s'
    elif field == 'LHFLX':
        if 'mer' in dataset.variables.keys():
            flux_down = dataset['mer']
            data = flux_down*-1.*2.501e6 #Scaled evap with latent heat of vaporization in shr_const_mod.F90
            data.attrs['long_name'] = 'Latent heat flux'
            data.attrs['units'] = 'W/m2'
    elif field == 'rad_heating_pdel':
        if 'QRS' in dataset.variables.keys() and 'QRL' in dataset.variables.keys():
            qrs = dataset['QRS']
            qrl = dataset['QRL']
            data = (qrs + qrl) * pdel
        else:
            raise RuntimeError('Cannot get rad_heating_pdel')
    elif field == 'LW_flux_up':
        data = get_data(dataset, 'FUL')
    elif field == 'LW_flux_dn':
        data = get_data(dataset, 'FDL')
    elif field == 'SW_flux_up':
        data = get_data(dataset, 'FUS')
    elif field == 'SW_flux_dn':
        data = get_data(dataset, 'FDS')
    elif field == 'LW_clrsky_flux_up':
        data = get_data(dataset, 'FULC')
    elif field == 'LW_clrsky_flux_dn':
        data = get_data(dataset, 'FDLC')
    elif field == 'SW_clrsky_flux_up':
        data = get_data(dataset, 'FUSC')
    elif field == 'SW_clrsky_flux_dn':
        data = get_data(dataset, 'FDSC')
    elif field == 'in_cloud_ice_path':
        data = 1e3 * get_data(dataset, 'ICLDIWP')
    elif field == 'in_cloud_water_path':
        d_tot = 1e3 * get_data(dataset, 'ICLDTWP')
        d_ice = get_data(dataset, 'ICLDIWP')
        data = (d_tot - d_ice)
        data.attrs = d_tot.attrs
        data['long_name'] = 'In-cloud water path'
    elif field == 'eff_radius_qc':
        data = get_data(dataset, 'REL')
    elif field == 'eff_radius_qi':
        data = get_data(dataset, 'REI')
    elif field == 'cldfrac_tot':
        data = get_data(dataset, 'CLOUD')
    elif field == 'qc':
        data = get_data(dataset, 'CLDLIQ')
    elif field == 'qi':
        data = get_data(dataset, 'CLDICE')
    elif field == 'qv':
        data = get_data(dataset, 'Q')
    elif field == 'SW_cloud_effect_up':
        dall = get_data(dataset, 'SW_flux_up')
        dclr = get_data(dataset, 'SW_clrsky_flux_up')
        data = dall - dclr
        data.attrs = dall.attrs
        data.attrs['long_name'] = 'SW cloud effect on upwelling flux'
    elif field == 'SW_cloud_effect_dn':
        dall = get_data(dataset, 'SW_flux_dn')
        dclr = get_data(dataset, 'SW_clrsky_flux_dn')
        data = dall - dclr
        data.attrs = dall.attrs
        data.attrs['long_name'] = 'SW cloud effect on downwelling flux'
    elif field == 'LW_cloud_effect_up':
        dall = get_data(dataset, 'LW_flux_up')
        dclr = get_data(dataset, 'LW_clrsky_flux_up')
        data = dall - dclr
        data.attrs = dall.attrs
        data.attrs['long_name'] = 'LW cloud effect on upwelling flux'
    elif field == 'LW_cloud_effect_dn':
        dall = get_data(dataset, 'LW_flux_dn')
        dclr = get_data(dataset, 'LW_clrsky_flux_dn')
        data = dall - dclr
        data.attrs = dall.attrs
        data.attrs['long_name'] = 'LW cloud effect on downwelling flux'
    elif field == 'cosine_solar_zenith_angle':
        data = get_data(dataset, 'COSZRS')

    # Automatically grab top or bottom by appending _toa or _sfc to field name
    elif field[-4:] == '_sfc':
        data = get_data(dataset, field[:-4])
        if 'ilev' in data.dims:
            data = data.isel(ilev=-1)
        elif 'lev' in data.dims:
            data = data.isel(lev=-1)
        else:
            raise RuntimeError(f'Not sure what to do with dims {data.dims}')
        data.attrs['long_name'] = data.attrs['long_name'] + ' at surface'
    elif field[-4:] == '_toa' or field[-4:] == '_tom':
        data = get_data(dataset, field[:-4])
        if 'ilev' in data.dims:
            data = data.isel(ilev=0)
        elif 'lev' in data.dims:
            data = data.isel(lev=0)
        else:
            raise RuntimeError(f'Not sure what to do with dims {data.dims}')
        data.attrs['long_name'] = data.attrs['long_name'] + ' at TOM'
    
    # Check if we were able to find or derive the requested field
    if data is None:
        raise ValueError(f'{field} not found in dataset and no variables found to derive') 

    # Adjust long_name
    if field == 'TREFHT':
        data.attrs['long_name'] = 'Reference height temperature'
    elif field == 'SHFLX':
        data.attrs['long_name'] = 'Sensible heat flux'
    elif field == 'LHFLX':
        data.attrs['long_name'] = 'Latent heat flux'
    elif field == 'TMQ':
        data.attrs['long_name'] = 'Total precipitable water'
    elif field == 'PS':
        data.attrs['long_name'] = 'Surface pressure'
    elif field == 'WINDSPD_10M':
        data.attrs['long_name'] = '10m wind speed'    
    

    # Adjust units if necessary
    if field in ('TMCLDLIQ', 'TMCLDICE', 'TGCLDLWP', 'TGCLDIWP', 'TGCLDWP'):
        if data.units == 'kg/m2':
            attrs = data.attrs
            data = 1e3 * data
            data.attrs = attrs
            data.attrs['units'] = 'g/m2'
    elif field in ('CLDTOT', 'CLDLOW', 'CLDMED', 'CLDHGH',
                   'cltcalipso', 'cltcalipso_liq', 'cltcalipso_ice',
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
    elif field in ('PS'):
        if data.attrs['units'].lower() == 'Pa':
            attrs = data.attrs
            data = 1e-2 * data
            data.attrs = attrs
            data.attrs['units'] = 'hPa'
        
    # Adjust long_name if necessary
    if field == 'PRECT':
        data.attrs['long_name'] = 'Total precipitation rate'

    return data

# Compute a cloud mask based on a threshold of liquid and ice water content
def get_liq_cld_mask(ds, threshold=1e-5):
    cldliq = get_data(ds, 'CLDLIQ')
    cld_mask = (cldliq > threshold).astype(numpy.float32)
    cld_mask.attrs = {
        'name': 'LIQ_CLD_MASK',
        'long_name': 'Liquid cloud mask',
        'units': 'none',
        'description': f'CLDLIQ > {threshold}',
    }
    return cld_mask

def get_ice_cld_mask(ds, threshold=1e-5):
    cldice = get_data(ds, 'CLDICE')
    cld_mask = (cldice > threshold).astype(numpy.float32)
    cld_mask.attrs = {
        'name': 'ICE_CLD_MASK',
        'long_name': 'Ice cloud mask',
        'units': 'none',
        'description': f'CLDICE > {threshold}',
    }
    return cld_mask

def get_tot_cld_mask(ds):
    liq_mask = get_liq_cld_mask(ds)
    ice_mask = get_ice_cld_mask(ds)
    cld_mask = ((liq_mask + ice_mask) > 0).astype(numpy.float32) #((liq_mask > 0) | (ice_mask > 0))
    cld_mask.attrs = {
        'name': 'TOT_CLD_MASK',
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
    cld_area = (cld_mask > 0).any(dim='lev').astype(numpy.float32)
    cld_area.attrs = {
        'name': 'LIQ_CLD_AREA',
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
    cld_area = (cld_mask > 0).any(dim='lev').astype(numpy.float32)
    cld_area.attrs = {
        'name': 'ICE_CLD_AREA',
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
    cld_area = (cld_mask > 0).any(dim='lev').astype(numpy.float32)
    cld_area.attrs = {
        'name': 'TOT_CLD_AREA',
        'long_name': 'Cloud area mask',
        'units': 'none',
        'description': 'any(cld_mask > 0)',
    }
    return cld_area


def can_retrieve_field(f, v):
    '''
    Check if named field can be retrieved/derived from file.

    Inputs:
        f: filename
        v: variable name
    Returns:
        True if variable can be retrieved/derived, False otherwise.
    '''
    result = False
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            with open_dataset(f) as ds:
                try: 
                    d = get_data(ds, v)
                    result = True
                except:
                    result = False
        except:
            result = False
    return result


def get_common_time_range(*datasets):
    t1 = max(ds.time[0].values for ds in datasets)
    t2 = min(ds.time[-1].values for ds in datasets)
    return t1,t2


def get_grid_name(ds, **kwargs):

    # If ds is a string, then we need to open the dataset still
    if isinstance(ds, str):
        ds = xarray.open_mfdataset(ds)

    # This is maybe...a little hacky and fragile. Determine grid based on number
    # of columns in the dataset. For well-annotated netcdf files, we could grab
    # this info from global attributes I think. For example, for EAM output, we
    # have a "ne" and a "fv_nphys" attribute in the global metadata we could
    # grab. But this works for now.
    grid_ncol_dict = {
         384:      'ne4pg2',
         866:      'ne4np4', 
         21600:    'ne30pg2',
         345600:   'ne120pg2',
         1572864:  'ne256pg2',
         25165824: 'ne1024pg2',
    }

    # Look for dims
    grid = None
    if 'ncol' in ds.dims: 
        grid = grid_ncol_dict[ds.sizes['ncol']]
    elif 'grid_dims' in ds:
        print(ds.grid_dims)
        if len(ds.grid_dims) == 2:
            grid = f'{ds.grid_dims[1]}x{ds.grid_dims[0]}'
    else:
        for (x, y) in zip(('lon', 'longitude'), ('lat', 'latitude')):
            if x in ds.dims and y in ds.dims:
                grid = '{}x{}'.format(ds.sizes[y], ds.sizes[x])

    # Check to make sure we found a grid
    if grid is None: raise ValueError('No valid grid found.')

    return grid


def check_ret(p):
    if p.returncode != 0:
        raise RuntimeError('Subprocess command failed: {}'.format(' '.join([str(s) for s in p.args])))


def get_scrip_grid_ds(ds, grid_root):
    grid_name = get_grid_name(ds)
    grid_file = f'{grid_root}/{grid_name}_scrip.nc'
    if not os.path.exists(grid_file):
        print(f'Grid {grid_file} not found, trying to create...', end=''); sys.stdout.flush()
        os.makedirs(grid_root, exist_ok=True)
        # For pg2 grids, we can create the grid file manually
        phys_grids = ('ne4pg2', 'ne30pg2', 'ne120pg2', 'ne256pg2', 'ne1024pg2')
        if grid_name in phys_grids:
            # infer ne and pg
            grid_ne, grid_pg = grid_name.split('ne')[1].split('pg')
            # Generate exodus element file
            check_ret(subprocess.run(f'GenerateCSMesh --alt --res {grid_ne} --file {grid_root}/ne{grid_ne}.g'.split(' ')))
            # Generate pg2 grid file
            check_ret(subprocess.run(f'GenerateVolumetricMesh --in {grid_root}/ne{grid_ne}.g --out {grid_root}/ne{grid_ne}pg{grid_pg}.g --np {grid_pg} --uniform'.split(' ')))
            # Convert exodus to scrip
            check_ret(subprocess.run(f'ConvertExodusToSCRIP --in {grid_root}/ne{grid_ne}pg{grid_pg}.g --out {grid_file}'.split(' ')))
        else:
            # Infer grid from coordinate data; use ncks to do this, so first we
            # need to write out the coordinate data to a temporary file
            tmp_file = None
            for (x, y) in zip(('lon', 'longitude', 'x'), ('lat', 'latitude', 'y')):
                if x in ds and y in ds:
                    tmp_file = f'{grid_root}/{ds.size[y]}x{ds.size[x]}_latlon.nc'
                    xarray.Dataset({x: ds[x], y: ds[y]}).to_netcdf(tmp_file)
            if tmp_file is not None:
                try:
                    os.remove(f'{grid_root}/foo.nc')
                except:
                    pass
                check_ret(subprocess.run(f'ncks --rgr infer --rgr scrip={grid_file} {tmp_file} {grid_root}/foo.nc'.split(' ')))
            else:
                raise RuntimeError('Failed to generate a coordinate file.')
        print('done.'); sys.stdout.flush()

    return grid_file


def get_scrip_grid(data_file, grid_root):

    with xarray.open_dataset(data_file) as ds: grid_name = get_grid_name(ds)
    grid_file = f'{grid_root}/{grid_name}_scrip.nc'
    if not os.path.exists(grid_file):
        print(f'Grid {grid_file} not found, trying to create...', end=''); sys.stdout.flush()
        os.makedirs(grid_root, exist_ok=True)
        # For pg2 grids, we can create the grid file manually
        phys_grids = ('ne4pg2', 'ne30pg2', 'ne120pg2', 'ne256pg2', 'ne1024pg2')
        if grid_name in phys_grids:
            # infer ne and pg
            grid_ne, grid_pg = grid_name.split('ne')[1].split('pg')
            # Generate exodus element file
            check_ret(subprocess.run(f'GenerateCSMesh --alt --res {grid_ne} --file {grid_root}/ne{grid_ne}.g'.split(' ')))
            # Generate pg2 grid file
            check_ret(subprocess.run(f'GenerateVolumetricMesh --in {grid_root}/ne{grid_ne}.g --out {grid_root}/ne{grid_ne}pg{grid_pg}.g --np {grid_pg} --uniform'.split(' ')))
            # Convert exodus to scrip
            check_ret(subprocess.run(f'ConvertExodusToSCRIP --in {grid_root}/ne{grid_ne}pg{grid_pg}.g --out {grid_file}'.split(' ')))
        else:
            # Infer grid from coordinate data; make sure foo.nc does not exist or script will require user interaction
            try:
                os.remove(f'{grid_root}/foo.nc')
            except:
                pass
            check_ret(subprocess.run(f'ncks --rgr infer --rgr scrip={grid_file} {data_file} {grid_root}/foo.nc'.split(' ')))
        print('done.'); sys.stdout.flush()

    return grid_file


def create_scrip_grid(x, y, grid_root):
    '''
    Create a new SCRIP-formatted grid file given lon/lat coordinates. This is to
    be used for mapping unstructured data to a lon/lat grid for easier analysis.
    '''

    # Construct name for scrip grid
    nx = len(x)
    ny = len(y)
    grid_file = f'{grid_root}/{ny}x{nx}_scrip.nc'

    # Create scrip grid using nco
    command = 'ncremap -G ttl=\'Equi-Angular_grid_{ny}x{nx}\'#latlon={ny},{nx}#lat_typ=uni#lon_typ=grn_ctr -g {grid_file}'.split(' ')
    print(command)
    check_ret(subprocess.run(f'ncremap -G ttl=\'Equi-Angular_grid_{ny}x{nx}\'#latlon={ny},{nx}#lat_typ=uni#lon_typ=grn_ctr -g {grid_file}'.split(' ')))

    return grid_file


def get_mapping_file(source_file, target_file, mapping_root, method='nco'):

    # Open files and get grids
    with xarray.open_dataset(source_file) as ds: source_grid = get_grid_name(ds)
    with xarray.open_dataset(target_file) as ds: target_grid = get_grid_name(ds)

    # Check if we have a mapping file from cntl to test grid. If not, we will create one.
    map_file = f'{mapping_root}/map_{source_grid}_to_{target_grid}_{method}.nc'
    if not os.path.exists(map_file):
        print(f'Creating mapping file {map_file}...'); sys.stdout.flush()
        source_grid_file = get_scrip_grid(source_file, mapping_root)
        target_grid_file = get_scrip_grid(target_file, mapping_root)
        check_ret(subprocess.run(f'ncremap --fl_fmt=64bit_data --alg_typ={method} --src_grd={source_grid_file} --dst_grd={target_grid_file} -m {map_file}'.split(' ')))
        print('done.'); sys.stdout.flush()

    return map_file


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


def infer_mapping_root():
    from os import uname
    nodename = os.uname().nodename
    if 'cori' in nodename:
        mapping_root = '/global/cfs/cdirs/e3sm/mapping'
    else:
        raise RuntimeError(f'No mapping root defined for {nodename}')
    return mapping_root


def infer_grid_file(ds, mapping_root=None):
    # Get mapping root where we can find grid files
    if mapping_root is None: mapping_root = infer_mapping_root()
    # Infer grid file name from ncol dim
    ncol = ds.dims['ncol']
    if   ncol == 21600:  # ne30pg2
        grid_file = f'{mapping_root}/grids/ne30pg2_scrip_20200209.nc'
    elif ncol == 48602:  # ne30np4
        grid_file = f'{mapping_root}/grids/ne30np4_pentagons.091226.nc'
    else:
        raise RuntimeError(f"Grid with ncol = {ncol} unknown.")
    return grid_file


def infer_grid_coords(ds, grid_file=None, mapping_root=None):
    if grid_file is None: grid_file = infer_grid_file(ds, mapping_root=mapping_root)
    # Get coordinate variables from file
    ds_grid = xarray.open_dataset(grid_file)
    xc = ds_grid['grid_center_lon']
    yc = ds_grid['grid_center_lat']
    xv = ds_grid['grid_corner_lon']
    yv = ds_grid['grid_corner_lat']
    return xc, yc, xv, yv


def get_area_weights(ds):
    # Get weights; either use pre-computed or cosine(latitude) weights
    if 'area' in ds.variables.keys():
        wgt = ds['area']
    elif 'grid_area' in ds.variables.keys():
        wgt = ds['grid_area']
    else:
        # Use xarray.ufuncs to work on lazily evaluated dask arrays
        print('WARNING: no area variable found in dataset, falling back to cosine(latitude) weights!')
        y = get_data(ds, 'latitude')
        wgt = xarray.ufuncs.cos(y * numpy.pi / 180.0)
    return wgt


def area_average(data, weights, dims=None, **kwargs):
    
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
    data_mean = (weights * data).sum(dim=dims, **kwargs) / weights.sum(dim=dims, **kwargs)
    
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


def is_latlon(f):
    if isinstance(f, xarray.Dataset):
        ds = f
    else:
        ds = xarray.open_mfdataset(f)
    if any(dim in ds.dims for dim in ('lat', 'latitude')) and any(dim in ds.dims for dim in ('lon', 'longitude')):
        return True
    else:
        return False


def mask_all_zero(d, dims=None):
    '''
    Mask data where ALL samples are zero. This is meant to handle the case
    where missing data was accidentally written as zeros instead of fillvalues.
    '''
    if dims is None: dims = list(set(d.dims) - set(['time',]))
    return d.where(abs(d).sum(dim=dims, keep_attrs=True))
