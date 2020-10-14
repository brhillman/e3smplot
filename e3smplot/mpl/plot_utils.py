#!/usr/bin/env python3
import xarray as xarray
import numpy as numpy

def open_files(*inputfiles):
    return xarray.open_mfdataset(
        sorted(inputfiles), combine='by_coords', 
        drop_variables=['P3_input_dim', 'P3_output_dim'],
        chunks={'time': 1}
    )


def fix_longitudes(lon):
    lon.values = numpy.where(lon > 180, lon - 360, lon)
    return lon


def get_data(ds, variable_name):
    if variable_name in ds.variables.keys():
        data = ds[variable_name]
    elif variable_name == 'PRECT':
        precc = get_data(ds, 'PRECC')
        precl = get_data(ds, 'PRECL')
        data = precc + precl
        data.attrs = precc.attrs
        data.attrs['long_name'] = 'Total precipitation rate'
    elif variable_name == 'CLDTOT_MISR':
        clmisr = get_data(ds, 'CLD_MISR')
        tau_min = 0.3
        data = clmisr.sel(cosp_tau=slice(tau_min, None)).sum(dim=('cosp_htmisr', 'cosp_tau'), keep_attrs=True, skipna=False)
        data.attrs['long_name'] = f'MISR total cloud area'
        data.attrs['description'] = f'sum(CLD_MISR(tau > {tau_min}))'
    elif variable_name == 'CLDLOW_MISR':
        clmisr = get_data(ds, 'CLD_MISR')
        tau_min = 0.3
        zbnd = (1, 3000)  # 0 level is no-retrieval bin
        data = clmisr.sel(
            cosp_tau=slice(tau_min, None),
            cosp_htmisr=slice(zbnd[0], zbnd[1])
        ).sum(dim=('cosp_htmisr', 'cosp_tau'), keep_attrs=True, skipna=False)
        data.attrs['long_name'] = f'MISR low-level cloud area'
        data.attrs['description'] = f'sum(CLD_MISR(tau > {tau_min}, {zbnd[0]} <= z < {zbnd[1]}))'
    elif variable_name == 'CLDMED_MISR':
        clmisr = get_data(ds, 'CLD_MISR')
        tau_min = 0.3
        zbnd = (3000, 7000)
        data = clmisr.sel(
            cosp_tau=slice(tau_min, None),
            cosp_htmisr=slice(zbnd[0], zbnd[1])
        ).sum(dim=('cosp_htmisr', 'cosp_tau'), keep_attrs=True, skipna=False)
        data.attrs['long_name'] = f'MISR mid-level cloud area'
        data.attrs['description'] = f'sum(CLD_MISR(tau > {tau_min}, {zbnd[0]} <= z < {zbnd[1]}))'
    elif variable_name == 'CLDHGH_MISR':
        clmisr = get_data(ds, 'CLD_MISR')
        tau_min = 0.3
        zbnd = (7000, None)
        data = clmisr.sel(
            cosp_tau=slice(tau_min, None),
            cosp_htmisr=slice(zbnd[0], zbnd[1])
        ).sum(dim=('cosp_htmisr', 'cosp_tau'), keep_attrs=True, skipna=False)
        data.attrs['long_name'] = f'MISR high-level cloud area'
        data.attrs['description'] = f'sum(CLD_MISR(tau > {tau_min}, {zbnd[0]} <= z))'
    elif variable_name == 'CLDLOW_ISCCP':
        clisccp = get_data(ds, 'FISCCP1_COSP')
        tmin = 0.3
        pbnd = (None, 68000)
        pmin = 68000
        pmax = numpy.inf
        data = clisccp.where(clisccp.cosp_tau>tmin, drop=True).where((clisccp.cosp_prs>pmin) & (clisccp.cosp_prs<=pmax), drop=True)
        data = data.sum(dim=('cosp_prs', 'cosp_tau'), keep_attrs=True, skipna=False)
        data.attrs['long_name'] = f'ISCCP low-level cloud area'
        data.attrs['description'] = f'sum(FISCCP1_COSP(tau > {tmin}, {pmin} < prs <= {pmax}))'
    elif variable_name == 'CLDMED_ISCCP':
        clisccp = get_data(ds, 'FISCCP1_COSP')
        tmin = 0.3
        pmin = 44000
        pmax = 68000
        data = clisccp.where(clisccp.cosp_tau>tmin, drop=True).where((clisccp.cosp_prs>pmin) & (clisccp.cosp_prs<=pmax), drop=True)
        data = data.sum(dim=('cosp_prs', 'cosp_tau'), keep_attrs=True, skipna=False)
        data.attrs['long_name'] = f'ISCCP mid-level cloud area'
        data.attrs['description'] = f'sum(FISCCP1_COSP(tau > {tmin}, {pmin} < prs <= {pmax}))'
    elif variable_name == 'CLDHGH_ISCCP':
        clisccp = get_data(ds, 'FISCCP1_COSP')
        tmin = 0.3
        pmin = -numpy.inf
        pmax = 44000
        data = clisccp.where(clisccp.cosp_tau>tmin, drop=True).where((clisccp.cosp_prs>pmin) & (clisccp.cosp_prs<=pmax), drop=True)
        data = data.sum(dim=('cosp_prs', 'cosp_tau'), keep_attrs=True, skipna=False)
        data.attrs['long_name'] = f'ISCCP high-level cloud area'
        data.attrs['description'] = f'sum(FISCCP1_COSP(tau > {tmin}, {pmin} < prs <= {pmax}))'
    else:
        raise NameError('%s not found in dataset'%variable_name)

    # Adjust units
    if variable_name in ('PRECC', 'PRECL', 'PRECT'):
        if data.attrs['units'].lower() == 'm/s':
            attrs = data.attrs
            data = 60 * 60 * 24 * 1e3 * data
            data.attrs = attrs
            data.attrs['units'] = 'mm/day'
    return data


def area_average(data, weights, dims=None):

    '''Calculate area-weighted average over dims.'''

    if dims is None: dims = data.dims

    # Need to broadcast weights to make sure they have the
    # same size/shape as data. For example, data is (lat, lon)
    # but we passed weights with shape (lat,), or data is
    # (time, ncol) but we passed weights with shape (ncol,).
    weights, *__ = xarray.broadcast(weights, data)

    # Mask weights consistent with data so we do not miscount
    # missing columns
    weights = weights.where(data.notnull())

    # Do the averaging
    data_mean = (weights * data).sum(dim=dims) / weights.sum(dim=dims)

    # Copy over attributes, which we lose in the averaging
    # calculation
    data_mean.attrs = data.attrs
    data_mean.name = data.name

    # Return averaged data
    return data_mean
