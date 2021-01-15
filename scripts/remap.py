#!/usr/bin/env python

import xarray
import numpy
import scipy.sparse
import dask.array
import sys
from e3smplot.e3sm_utils import get_data
from e3smplot.utils import update_progress

def get_weights(ds_map):

    return scipy.sparse.coo_matrix((
        ds_map['S'].values, (ds_map['row'].values-1, ds_map['col'].values-1)
    ))


def dot_prod(weights, data_flat):
    return weights.dot(data_flat)

def apply_map(ds_map, data, ndim=1):

    # Input dimensions
    non_horiz_dims = data.shape[0:-ndim]
    horiz_dims = data.shape[-ndim:]

    weights = get_weights(ds_map)

    # Output dimensions
    x = ds_map.xc_b.values.reshape(ds_map.dst_grid_dims.values[::-1])[0,:]
    y = ds_map.yc_b.values.reshape(ds_map.dst_grid_dims.values[::-1])[:,0]
    horiz_dims_out = ds_map.dst_grid_dims.values[::-1]

    # ndim is number of horizontal dimensions to regrid over (1-2)
    horiz_flat = numpy.prod([i for i in data.shape[-ndim:]])

    print('non_horiz_dims: ', non_horiz_dims)
    print('horiz_dims    : ', horiz_dims)
    print('horiz_dims_out: ', horiz_dims_out)
    print('horiz_flat    : ', horiz_flat)

    # Do the remap
    data_flat = data.data.reshape(
        [numpy.prod(non_horiz_dims), numpy.prod(horiz_dims)]
    ).rechunk([1, numpy.prod(horiz_dims)])
    dout_flat = dask.array.zeros(
        [numpy.prod(non_horiz_dims), numpy.prod(horiz_dims_out)], 
        chunks=(1, numpy.prod(horiz_dims_out))
    )
    print('Do remap...'); sys.stdout.flush()
    dout_flat = dask.array.map_blocks(
        dot_prod, weights, data_flat, dtype=numpy.float32,
        chunks=(1, numpy.prod(horiz_dims_out)),
        meta=dout_flat
    )
    print('dout_flat.shape: ', dout_flat.shape, '; dout_flat.chunks: ', dout_flat.chunks)

    # Unflatten output array
    print('Reshape...'); sys.stdout.flush()
    dout = dout_flat.reshape(
        [*non_horiz_dims, *horiz_dims_out]
    ).rechunk([*[1 for i in range(len(non_horiz_dims))], *horiz_dims_out])
    print('dout.shape: ', dout.shape, '; dout.chunks: ', dout.chunks)
    #print('dout_flat.max(): ', dout_flat.max().compute())
    #print('dout.max(): ', dout.max().compute())

    da_out = xarray.DataArray(dout, dims=('time', 'lev', 'lat', 'lon'))

    return da_out#dout


def main(vname, mapfile, inputfile, outputfile, **kwargs):

    # Read mapping file
    ds_map = xarray.open_mfdataset(mapfile)

    # Read data
    ds = xarray.open_mfdataset(
        inputfile, drop_variables=('P3_input_dim', 'P3_output_dim'),
        chunks={'time': 1, 'lev': 1}, engine='netcdf4',
    )
    da = get_data(ds, vname)

    # Remap
    da_out = apply_map(ds_map, da, ndim=1)
    print('Write to file...'); sys.stdout.flush()
    da_out.to_netcdf(outputfile)


if __name__ == '__main__':
    import plac; plac.call(main)
