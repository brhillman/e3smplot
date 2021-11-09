import numpy
import xarray
import scipy.sparse
import sys

################################################################
# This function is one of two codes contributed by Lou Wicker
# of NOAA:
#
#    _nearlyequal
#    nice_cntr_levels
#
# There was a "nice_mnmxintvl", but we decided to combine this
# into one function, "nice_cntr_levels".
#
def nearlyequal(a, b, sig_digit=None):
    """ Measures the equality (for two floats), in unit of decimal significant 
        figures.  If no sigificant digit is specified, default is 7 digits. """

    if sig_digit is None or sig_digit > 7:
        sig_digit = 7
    if a == b:
        return True
    difference = abs(a - b)
    avg = abs((a + b)/2)
    
    return numpy.log10(avg / difference) >= sig_digit
    

################################################################
# This function is one of two codes contributed by Lou Wicker
# of NOAA:
#
#    _nearlyequal
#    nice_cntr_levels
#
# There was a "nice_mnmxintvl", but we decided to combine this
# into one function, "nice_cntr_levels".
#
def nice_cntr_levels(lmin, lmax, outside=True, max_steps=15, cint=None, returnLevels=False, aboutZero=False):
    """ Description: Given min and max values of a data domain and the maximum
                     number of steps desired, determines "nice" values of 
                     for endpoints, the contour value, and if requested, an numpy
                     array containing the individual contour levels.
                     through the data domain. A flag controls whether the max 
                     and min are inside or outside the data range.  Another
                     flag can make sure that the 0 contour is in the array.
  
        In Args: float   lmin         the minimum value of the domain
                 float   lmax         the maximum value of the domain
                 int     max_steps    the maximum number of steps desired
                 logical outside      controls whether return min/max fall just
                                      outside or just inside the data domain.
                     if outside: 
                         min_out <= min < min_out + step_size
                                         max_out >= max > max_out - step_size
                     if inside:
                         min_out >= min > min_out - step_size
                                         max_out <= max < max_out + step_size
      
                 float    cint      if specified, the contour interval is set 
                                    to this, and the max/min bounds, based on 
                                    "outside" are returned.
                 returnLevels:  if True, an additional argument is returned that
                                is a numpy array containing the contour levels
                  aboutZero:    if True, makes sure that the contour interval will
                                be centered about zero.
      
      
        Out Args: min_out     a "nice" minimum value
                  max_out     a "nice" maximum value  
                  step_size   a step value such that 
                                     (where n is an integer < max_steps):
                                      min_out + n * step_size == max_out 
                                      with no remainder 
                  clevels     if returnLevels=True, a numpy array containing the contour levels
      
        If max==min, or a contour interval cannot be computed, returns "None"
     
        Algorithm mimics the NCAR NCL lib "nice_cntr_levels"; code adapted from 
        "nicevals.c" however, added the optional "cint" arg to facilitate user 
        specified specific interval.
     
        Lou Wicker, NSSL, August 2010 """

    table = numpy.array([1.0,2.0,2.5,4.0,5.0,10.0,20.0,25.0,40.0,50.0,100.0,200.0,
                      250.0,400.0,500.0])

    if nearlyequal(lmax,lmin):
        return None
    
    # Help people like me who can never remember - flip max/min if inputted reversed
    if lmax < lmin:
        amax = lmin
        amin = lmax
    else:
        amax = lmax
        amin = lmin

# If aboutZero == True, adjust the max/mins so that they symmetrically straddle zero

    if aboutZero:
        vmax = max(abs(lmax), abs(lmin))
        amax =  vmax
        amin = -vmax

    d = 10.0**(numpy.floor(numpy.log10(amax - amin)) - 2.0)

    # List of candidate step sizes
    if cint is None or cint == 0.0:
        t = table * d
    else:
        t = cint

    # Get candidate min/maxes
    if outside:
        am1 = numpy.floor(amin/t) * t
        ax1 = numpy.ceil(amax/t)  * t
    else:
        am1 = numpy.ceil(amin/t) * t
        ax1 = numpy.floor(amax/t)  * t
 
    if cint is None or cint == 0.0:

        # Candidate nsteps
        cints = (ax1 - am1) / t

        # DEBUG LINE BELOW
        #print(t, am1, ax1, cints)

        try:
            index = numpy.where(cints < max_steps)[0][0]
        except IndexError:
            return None

        if returnLevels:
            return am1[index], ax1[index], t[index], numpy.arange(am1[index], ax1[index]+t[index], t[index])
        else:
            return am1[index], ax1[index], t[index]

    else:

        if returnLevels:
            return am1, ax1, cint, numpy.arange(am1, ax1+cint, cint)
        else:
            return am1, ax1, cint
        

def myprint(*args, **kwargs):
    print(*args, **kwargs); sys.stdout.flush()


def apply_map(da, map_file, template=None, verbose=False):

    # Allow for passing either a mapping file name or a xarray.Dataset
    if isinstance(map_file, xarray.Dataset):
        ds_map = map_file
    else:
        if verbose: myprint('Open map file as xarray.Dataset object')
        ds_map = xarray.open_mfdataset(map_file)

    # Do the remapping
    if verbose: myprint('Create weights as a COO matrix')
    weights = scipy.sparse.coo_matrix((ds_map['S'].values, (ds_map['row'].values-1, ds_map['col'].values-1)))
    if verbose: myprint('Flatten data array')
    if len(da.shape) == 1:
        da_flat = da.data
    elif len(da.shape) == 2:
        da_flat = da.data.reshape([da.shape[0]*da.shape[1]])
    if verbose: myprint('Apply weights via matrix multiply')
    da_regrid = weights.dot(da_flat)

    # Figure out coordinate variables and whether or not we should reshape the
    # output before returning
    if verbose: myprint('Reshape output')
    if isinstance(template, xarray.DataArray):
        da_regrid = xarray.DataArray(
            da_regrid.reshape(template.shape), dims=template.dims, coords=template.coords,
            attrs=da.attrs
        )
        x = da_regrid.lon
        y = da_regrid.lat
    elif len(ds_map.dst_grid_dims) == 2:
        # Get lon and lat coordinates from mapping file
        x = ds_map.xc_b.values.reshape(ds_map.dst_grid_dims.values[::-1])[0,:]
        y = ds_map.yc_b.values.reshape(ds_map.dst_grid_dims.values[::-1])[:,0]

        # Reshape to expected output
        da_regrid = xarray.DataArray(
            da_regrid.reshape(ds_map.dst_grid_dims.values[::-1]),
            dims=('lat', 'lon'),
            coords={'lat': y, 'lon': x},
            attrs=da.attrs,
        )
    else:
        x = ds_map.xc_b
        y = ds_map.yc_b

    # Return remapped array and coordinate variables
    if verbose: myprint('Return remapped data, x, y')
    return da_regrid, x, y


# update_progress() : Displays or updates a console progress bar
def update_progress(iteration, num_iterations, bar_length=10):
    '''
    Display and update a console progress bar
    '''

    # Get progress as a fraction and compute size of "block" filled to visually
    # represent fraction completed.
    progress = 1.0 * iteration / num_iterations
    block = int(round(bar_length * progress))

    # Get status; if status < 1 there will be no newline, so we need to add one
    # explicitly on the last iteration.
    if progress >= 1:
        status = "\r\n"
    else:
        status = ""

    # Display appropriate text to build status bar
    text = "\rProgress: [%s] %i of %i (%.2f%%)%s"%(
        "#"*block + "-"*(bar_length-block),
        iteration, num_iterations,
        progress*100, status
    )
    sys.stdout.write(text)
    sys.stdout.flush()


def interp_along_axis(xout, xin, yin, axis=0, *args, **kwargs):
    """
    Perform 1D interpolation along 1D slices of a ND array.
    """

    # Make sure input arrays are sized appropriately
    assert xin.shape == yin.shape

    # Setup output array
    shape_in  = yin.shape
    shape_out = list(shape_in)
    shape_out[axis] = len(xout)
    yout = numpy.empty(shape_out)

    # Apply function along axis specified by looping over
    # indices across the other axes and applying function
    # to 1d slices; equivalent to, e.g.,
    # for ii in yin.shape[0]:
    #     for jj in yin.shape[1]:
    #        for kk in yin.shape[2]:
    #            yout[ii,jj,kk,:] = numpy.interp(xout, xin[ii,jj,kk,:], yin[ii,jj,kk,:])
    Ni, Nk = yin.shape[:axis], yin.shape[axis+1:]
    for ii in numpy.ndindex(Ni):
        for kk in numpy.ndindex(Nk):
            yout[ii + numpy.s_[...,] + kk] = numpy.interp(
                xout, xin[ii + numpy.s_[:,] + kk], yin[ii + numpy.s_[:,] + kk], *args, **kwargs
            )
    return yout
