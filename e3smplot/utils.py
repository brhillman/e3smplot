import numpy

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
        
