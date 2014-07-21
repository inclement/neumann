'''
Cython functions for finding Neumann domains.
'''

import numpy as n
cimport numpy as n

cimport cython
from libc.math cimport sin

cpdef double random_wave_function(double [:, :] wvs, double [:] amps,
                                  double [:] phases, double x, double y):
    '''Calculates the random wave function at the given x, y position with
    the given coefficients.'''

    cdef double result = 0.0
    cdef int i

    for i in range(0, len(wvs)):
        result += (amps[i] *
                   sin(dotprod(wvs[i, 0],
                               wvs[i, 1],
                               x, y) +
                       phases[i])
                       )
    return result


cdef double dotprod(double wvx, double wvy, double x, double y):
    '''Fast dot product meeting our exact needs.'''
    return wvx * x + wvy * y
