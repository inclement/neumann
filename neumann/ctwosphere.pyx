from __future__ import print_function
import numpy as n
cimport numpy as n

from libc.math cimport atan2, abs, pow, acos, asin, tan, sin, cos, sqrt
cdef double pi = 3.14159265358979323846

cpdef complex spherical_harmonic_2d(long l, long m, double theta, double phi):
    absm = abs(m)
    cdef double lpart = assoc_legendre(l, absm, cos(theta))
    cdef complex phipart = cos(absm*phi) + 1j*sin(absm*phi)

    cdef complex result = lpart * phipart
    if m < 0:
        result = (-1.)**absm * n.conj(result)
    return result

cpdef assoc_legendre(long l, long m, double x):
    '''Computes the *renormalised* (including prefactor) associated
    legendre polynomial.'''
    cdef long i, ll
    cdef double fact, oldfact, pll, pmm, pmmp1, omx2

    pmm = 1.0

    if m > 0:
        omx2 = (1.-x) * (1.+x)
        fact = 1.0
        for i in range(1, m+1):
            pmm *= omx2*fact / (fact+1.)
            fact += 2.

    pmm = sqrt((2*m+1) * pmm / (4*pi))

    if m & 1:  # ? <- some bit operation?
        pmm = -pmm
    if l == m:
        return pmm
    else:
        pmmp1 = x*sqrt(2.0*m + 3.0) * pmm
        if (l == (m+1)):
            return pmmp1
        else:
            oldfact = sqrt(2.0*m + 3.0)
            for ll in range(m+2, l+1):
                fact = sqrt((4.0*ll*ll - 1.0) / (ll*ll - m*m))
                pll = (x*pmmp1 - pmm/oldfact)*fact
                oldfact = fact
                pmm = pmmp1
                pmmp1 = pll
            return pll

cpdef assoc_legendre_norm(long l, long m, double x):
    'Computes the non-renormalisd associated legendre polynomial.'''
    cdef double prod = 1.0
    cdef long j
    for j in range(l-m+1, l+m+1):
        prod *= j
    return sqrt(4.0*pi*prod / (2*l + 1)) * assoc_legendre(l, m, x)
