'''
Cython functions for finding Neumann domains.
'''

import numpy as n
cimport numpy as n

cimport cython
from libc.math cimport sin, atan2, sqrt, pow, round

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


cdef inline double dotprod(double wvx, double wvy, double x, double y):
    '''Fast dot product meeting our exact needs.'''
    return wvx * x + wvy * y


cdef long [:, :] even_adj_indices = n.array(
    [[-1, -1], [-1, 0], [0, 1], [1, 0], [1, -1], [0, -1]])
cdef long [:, :] odd_adj_indices =  n.array(
    [[-1, 0], [-1, 1], [0, 1], [1, 1], [1, 0], [0, -1]])
cdef long [:, :] all_adj_indices = n.array(
    [[-1, -1], [-1, 0], [-1, 1], [0, 1],
     [1, 1], [1, 0], [1, -1], [0, -1]])
cdef long [:, :] safe_adj_indices = n.array(
    [[-1, 0], [0, 1], [1, 0], [0, -1]])

cpdef get_critical_points(double [:, :] arr,
                          to_edges=False,
                          verbose=True):
    '''Cython duplicate of function in neuman.py'''

    cdef long lx = arr.shape[0]
    cdef long ly = arr.shape[1]
    cdef double [:] adjs = n.zeros(6, dtype=n.double)

    cdef list maxima = []
    cdef list minima = []
    cdef list saddles = []
    cdef list degenerate = []

    cdef double [:] border_mult = n.ones(6, dtype=n.double)

    cdef long prevx = -1
    cdef long x, y, i

    cdef long [:, :] ais

    cdef double val

    for x in range(lx):
        for y in range(ly):
            if x != prevx:
                prevx = x
                if verbose:
                    pass
            if to_edges or ((0 < x < lx - 1) and (0 < y < ly - 1)):
                val = arr[x, y]
                if x % 2 == 0:
                    ais = even_adj_indices.copy()
                else:
                    ais = odd_adj_indices.copy()

                for i in range(6):
                    ais[i, 0] += x
                    ais[i, 1] += y
                    if to_edges == 'periodic':
                        ais[i, 0] = ais[i, 0] % lx
                        ais[i, 1] = ais[i, 1] % ly

                for i in range(6):
                    adjs[i] = arr[ais[i, 0] % lx, ais[i, 1] % ly]

                for i in range(6):
                    adjs[i] -= val
                point_type = classify_point(adjs)

                if point_type == 'maximum':
                    maxima.append((x, y))
                elif point_type == 'minimum':
                    minima.append((x, y))
                elif point_type == 'saddle':
                    saddles.append((x, y))
                elif point_type == 'degenerate':
                    degenerate.append((x, y))
                elif point_type == 'fail':
                    print
                    print 'A failure occurred at', x, y
                    print ('=> odd number of sign changes, perhaps the '
                           'function is symmetrical about this point.')

    return (maxima, minima, saddles, degenerate)


cpdef classify_point(double [:] ds):
    cdef double first, second
    cdef long changes

    cdef long i

    for i in range(6):
        if ds[i] <= 0.:
            break
    else:  # Found a place to use for/else!
        return 'minimum'

    for i in range(6):
        if  ds[i] >= 0.:
            break
    else:
        return 'maximum'

    changes = 0
    for i in range(6):
        first = ds[i]
        second = ds[(i+1) % 6]
        if n.sign(first) != n.sign(second):
            changes += 1

    if changes == 2:
        return 'regular'
    elif changes == 4:
        return 'saddle'
    elif changes == 6:
        return 'degenerate'
    else:
        print 'changes', changes
        return 'fail'


cpdef trace_gradient_line(double sx, double sy, double dx, double dy,
                          double xnum, double ynum, func,
                          dict critdict, double [:] start_point,
                          bytes direction, to_edges,
                          tuple func_params=()):
    '''Traces the line of maximal gradient from the given position until
    reaching a critical point or until not really moving any more.'''

    cdef bint to_edges_bool = 1 if to_edges else 0  # to_edges can
                                                    # effectively only be
                                                    # true/false

    cdef long ixnum = n.int64(xnum), iynum = n.int64(ynum)
    cdef double cx = sx, cy = sy
    cdef double startx = start_point[0]
    cdef double starty = start_point[1]

    cdef double dirfac
    if direction == 'down':
        dirfac = -1.
    else:
        dirfac = 1.

    cdef list points = [[cx, cy]]

    cdef tuple gradient
    cdef double angle
    cdef long nearx, neary

    cdef bint use_func_params = 0
    cdef double [:, :] wvs
    cdef double [:] amps
    cdef double [:] phases
    if func_params:
        use_func_params = 1
        params_type, wvs, amps, phases = func_params
    while True:
        if use_func_params:
            gradient = grad_rwm(wvs, amps, phases, startx + cx*dx, starty + cy*dy,
                                dx, dy)
        else:
            gradient = grad(func, startx+cx*dx, starty+cy*dy, dx, dy)
        angle = atan2(gradient[1] * dirfac, gradient[0] * dirfac)

        cx += 0.25*n.cos(angle)
        cy += 0.25*n.sin(angle)

        if len(points) > 20:
            if magdiff(cx, cy, points[-20][0], points[-20][1]) < 0.75:
                return (points, [int(n.round(cx)), int(n.round(cy))])

        points.append([cx, cy])

        if cx < 0 or cx > xnum or cy < 0 or cy > ynum:
            if to_edges_bool:
                cx %= xnum
                cy %= ynum
            else:
                return (points, None)

        nearx = n.int64(round(cx))
        neary = n.int64(round(cy))
        if to_edges_bool:
            nearx %= ixnum
            neary %= iynum
        if (nearx, neary) in critdict:
            crit_type = critdict[nearx, neary]
            if ((crit_type == 'maximum' and direction == 'down') or
                (crit_type == 'minimum' and direction == 'up')):
            # if crit_type in ['maximum','minimum']:
                #print (nearx, neary), crit_type, direction
                points.append((nearx, neary))
                return (points, (nearx, neary))
        else:
            for indices in safe_adj_indices:
                if (nearx + indices[0], neary + indices[1]) in critdict:
                    coords = (nearx + indices[0], neary + indices[1])
                    crit_type = critdict[coords]
                    if ((crit_type == 'maximum' and direction == 'down') or
                        (crit_type == 'minimum' and direction == 'up')):
                    # if crit_type in ['maximum','minimum']:
                        #print (nearx, neary), crit_type, direction
                        points.append(coords)
                        return (points, coords)


cdef inline tuple grad(func, double x, double y, double dx, double dy):
    '''Local gradient of given function at given postion and jump.'''
    cdef double dfdx, dfdy
    dfdx = (func(x, y)-func(x+0.015*dx, y))/(0.015*dx)
    dfdy = (func(x, y)-func(x, y+0.015*dy))/(0.015*dy)
    return dfdx, dfdy

cdef inline tuple grad_rwm(double [:, :] wvs, double [:] amps, double [:] phases,
                     double x, double y, double dx, double dy):
    '''Local gradient of the rwm with the given parameters at the given position
    and jump.'''
    cdef double dfdx, dfdy
    dfdx = (random_wave_function(wvs, amps, phases, x, y) -
            random_wave_function(wvs, amps, phases, x+0.015*dx, y)) / (0.015 * dx)
    dfdy = (random_wave_function(wvs, amps, phases, x, y) -
            random_wave_function(wvs, amps, phases, x, y+0.015*dy)) / (0.015 * dy)
    return dfdx, dfdy

cdef inline magdiff(double a, double b, double c, double d):
    '''Magnitude of the distance from (a, b) to (c, d).'''
    return sqrt(pow(c-a, 2) + pow(d-b, 2))
