'''
Cython functions for finding Neumann domains.
'''

import numpy as n
cimport numpy as n
int64 = n.int64

cimport cython
from libc.math cimport sin, cos, atan2, sqrt, pow, round

import sys
def lineprint(s='', newline=True):
    sys.stdout.write(s)
    if newline:
        sys.stdout.write('\n')
    sys.stdout.flush()

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


cpdef fill_arr(func, double sx, double sy, double dx, double dy, long xnum,
               long ynum, double [:, :] arr, verbose=True):
    cdef long x, y
    for x in range(xnum):
        if x % 100 == 0 and verbose:
            lineprint('\r\tx = {0} / {1}'.format(x, xnum), False)
        for y in range(ynum):
            arr[x, y] = func(sx + x*dx, sy + y*dy)


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
        if x % 100 == 0 and verbose:
            lineprint('\r\tx = {} / {}'.format(x, lx))
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
    if verbose:
        lineprint()

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

    cdef long ixnum = <long>xnum, iynum = <long>ynum
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
    cdef long index, indx, indy

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

        cx += 0.25*cos(angle)
        cy += 0.25*sin(angle)

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

        nearx = <long>round(cx)
        neary = <long>round(cy)
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
            for index in range(len(safe_adj_indices)):
                indx = safe_adj_indices[index, 0]
                indy = safe_adj_indices[index, 1]
                if (nearx + indx, neary + indy) in critdict:
                    coords = (nearx + indx, neary + indy)
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


cdef class UpsampleTracer:
    '''Simple Neumann tracer that can only find critical points,
    with everything possible implemented in cython for speed.'''

    cdef public double dx, dy, sx, sy
    cdef public double [:, :] arr
    cdef long xnum, ynum
    cdef public list maxima, minima, saddles, degenerate
    cdef public tuple crits
    cdef public dict crits_dict
    cdef object func

    def __init__(self, long xnum, long ynum, double dx, double dy, func,
                   start_point=(0.00123, 0.00123), verbose=False):
        self.arr = n.zeros((xnum, ynum), dtype=n.float64)
        self.xnum = xnum
        self.ynum = ynum
        self.dx = dx
        self.dy = dy
        self.sx = start_point[0]
        self.sy = start_point[1]
        self.func = func
        self.fill_arr()

    cpdef contains(self, double x, double y, double border=0.):
        '''Return True if the x, y point is within self (sx -> sx + xnum*dx etc.),
        else False.'''
        cdef double sx = self.sx
        cdef double sy = self.sy
        cdef double dx = self.dx
        cdef double dy = self.dy
        cdef long xnum = self.xnum
        cdef long ynum = self.ynum
        if ((sx + border < x < sx + xnum*dx - border) and
            (sy + border < y < sy + ynum*dy - border)):
            return True
        return False

    cdef fill_arr(self):
        '''
        Sample the function on a (self.xnum, self.ynum) array.

        Result stored in self.arr.
        '''
        arr = self.arr
        cdef double sx = self.sx
        cdef double sy = self.sy
        cdef double dx = self.dx
        cdef double dy = self.dy
        cdef long xnum = self.xnum
        cdef long ynum = self.ynum
        fill_arr(self.func, sx, sy, dx, dy, xnum, ynum, arr, verbose=False)

    cpdef find_critical_points(self):
        '''
        Iterate over self.arr walking about each pixel and checking the
        number of sign changes. Bins the result appropriately as a
        maximum, minimum, saddle or regular point, and store in
        self.crits.
        '''
        maxima, minima, saddles, degenerate = get_critical_points(
            self.arr, False, False)
            
        self.crits = (maxima, minima, saddles, degenerate)
        self.maxima = self.crits[0]
        self.minima = self.crits[1]
        self.saddles = self.crits[2]
        self.degenerate = self.crits[3]

        self.crits_dict = critical_points_to_index_dict(self.crits)


def critical_points_to_index_dict(crits):
    maxima, minima, saddles, degenerate = crits
    d = {}
    for entry in maxima:
        d[tuple(entry)] = 'maximum'
    for entry in minima:
        d[tuple(entry)] = 'minimum'
    for entry in saddles:
        d[tuple(entry)] = 'saddle'
    for entry in degenerate:
        d[tuple(entry)] = 'degenerate'
    return d

cpdef hessian(func, double x, double y, double dx, double dy):
    cdef tuple grads
    grads = grad(func, x, y, dx, dy)
    dfdx = grads[0]
    dfdy = grads[1]

    dfdxdx = (grad(func, x+0.05*dx, y, dx, dy)[0] - dfdx) / (0.05*dx)
    dfdydy = (grad(func, x, y+0.05*dy, dx, dy)[1] - dfdy) / (0.05*dy)
    dfdxdy = (grad(func, x+0.05*dx, y, dx, dy)[1] - dfdy) / (0.05*dx)
    dfdydx = (grad(func, x, y+0.05*dy, dx, dy)[0] - dfdx) / (0.05*dy)

    return n.array([[dfdxdx, dfdxdy], [dfdydx, dfdydy]])

cpdef hessian_det(func, double x, double y, double dx, double dy):
    cdef double [:, :] hess_mat
    hess_mat = hessian(func, x, y, dx, dy)
    return hess_mat[1, 1] * hess_mat[0, 0] - hess_mat[0, 1] * hess_mat[1, 0]
