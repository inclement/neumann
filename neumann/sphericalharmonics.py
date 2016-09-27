'''Python module for tracing Neumann domains of spherical harmonics,
by subclassing a normal (planar) NeumannTracer.'''

from neumann.neumann import NeumannTracer
from functools import partial
import numpy as n

import matplotlib.pyplot as plt

try:
    import mayavi.mlab as may
except ImportError:
    pass
    # print('Failed to import mayavi. 3d plotting will not work '
    #       '(but other stuff will work fine).')


class SphericalNeumannTracer(NeumannTracer):
    '''(internal)

    A NeumannTracer subclass with modified boundary behavior; instead
    of wrapping periodically, it calls back to the parent to continue
    tracing in an adjacent tracer.

    This should be used *only* with a NeumannCubeTracer!'''

    def __init__(self, *args, **kwargs):
        self.face = kwargs.pop('face')
        self.handler = kwargs.pop('handler')
        super(SphericalNeumannTracer, self).__init__(*args, **kwargs)


class NeumannCubeHandler(object):
    '''Wrapper class for tracing Neumann lines of spherical harmonics via
    a coordinate system based on stereographic projection of the
    hypercube. Automatically takes care of wrapping the boundaries of
    NeumannTracers on each face.

    Args:

    * shape: xnum,ynum for NeumannTracer (these are passed through).
    * dn: dx,dy for NeumannTracer
    * func: Spherical harmonic function of theta and phi.

    Other NeumannTracer arguments are generated internally as appropriate.
    '''

    def __init__(self, shape=10, func=None):
        if func is None:
            raise Exception('func must not be None')

        self.shape = shape

        self.func = func

        # Line tracing milestones
        self.traced_lines = False

        self.create_tracers()

    def create_tracers(self):
        '''Creates 6 NeumannTracer objects, one for each cube face.'''
        xnum = self.shape
        ynum = self.shape
        dx = 1. / xnum
        dy = 1. / ynum

        self.tracers = []
        for i in range(6):
            self.tracers.append(
                SphericalNeumannTracer(xnum, ynum, dx, dy,
                                       partial(self.cell_side_function, i),
                                       face=i, handler=self))

    def find_critical_points(self):
        '''Ask each tracer in self.tracers to find its critical points.'''
        for tracer in self.tracers:
            tracer.find_critical_points()

    def trace_lines(self):
        '''Ask each tracer in self.tracers to find its Neumann lines. Ignores
        boundaries for the time being.'''
        for tracer in self.tracers:
            tracer.trace_neumann_lines()

    def cell_side_function(self, side, x, y):
        '''Converts x and y to appropriate spherical coordinates before
        calling self.func.'''
        x = x - 0.5  # Individual tracers are in 0-1 range
        y = y - 0.5
        position = side_xy_to_xyz(side, x, y)
        theta, phi = cube_to_angles(*position)
        return self.func(theta, phi)

    def plot_3d_mayavi(self, stereographic=True, plot_criticals=False,
                       trace_lines=False, clf=True):
        '''Plot self and lines via stereographic projection.'''
        if clf:
            may.clf()
        if plot_criticals:
            self.find_critical_points()
        if trace_lines:
            self.trace_lines()

        arrs = [n.rot90(tracer.arr[::-1], 3)*-1 for tracer in self.tracers]

        arr_max = n.max([n.max(n.abs(arr)) for arr in arrs])

        shape = self.shape

        if stereographic:
            for side, arr in enumerate(arrs):
                arr[0, 0] = arr_max
                arr[0, 1] = -1*arr_max
                if side != 1:
                    xs, ys, zs = side_to_xs_ys_zs(side, shape)
                    stxs, stys = vcube_to_stereographic_projection(xs, ys, zs)
                    real_zs = n.zeros(xs.shape)
                    may.mesh(stxs, stys, real_zs, scalars=arr*-1,
                             colormap='RdYlBu')

            if plot_criticals:
                crit_sets = [tracer.crits for tracer in self.tracers]
                for side, crit_set in enumerate(crit_sets):
                    maxima, minima, saddles, degenerate = crit_set
                    maxima = crits_to_stereographic_projection(
                        side, shape, maxima)
                    print('maxima are', maxima)

                    if len(maxima) > 0:
                        scalars = 1. * n.ones(len(maxima))
                        may.points3d(maxima[:, 0], maxima[:, 1],
                                     n.zeros(len(maxima)), scalars,
                                     color=(1, 0, 0))

        else:  # not stereographic, plot full sphere
            for side, arr in enumerate(arrs):
                arr[0, 0] = arr_max
                arr[0, 1] = -1*arr_max
                xs, ys, zs = side_to_xs_ys_zs(side, shape)
                stxs, stys, stzs = vcube_to_sphere(xs, ys, zs)
                may.mesh(stxs, stys, stzs, scalars=arr*-1, colormap='RdYlBu')

            if plot_criticals:
                crit_sets = [tracer.crits for tracer in self.tracers]
                for side, crit_set in enumerate(crit_sets):
                    maxima, minima, saddles, degenerate = crit_set
                    print('maxima are', maxima)
                    maxima = crits_to_sphere(side, shape, maxima)
                    print('maxima are', maxima)

                    if len(maxima[0]) > 0:
                        print('plotting')
                        may.points3d(maxima[0], maxima[1], maxima[2],
                                     color=(1, 0, 0))

    def plot_3d_vispy(self, stereographic=True, plot_criticals=False,
                      trace_lines=False,
                      clf=True):
        '''Plot self and lines via stereographic projection.'''
        from pyknot2 import visualise as pvis
        pvis.ensure_vispy_canvas()
        import vispy.scene as vs
        if clf:
            pvis.clear_vispy_canvas()
        if plot_criticals:
            self.find_critical_points()
        if trace_lines:
            self.trace_lines()

        arrs = [n.rot90(tracer.arr[::-1], 3)*-1 for tracer in self.tracers]

        arr_max = n.max([n.max(n.abs(arr)) for arr in arrs])

        shape = self.shape

        if stereographic:
            for side, arr in enumerate(arrs):
                arr[0, 0] = arr_max
                arr[0, 1] = -1*arr_max
                if side != 1:
                    xs, ys, zs = side_to_xs_ys_zs(side, shape)
                    stxs, stys = vcube_to_stereographic_projection(
                        xs, ys, zs)
                    real_zs = n.zeros(xs.shape)
                    mesh = vs.GridMesh(stxs, stys, real_zs)
                    pvis.vispy_canvas.view.add(mesh)
                    # may.mesh(stxs, stys, real_zs, scalars=arr*-1,
                    #          colormap='RdYlBu')

            if plot_criticals:
                raise ValueError('criticals not supported yet')
                crit_sets = [tracer.crits for tracer in self.tracers]
                for side, crit_set in enumerate(crit_sets):
                    maxima, minima, saddles, degenerate = crit_set
                    maxima = crits_to_stereographic_projection(
                        side, shape, maxima)
                    print('maxima are', maxima)

                    if len(maxima) > 0:
                        scalars = 1. * n.ones(len(maxima))
                        may.points3d(maxima[:, 0], maxima[:, 1],
                                     n.zeros(len(maxima)), scalars,
                                     color=(1, 0, 0))

        else:  # not stereographic, plot full sphere
            for side, arr in enumerate(arrs):
                arr[0, 0] = arr_max
                arr[0, 1] = -1*arr_max
                xs, ys, zs = side_to_xs_ys_zs(side, shape)
                stxs, stys, stzs = vcube_to_sphere(xs, ys, zs)
                mesh = vs.GridMesh(stxs, stys, stzs)
                pvis.vispy_canvas.view.add(mesh)
                # may.mesh(stxs, stys, stzs, scalars=arr*-1, colormap='RdYlBu')

            if plot_criticals:
                raise ValueError('criticals not yet supported')
                crit_sets = [tracer.crits for tracer in self.tracers]
                for side, crit_set in enumerate(crit_sets):
                    maxima, minima, saddles, degenerate = crit_set
                    print('maxima are', maxima)
                    maxima = crits_to_sphere(side, shape, maxima)
                    print('maxima are', maxima)

                    if len(maxima[0]) > 0:
                        print('plotting')
                        may.points3d(maxima[0], maxima[1], maxima[2],
                                     color=(1, 0, 0))

        pvis.vispy_canvas.view.camera = vs.ArcballCamera(fov=30)
        pvis.vispy_canvas.show()
            

    def plot_net(self, plot_criticals=True, trace_lines=False, cmap='RdYlBu_r',
                 normalise=True):
        '''Plot component tracers via the cube net.'''
        if plot_criticals:
            self.find_critical_points()
        if trace_lines:
            self.trace_lines()
        
        arrs = [n.rot90(tracer.arr[::-1], 3) for tracer in self.tracers]

        shape = self.shape

        full_arr = n.zeros((4*shape, 3*shape))

        # Top, bottom
        full_arr[shape:(shape+shape), shape:(shape+shape)] = arrs[0]
        full_arr[3*shape:(3*shape+shape), shape:(shape+shape)] = arrs[1]

        # x sides
        full_arr[0:(shape), shape:(shape+shape)] = arrs[4]
        full_arr[2*shape:(2*shape+shape), shape:(shape+shape)] = arrs[5]

        # y sides
        full_arr[shape:(shape+shape), 0:(shape)] = arrs[3]
        full_arr[shape:(shape+shape), 2*shape:(2*shape+shape)] = arrs[2]


        abs_arr = n.abs(full_arr)
        if normalise:
            full_arr[0, 0] = n.max(abs_arr)
            full_arr[0, 1] = full_arr[0, 0] * -1

        boundary_lines = [[[shape, 0],
                          [shape, shape],
                          [0, shape]],
                         [[2*shape, 0],
                          [2*shape, shape],
                          [3*shape, shape]],
                         [[0, 2*shape],
                          [shape, 2*shape],
                          [shape, 4*shape]],
                         [[3*shape, 2*shape],
                          [2*shape, 2*shape],
                          [2*shape, 4*shape]]]
        boundary_lines = [n.array(line) - 0.5 for line in boundary_lines]

        fig, ax = plt.subplots()

        ax.imshow(full_arr, interpolation='none',
                  cmap=cmap, alpha=1.0)

        if plot_criticals:
            crit_sets = [tracer.crits for tracer in self.tracers]
            for side, crit_set in enumerate(crit_sets):
                maxima, minima, saddles, degenerate = crit_set
                maxima = net_representation_shift(side, shape, maxima)
                minima = net_representation_shift(side, shape, minima)
                saddles = net_representation_shift(side, shape, saddles)

                if len(maxima) > 0:
                    ax.scatter(maxima[:, 0], maxima[:, 1], 60, c='r')
                if len(minima) > 0:
                    ax.scatter(minima[:, 0], minima[:, 1], 60, c='b')
                if len(saddles) > 0:
                    ax.scatter(saddles[:, 0], saddles[:, 1], 60, color='yellow')

        for line in boundary_lines:
            ax.plot(line[:, 0], line[:, 1], color='black',
                    linewidth=2.0)

        if trace_lines:
            line_sets = [tracer.lines for tracer in self.tracers]
            for side, lines in enumerate(line_sets):
                print(side)
                for line in lines:
                    shifted_line = net_representation_shift(side, shape, line)
                    ax.plot(shifted_line[:, 0], shifted_line[:, 1], '-',
                            color='purple')
                    
        ax.set_xticks([])
        ax.set_yticks([])

        ax.set_xlim(0, 3*shape-1)
        ax.set_ylim(0, 4*shape-1)

        fig.tight_layout()

        fig.show()
            
        return fig, ax

def net_representation_shift(side, shape, line):
    '''Given a side, a shape, and a line from that side's array, shifts
    the line in the x-y plane so that it will be in the right place in
    a net diagram.
    '''
    line = n.array(line)
    if side == 0:
        line[:, 0] += shape 
        line[:, 1] += shape
    elif side == 1:
        line[:, 0] += shape
        line[:, 1] += 3*shape
    elif side == 2:
        line[:, 0] += 0
        line[:, 1] += shape
    elif side == 3:
        line[:, 0] += 2*shape
        line[:, 1] += shape
    elif side == 4:
        line[:, 0] += shape
        line[:, 1] += 0
    elif side == 5:
        line[:, 0] += shape 
        line[:, 1] += 2*shape 
    return line

def side_xy_to_xyz(side, x, y):
    '''Takes a cube side, and x,y values, and returns the x,y,z positions.'''
    if side == 0:
        return x, y, 0.5
    elif side == 1:
        return x, y, -0.5
    elif side == 2:
        return 0.5, -y, -x
    elif side == 3:
        return -0.5, y, x 
    elif side == 4:
        return x, 0.5, y
    elif side == 5:
        return x, -0.5, -y
vside_xy_to_xyz = n.vectorize(side_xy_to_xyz)

def cube_to_angles(x, y, z):
    '''Takes points on a face of the cube, and maps them to angles on the
    sphere.'''
    r = n.sqrt(x**2 + y**2 + z**2)
    phi = n.arctan2(y, x)
    theta = n.arccos(z/r)
    return theta, phi

def angles_to_sphere(theta, phi):
    '''Takes angles on the 2-sphere, and gives their x,y positions on that
    sphere.
    '''
    return n.sin(theta)*n.cos(phi), n.sin(theta)*n.sin(phi), n.cos(theta)




def angles_to_stereographic_projection(theta, phi):
    '''Takes angles on the 2-sphere, and gives their x,y positions under
    stereographic projection from z=-1?'''
    chi = theta / 2.  # Projection from south pole
    x = 2*n.tan(chi)*n.cos(phi)
    y = 2*n.tan(chi)*n.sin(phi)
    return x, y

def crits_to_stereographic_projection(side, shape, crits):
    scaled_crits = []
    for crit in crits:
        scaled_crits.append(n.array(crit)/float(shape) - 0.5)
    print('crits', crits)
    print('scaled', scaled_crits)
    return(side_line_to_stereographic_projection(side, scaled_crits))

def crits_to_sphere(side, shape, crits):
    scaled_crits = []
    for crit in crits:
        scaled_crits.append(n.array(crit)/float(shape) - 0.5)
    xyzs = []
    for crit in scaled_crits:
        xyzs.append(side_xy_to_xyz(side, crit[0], crit[1]))
    xyzs = n.array(xyzs)
    print('xyzs are', xyzs)

    print('xyzs are', xyzs)
    sphere_xyzs = vcube_to_sphere(xyzs[:, 0], xyzs[:, 1], xyzs[:, 2])
    print('and on sphere', sphere_xyzs)
    return sphere_xyzs

def side_line_to_stereographic_projection(side, line):
    points = []
    for point in line:
        points.append(
            side_xy_to_stereographic_projection(side, point[0], point[1]))
    return n.array(points)

def side_line_to_sphere(side, line):
    points = []
    for point in line:
        points.append(side_xy_to_sphere(side, point[0], point[1]))
    return n.array(points)

def side_line_to_stereographic_projection(side, line):
    points = []
    for point in line:
        points.append(
            side_xy_to_stereographic_projection(side, point[0], point[1]))
    return n.array(points)

def side_xy_to_stereographic_projection(side, x, y):
    x, y, z = side_xy_to_xyz(side, x, y)
    theta, phi = cube_to_angles(x, y, z)
    px, py = angles_to_stereographic_projection(theta, phi)
    return [px, py]

def side_to_xs_ys_zs(side, shape):
    xs, ys = n.mgrid[-0.5:0.5:shape*1j, -0.5:0.5:shape*1j]
    zs = n.ones(xs.shape) * 0.5

    if side == 0:
        return xs, ys, zs
    elif side == 1:
        return xs, ys, -1*zs
    elif side == 5:
        return zs, ys, -1*xs
    elif side == 4:
        return -1*zs, ys, xs
    elif side == 2:
        return -1*xs, zs, -1*ys
    elif side == 3:
        return xs, -1*zs, ys

def cube_to_stereographic_projection(x, y, z):
    angles = cube_to_angles(x, y, z)
    return angles_to_stereographic_projection(*angles)
vcube_to_stereographic_projection = n.vectorize(
    cube_to_stereographic_projection)

def cube_to_sphere(x, y, z):
    angles = cube_to_angles(x, y, z)
    return angles_to_sphere(*angles)
vcube_to_sphere = n.vectorize(cube_to_sphere)
    
