'''Python module for tracing Neumann domains of spherical harmonics,
by subclassing a normal (planar) NeumannTracer.'''

from neumann import NeumannTracer
from functools import partial
import numpy as n

from matplotlib import interactive as mpl_interactive
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

try:
    import mayavi.mlab as may
except ImportError:
    print ('Failed to import mayavi. 3d plotting will not work '
           '(but other stuff will work fine).')

import sys
try:
    sys.path.append('/home/asandy')
    from pyknot.rewrite.twosphere import get_random_spherical_harmonic
except ImportError:
    raise ImportError('Error: This module currently depends on c libraries on Sandy\'s computer!')
    

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

    def plot_stereographic(self):
        '''Plot self and lines via stereographic projection.'''
        pass

    def plot_net(self, trace_lines=False):
        '''Plot component tracers via the cube net.'''
        if trace_lines:
            self.trace_lines()
        
        arrs = []
        for tracer in self.tracers:
            arrs.append(n.rot90(tracer.arr[::-1], 3))

        shape = self.shape

        full_arr = n.zeros((4*shape, 3*shape))

        # Top, bottom
        full_arr[shape:(shape+shape), shape:(shape+shape)] = arrs[0]
        full_arr[3*shape:(3*shape+shape), shape:(shape+shape)] = arrs[1]

        # x sides
        full_arr[0:(shape), shape:(shape+shape)] = arrs[4]
        full_arr[2*shape:(2*shape+shape), shape:(shape+shape)] = arrs[5]

        # y sides
        full_arr[shape:(shape+shape), 0:(shape)] = arrs[2]
        full_arr[shape:(shape+shape), 2*shape:(2*shape+shape)] = arrs[3]

        lines = [[[shape, 0],
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
        lines = [n.array(line) - 0.5 for line in lines]
        lines = map(n.array, lines)

        fig, ax = plt.subplots()

        ax.imshow(full_arr, interpolation='none',
                  cmap='RdYlBu_r')

        for line in lines:
            ax.plot(line[:, 0], line[:, 1], color='black',
                    linewidth=2.0)

        ax.set_xlim(0, 3*shape-1)
        ax.set_ylim(0, 4*shape-1)
            
        return fig, ax

    

def side_xy_to_xyz(side, x, y):
    '''Takes a cube side, and x,y values, and returns the x,y,z angles.'''
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

def cube_to_angles(x, y, z):
    '''Takes points on a face of the cube, and maps them to angles on the
    sphere.'''
    r = n.sqrt(x**2 + y**2 + z**2)
    phi = n.arctan2(y, x)
    theta = n.arccos(z/r)
    return theta, phi

def angles_to_stereographic_projection(theta, phi):
    '''Takes angles on the 2-sphere, and gives their x,y positions under
    stereographic projection from z=-1?'''
    chi = theta / 2.  # Projection from south pole
    r = n.tan(chi)
    x = r*n.sin(theta)*n.cos(phi)
    y = r*n.sin(theta)*n.sin(phi)
    return x, y
