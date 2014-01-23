'''Python module for tracing Neumann domains of spherical harmonics,
by subclassing a normal (planar) NeumannTracer.'''

from neumann import NeumannTracer
from functools import partial

from matplotlib import interactive as mpl_interactive
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

try:
    import mayavi.mlab as may
except ImportError:
    print ('Failed to import mayavi. 3d plotting will not work '
           '(but other stuff will work fine).')


class SphericalNeumannTracer(neu.NeumannTracer):
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

    def __init__(self, shape, dn, func):
        self.shape = shape
        self.dn = dn

        self.func = func

        self.create_tracers()
        
    def create_tracers(self):
        '''Creates 6 NeumannTracer objects, one for each cube face.'''
        xnum = self.shape
        ynum = self.shape
        dx = self.dn
        dy = self.dn

        self.tracers = []
        for i in range(6):
            self.tracers.append(SphericalNeumannTracer(xnum, ynum, dx, dy,
                                                       partial(self.cell_side_function, i),
                                                       face=i, handler=self))

    def cell_side_function(self, side, x, y):
        '''Converts x and y to appropriate spherical coordinates before
        calling self.func.'''
        pass  # TODO
