'''Equations
=========

.. math::
    \\psi(x,y) = \\sum_n^N a_n \\sin\\left( \\frac{2\\pi}{\\sqrt{|k|^2/2}}
     \\boldsymbol{k_n}\\cdot\\boldsymbol{r} + \\theta_n \\right)

Where:

* :math:`a_n` is a gaussian random amplitude.
* :math:`\\boldsymbol k_n` is a random wavevector with *integer* coefficients,
   and where all :math:`k_n` have the same magnitude.
* :math:`\\theta_n` is a random phase.

An 'energy' parameter in the code is always :math:`|k|^2`, and does
*not* include the factors of :math:`2\pi`. This is the root of the
normalisation problems!

The funny prefactor (:math:`(2\pi)/\sqrt{|k|^2/2}`) is to enforce
periodicity when generating a torus eigenfunction. It means that as
any component of
:math:`r` runs from :math:`0` to :math:`\sqrt{|k|^2/2}`, the argument
of :math:`\sin` runs
from :math:`0` to :math:`1\\times` some integer, where the integer is one of the
components of :math:`k`.

Code
====

Wherever there is a scale or energy parameter, this refers to
:math:`|k|^2`, and does *not* include the prefactor of :math:`2\pi` in
the equation.

Everything here assumes you already did::

    import numpy as n
    import neumann as neu

You can see the documentation of any python function or class using
`help(functionname)`. Not all of my code is documented this way, but
some of it is. Any functions with docstrings are also documented on
this page.

Torus eigenfunctions
--------------------

To get a torus eigenfunction::

    a = neu.get_periodic_tracer(17, downscale=3)

The first argument (`17`) is the scale :math:`|k|^2`.

The second argument
(`downscale=3`) controls the numerical resolution of the sampled
function - this is automatically increased as the energy increases,
so that the number of samples per domain remains roughly
constant. You can mostly ignore this argument entirely.

The function returns a :class:`NeumannTracer` object, which is a python class
implementing all the critical point detection, domain area
calculation, degree counting etc.

Random wave models
------------------

You can sample a section of a true random wave (no integer limitation on
wavevectors) with::

    func = neu.random_wave_function(number=50, mag=10)
    a = neu.NeumannTracer(100, 100, n.pi/50, n.pi/50, func)

The first line generates a function (func) describing a random wave
with 'number' superposed sine waves of wavevector magnitude 'mag'.

The second line creates a NeumannTracer object. The first two
parameters (100, 100) are the X and Y sizes of the sampled array. The
second two numbers (n.pi/50, n.pi/50) are the step size when sampling
the random wave function. The final argument (func) is obviously the
function generated in the first line.

This isn't a very neat interface to creating random wave functions,
and there's currently no equivalent of the hand 'get_periodic_tracer'
function for generating eigenfunctions. I may add one.

Statistics
----------

Once you have a NeumannTracer (e.g. from neu.get_periodic_tracer),
you can retrieve various statistics from it.

    a = neu.get_periodic_tracer(17, downscale=3)

If you now type a and press tab, you can see all the available
methods of a (this is one of ipython's nice features). Most are not
important, but the ones named 'get_...' retrieve interesting
statistics. Specifically:

- `a.get_domain_areas()` returns a list
  of all recognised domain areas.
- `a.get_domain_perimeters()` returns a list
  of all recognised domain perimeters.
- `a.get_domain_rhos()` returns a list
  of the dimensionless rho parameter for each domain.
- `a.get_critical_degrees()` returns a tuple
  of lists, the first containing a list of degrees of maxima, and the
  second containing a list of degrees of minima.
- `a.get_critical_degree_dists()` returns an array
  of critical degrees along with the fraction of critical points with
  this number.

At the time of writing, the areas/perimeters/rhos are *not*
normalised properly.

Compilation
-----------

The code now includes a cython version of some routines, which is a lot faster
(even though the rwm functions were already using numpy).

If cneumann.pyx exists and is compiled (to cneumann.so), the main code
will automatically load and try to use it. However, if the cython code isn't
available, it will fall back to pure python rather than failing (even if
you specify :code:`compiled=True` to functions that take this argument).

Plotting
========

Basic domain plots
------------------

Once you have a NeumannTracer, you can plot it
with `a.plot()`. This creates a basic
visualisation showing the critical points and Neumann lines.

You can do `help(a.plot)` to see the
available arguments::

    plot(self, trace_lines=True, plot_hessian=False,
         show_saddle_directions=False, show_domain_patches=False,
         print_patch_areas=False, figsize=None, show_sample_directions=False,
         save=False, figax=None)

You can toggle most of these to see what they do. I mostly use just
the basic plot (using the default arguments, so just `a.plot()`) or
`a.plot(show_domain_patches=True, print_patch_areas=True)` which plots
the different colours for each domain along with the areas.

You can also now set some more options; :code:`maxima_style`,
:code:`minima_style` and :code:`saddle_style`. These should be
dictionaries that may contain any normal matplotlib code

Histogram of :math:`\\rho`
------------------------

First, this assumes you have some data generated by :func:`get_statistics_at` or
similar::

    results = neu.get_statistics_at([50], 200)

:math:`\\rho` data is the third entry in `results[50]`, and is returned as a
list. To plot a histogram, you can do::

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots() # Do help(plt.subplots) to see options
                             # for (e.g.) multiple plots in one figure
    ax.histogram(results[50][2], bins=20) # results[50][2] is the rho data

Here, the `bins` is the number of bins to use (you can also supply a
list of bin minima). You can do `help(ax.histogram)` or check the
matplotlib docs online to see the (many) other options. One useful
addition is `normed=True`, which normalises the histogram.

Setting the axes etc. is done via::

    ax.set_xlabel('x label')
    ax.set_ylabel('y label')

You can use LaTeX code for equations by enclosing it in single dollar signs.

Save the figure using::

    fig.set_size_inches((x,y)) # The x and y size in inches
    fig.savefig(filename, dpi=300) # dpi doesn't need to be that high

Degree distributions
--------------------

As above, this assumes you have some data generated by
:func:`get_statistics_at` or similar::

    results = neu.get_statistics_at([50], 200)

You can also retrieve results from :func:`save_results_at`, but I'll write
how later.

The degree statistics are the fourth entry in results[50], so you can do::

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    degrees = results[3]
    ax.plot(degrees[:,0], degrees[:,1])  # See `help(ax.plot)` for many
                                         # other options

This makes a simple line plot with degree on the x-axis and proportion on
the y-axis.


Plotting arbitrary functions
----------------------------

You can plot any function you like by defining a python function::

    import numpy as n
    def my_function(x, y): # x and y will be the arguments
        return n.sin(4*x) * n.cos(2*y)  # Returning an example
                                        # trigonometric function

To trace the Neumann lines, you just pass your function into a
:class:`NeumannTracer`::

    import neumann as neu
    tracer = neu.NeumannTracer(100, 100, n.pi/50, n.pi/50, func)

Or, if the function is periodic over the input domain (see below), you
should do::

    tracer = neu.NeumannTracer(100, 100, n.pi/50, n.pi/50, func,
                               to_edges='periodic')

The first two arguments (100 and 100) are the x and y
sampling. The second two arguments (n.pi/50) are the dx and dy values
to sample on. For instance, the example code would have a domain of
:math:`0` to :math:`2\\pi` :math:`(=100*\\pi/50)`.

The tracer then just works as normal, so you can to (for instance)::

    tracer.get_domain_areas()

to get a list of Neumann domain areas in the function you passed over the
range of the function.




Module documentation
====================

The individual classes and functions of neumann.py are documented below.

'''

import numpy as n
from itertools import product, chain, permutations
from colorsys import hsv_to_rgb

from matplotlib import interactive as mpl_interactive
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.cm import jet, hsv

try:
    from scipy.spatial import (Voronoi, voronoi_plot_2d,
                               Delaunay, delaunay_plot_2d)
    scipy_spatial_import = True
except ImportError:
    scipy_spatial_import = False
    print 'Failed to import Delaunay and Voronoi tools.'
from scipy.misc import factorial

import random
import os
import sys
import cPickle
from functools import partial

try:
    import cneumann as cneu
except ImportError:
    print ('Failed to import cneumann. Everything will work fine, but if '
           'fixed this will make things much faster!')
    cneu = None
    
# try:
#     import mayavi.mlab as may
# except ImportError:
#     may = None
#     print 'Failed to import mayavi. 3d plotting will not work.'

try:
    import igraph as ig
except ImportError:
    ig = None
    print 'Failed to import igraph. Using igraph tools may crash the program.'

mpl_linestyles = ['', ' ', 'None', '--', '-.', '-', ':']
patch_linestyles = ['solid','dashed','dashdot','dotted']


# Some default styles for plotting
# These may be passed as options in the NeumannTracer plot function
maxima_style_old = {'c': 'r'}
minima_style_old = {'c': 'b'}
saddle_style_old = {'color': 'yellow'}
saddle_style_rami = {'color': 'purple', 'marker': 'd'}
saddle_style_sandy = {'color': 'green', 'marker': 'd'}


def rotation_matrix(angle):
    return n.array([[n.cos(angle), -n.sin(angle)],
                    [n.sin(angle), n.cos(angle)]])

def angle_between_vectors(a, b):
    '''Returns the angle between the two vectors, calculated via the dot
    product.'''
    return n.arccos(a.dot(b) / (mag(a)*mag(b)))


class CriticalGraph(dict):
    '''A dict subclass representing the graph of Neumann nodes (critical
    points) and their connections (Neumann gradient lines). Nodes are
    stored via their coordinates, with node information as [node type,
    [list of Neumann lines]].

    Nodes should be added with add_or_edit_node rather than by direct
    dict setting.

    '''
    def __init__(self, xnum=0, ynum=0, dr=1.0, *args):
        '''Arguments xnum and ynum should be the array size in each direction,
        and are used when calculating angles between critical
        points.

        dr is the length normalisation, applied when retrieving any
        distance quantity.

        '''
        super(CriticalGraph, self).__init__(*args)
        self.xnum = xnum
        self.ynum = ynum
        self.dr = dr  # Assumes equal dx, dy!
        self.nodes_aligned = False
        self.closed_domains = []
        self.got_closed_domains = False

    def dual_graph(self):
        '''Returns the dual graph of self, where neumann domains are the new
        vertices, and edges exist between domains bounding one another.'''

        self.align_nodes()
        nodes = self.keys()

        labels_by_crit = {}
        for node in nodes:
            labels_by_crit[node] = []

        # For each node, label domains as the clockwise neighbours of
        # incoming lines
        current_index = 0
        for critical_point, data in self.iteritems():
            for line in data[1]:
                labels_by_crit[critical_point].append(current_index)
                current_index += 1

        # Build a dictionary of equalities between node labels
        equivalent_labels = {}
        for critical_point, data in self.iteritems():
            lines = data[1]
            for i in range(len(lines)):
                line = lines[i]
                current_node_label = labels_by_crit[critical_point][i]

                other_end = line.end
                other_line_index = self.incoming_index(other_end,
                                                       critical_point)
                other_possibilities = labels_by_crit[other_end]
                matching_index = (other_line_index-1) % len(other_possibilities)
                other_node_label = other_possibilities[matching_index]

                if current_node_label not in equivalent_labels:
                    equivalent_labels[current_node_label] = set()
                if other_node_label not in equivalent_labels:
                    equivalent_labels[other_node_label] = set()
                equivalent_labels[current_node_label].add(other_node_label)
                equivalent_labels[other_node_label].add(current_node_label)

        for i in range(10):
            for label, equivs in equivalent_labels.iteritems():
                equivs.add(label)
                for equiv in equivs:
                    for equiv2 in equivs:
                        equivalent_labels[equiv].add(equiv2)

        real_labels = {}
        # Replace all the labels with one number
        for label, equivs in equivalent_labels.iteritems():
            if real_labels.has_key(label):
                pass
            else:
                for equiv in equivs:
                    real_labels[equiv] = label

        real_labels_by_crit = {}
        for label, domains in labels_by_crit.iteritems():
            real_labels_by_crit[label] = map(lambda j: real_labels[j], domains)

        # Create the connections list
        connections = {}
        for label, reallabel in real_labels.iteritems():
            if connections.has_key(reallabel):
                pass
            else:
                connections[reallabel] = set()

        for crit, domains in real_labels_by_crit.iteritems():
            for i in range(len(domains)):
                prev_domain = domains[(i-1) % len(domains)]
                cur_domain = domains[i]
                next_domain = domains[(i+1) % len(domains)]
                if next_domain != cur_domain:
                    connections[cur_domain].add(next_domain)
                if prev_domain != cur_domain:
                    connections[cur_domain].add(prev_domain)

        return connections

        # Build a dictionary of connections between the new nodes

        # Build an igraph representation

    def add_or_edit_node(self, node_coords, node_type, new_line):
        '''
        Add the line new_line to the node with coords node_coords. If this
        node does not already exist, create it first.

        node_type is the critical point type, i.e. maximum, minimum or saddle.
        '''
        assert isinstance(new_line, NeumannLine), "Didn't receive NeumannLine!"
        if not self.has_key(node_coords):
            self[node_coords] = [node_type, []]
        self[node_coords][1].append(new_line)

    def align_nodes(self):
        '''Goes through every node of the graph and rearranges the lines so
        that they are listed in increasing order of approach angle about the
        node.

        This is called automatically (and is essential) before tracing
        neumann domains.
        '''
        for key in self:
            lines = self[key][1]
            angles = n.zeros(len(lines), dtype=n.float64)
            for i in range(len(lines)):
                line = lines[i]
                xj, yj = line.jumps()
                first = line[0]
                second = line[-1]  # Strictly could fail in *very*
                                   # highly pathological case
                distance = second - first
                if xj:
                    if distance[0] > 0:
                        distance[0] -= self.xnum
                    else:
                        distance[0] += self.xnum
                if yj:
                    if distance[1] > 0:
                        distance[1] -= self.ynum
                    else:
                        distance[1] += self.ynum
                angle = n.arctan2(distance[1], distance[0])
                angles[i] = angle
            order = n.argsort(angles)
            newlines = []
            for i in range(len(order)):
                newlines.append(lines[order[i]])
            self[key][1] = newlines
        self.nodes_aligned = True

    def get_closed_domains(self, recreate=False):
        '''
        Returns a list of all closed domains, recognised by walking
        clockwise or anticlockwise from every saddle point.

        Explicitly fails when a node passes through a cell boundary or
        includes a line that returns to its own node.
        '''
        if not self.nodes_aligned:
            self.align_nodes()
        if self.got_closed_domains and not recreate:
            return self.closed_domains
        blighted_starts = []
        domains = []
        for start in self:
            crit_type, lines = self[start]
            if crit_type == 'saddle' and start not in blighted_starts:
                blighted_starts.append(start)
                for line in self[start][1]:
                    dom = self.get_domain_from(line)
                    if dom is not None:
                        domains.append(dom)
                        # for nl in dom.lines:
                        #     if nl.end not in blighted_starts:
                        #         blighted_starts.append(nl.end)

        domain_keys = []
        unique_domains = []
        for domain in domains:
            vertices = domain.vertices()
            if vertices not in domain_keys:
                domain_keys.append(vertices)
                unique_domains.append(domain)

        self.closed_domains = unique_domains
        self.got_closed_domains = True

        return unique_domains

    def get_domain_areas(self):
        domains = self.get_closed_domains()
        areas = []
        i = 0
        for domain in domains:
            if i % 50 == 0:
                lineprint('\rGetting area of domain {} / {}'.format(
                    i, len(domains)),
                          False)
            i += 1
            area = domain.area()
            areas.append(area)
        print # Lineprint newline
        return areas

    def get_domain_perimeters(self):
        domains = self.get_closed_domains()
        perimeters = []
        i = 0
        for domain in domains:
            if i % 50 == 0:
                lineprint('\rGetting perimeter of domain {} / {}'.format(
                    i, len(domains)),
                          False)
            i += 1
            perimeter = domain.perimeter()
            perimeters.append(perimeter)
        print  # Lineprint newline
        return perimeters

    def get_domain_diameters(self):
        domains = self.get_closed_domains()
        diameters = []
        i = 0
        for domain in domains:
            if i % 50 == 0:
                lineprint('\rGetting diameter of domain {} / {}'.format(
                    i, len(domains)),
                          False)
            i += 1
            diameter = domain.diameter()
            if diameter is not None:
                diameters.append(diameter)
        print # Lineprint newline
        return diameters

    def get_domain_rhos(self, eigenvalue=None):
        domains = self.get_closed_domains()
        rhos = []
        i = 0
        for domain in domains:
            if i % 50 == 0:
                lineprint('\rGetting rho of domain %d / %d' % (i, len(domains)),
                          False)
            i += 1
            rho = domain.rho()
            rhos.append(rho)
        print  # Lineprint newline
        if eigenvalue is not None:
            root_eig = n.sqrt(eigenvalue)
            rhos = [rho * root_eig for rho in rhos]
        return rhos

    def get_domain_rhos_by_type(self, eigenvalue=None):
        domains = self.get_closed_domains()
        rhos = {'lens': [], 'star': [], 'wedge': [], 'bad': []}
        i = 0
        for domain in domains:
            if i % 50 == 0:
                lineprint('\rGetting rho of domain %d / %d' % (i, len(domains)),
                          False)
            i += 1
            rho = domain.rho()
            domain_type = domain.guess_type()
            if domain_type is None:
                domain_type = 'bad'
            rhos[domain_type].append(rho)
        print  # Lineprint newline
        if eigenvalue is not None:
            root_eig = n.sqrt(eigenvalue)
            for key, rho_vals in rhos.items():
                rhos[key] = [rho * root_eig for rho in rho_vals]
        return rhos

    def get_domain_from(self, line, dir='clockwise'):
        '''
        Returns the closed domain (or None if the algorithm fails)
        obtained by walking along line until reaching the start point
        of line again.
        '''
        if not self.nodes_aligned:
            self.align_nodes()
        lines = []
        lines.append(line)
        closing_node = line.start
        start = line.start
        end = line.end
        curline = line
        if end == start:
            return None
        while end != closing_node:
            if len(lines) > 10 or end is None or end == start:
                return None
            node = curline.end
            crit_type, node_lines = self[node]
            lines_come_from = map(lambda j: j.end, node_lines)
            if start in lines_come_from:
                in_line_index = lines_come_from.index(start)
            else:
                print 'Start not found'
                return None
            if dir == 'clockwise':
                out_line_index = (in_line_index + 1) % len(node_lines)
            else:
                out_line_index = (in_line_index - 1) % len(node_lines)
            curline = node_lines[out_line_index]
            lines.append(curline)
            start = curline.start
            end = curline.end
        return NeumannDomain(lines, self.dr)

    def get_crit_degree_dists(self):
        '''Returns the dimension (number of lines going in/out) of each
        critical point in self. Order is maxima, minima, saddles.

        '''
        sadnums = []
        maxnums = []
        minnums = []
        for node in self:
            #print self[node]
            crit_type, lines = self[node]
            if crit_type == 'saddle':
                sadnums.append(len(lines))
            elif crit_type == 'maximum':
                maxnums.append(len(lines))
            elif crit_type == 'minimum':
                minnums.append(len(lines))
        return maxnums, minnums, sadnums

    def incoming_index(self, node, coords):
        '''Returns the index (if any) of incoming lines from coords to
        node.'''
        lines = self[node][1]
        matches = filter(lambda j: j.end == coords, lines)
        if len(matches) == 1:
            return lines.index(matches[0])
        elif len(matches) > 1:
            return lines.index(matches[0])
            raise Exception('Two incoming lines from the same node!?')
        else:
            return None

class DelaunayGraph(object):
    '''Takes a set of points, and provides some basic methods for
    working with their Delaunay triangulation.'''
    def __init__(self, points, types, boundary=None):
        assert len(types) == len(points)
        self.triangulation = triangulation = Delaunay(n.array(points))
        self.graph = graph = {}

        points = triangulation.points.astype(n.int64)
        points = [tuple(point) for point in points]
        simplices = triangulation.simplices

        for simplex in simplices:
            simplex_points = n.array([points[i] for i in simplex])
            for i in range(3):
                cur_roll = n.roll(simplex_points, i, axis=0)
                cur_crit = tuple(cur_roll[0])
                nex_crit = tuple(cur_roll[1])
                cur_type = types[cur_crit]
                nex_type = types[nex_crit]
                if (cur_type == 'saddle' and nex_type == 'saddle' or
                    cur_type == 'minimum' and nex_type == 'maximum' or
                    cur_type == 'maximum' and nex_type == 'minimum'):
                    continue  
                if cur_crit not in graph:
                    graph[cur_crit] = set()
                adjacency_set = graph[cur_crit]
                adjacency_set.add(nex_crit)

    def mark_boundary_points(self):
        pass

    def get_degree_dists(self):
        degrees = {}
        for key, value in self.graph.items():
            num = len(value)
            if num not in degrees:
                degrees[num] = 0
            degrees[num] += 1
        return degrees

    def get_degree_dists_array(self):
        degrees = self.get_degree_dists()

        arr = n.zeros((8, 2), dtype=n.float64)
        arr[:, 0] = n.arange(2, 10)

        total = 0.
        for degree, number in degrees.items():
            total += number

        for degree, number in degrees.items():
            arr[degree-2, 1] = number / total

        return arr 








            

        

class NeumannDomain(object):
    '''Represents a Neumann domain by storing a list of boundary
    NeumannLines. Also stores the number of saddles, maxima and minima
    participating in the domain.

    Argument dr is the normalisation for distance.

    '''
    def __init__(self, lines, dr=1.0):
        self.lines = lines
        self.dr = dr
        maxima, minima, saddles = 0, 0, 0
        for line in lines:
            start = line.start_type
            if start == 'maximum':
                maxima += 1
            elif start == 'minimum':
                minima += 1
            elif start == 'saddle':
                saddles += 1
        self.maxnum = maxima
        self.minnum = minima
        self.sadnum = saddles

    def diameter(self):
        '''Returns the 'diameter' of the domain, the distance between its
        critical points.'''
        crits = []
        for line in self.lines:
            if line.end_type in ['maximum', 'minimum']:
                crits.append(line.end)
        if len(crits) != 2:
            return None  # Domain detected incorrectly!
        else:
            return mag(n.array(crits[1]) - n.array(crits[0])) * self.dr

    def as_closed_curve(self):
        '''Joins the component lines and returns a single 2d array of points
        making up the domain boundary.

        '''
        points = []
        for line in self.lines:
            points.append(line.points)
        return n.vstack(points)

    def as_sanitised_curves(self):
        '''Not implemented.'''
        pass

    def as_closed_curves(self):
        '''Joins the component lines and returns a list of 2d arrays
        describing sections of the domain cut by any boundary it
        passes through.

        '''
        points = []
        for line in self.lines:
            points.append(line.points)
        arr = n.vstack(points)
        segs = []
        curcut = 0

        for i in range(len(arr)-1):
            next = arr[i+1]
            cur = arr[i]
            if n.abs(next[0]-cur[0]) > 5 or n.abs(next[1]-cur[1]) > 5:
                segs.append(arr[curcut:(i+1)])
                curcut = i+1
        if curcut < len(arr):
            segs.append(arr[curcut:])
        if len(segs)>2:
            lastseg = segs.pop(-1)
            segs[0] = n.vstack((lastseg, segs[0]))
        return segs

    def as_sanitised_curve(self):
        '''Joins the component lines and shifts by the width of the torus if
        the cell crosses a boundary.

        '''
        return sanitise_domain(self.as_closed_curve())

    def number_of_sides(self):
        '''Returns the number of domain walls.'''
        return len(self.lines)

    def vertices(self):
        vertex_set = set()
        for line in self.lines:
            vertex_set.add(line.end)
        return vertex_set

    def area(self):
        return area_from_border(self.as_sanitised_curve()) * self.dr**2

    def perimeter(self):
        points = self.as_sanitised_curve()
        diffs = n.roll(points,-1,axis=0) - points
        return n.sum(n.sqrt(n.sum(diffs*diffs,axis=1))) * self.dr

    def rho(self):
        return self.area() / self.perimeter()

    def guess_type(self):
        '''Tries to work out if self is a star, lens or wedge.'''
        lines = self.lines
        lens_cusps = 0
        star_cusps = 0
        for i, line in enumerate(lines):
            if line.end_type not in ('maximum', 'minimum'):
                continue
            ni = (i+1) % len(lines)
            other_line = lines[ni]
            angle1 = angle_of(line[-1] - line[-2])
            angle2 = angle_of(other_line[1] - other_line[0])
            diff = angle2 - angle1
            if diff > n.pi:
                diff -= 2*n.pi
            if diff < -n.pi:
                diff += 2*n.pi
            if n.abs(diff) < n.pi/2.:
                lens_cusps += 1
            else:
                star_cusps += 1
        if lens_cusps + star_cusps != 2:
            return None
        if lens_cusps == 2:
            return 'lens'
        if star_cusps == 2:
            return 'star'
        if lens_cusps == 1 and star_cusps == 1:
            return 'wedge'
        return None
            
            
            
                

    def crude_area(self):
        '''Crude and slow...do not use!'''
        return crude_area_from_border([line.points for line in self.lines])

    def __str__(self):
        return ('<Neumann domain with {0} saddles,'
                '{1} maxima, {2} minima>').format(self.sadnum,
                                                  self.maxnum,
                                                  self.minnum)

    def __repr__(self):
        return self.__str__()


class NeumannLine(object):
    '''
    Represents a Neumann line. Points may be indexed using array notation.

    Args:
    - start: Coordinates of start node.
    - end: Coordinates of end node.
    - start_type: Start critical point type (maximum, minimum, saddle).
    - end_type: End critical point type.
    - points: 2d array of points making up the line.
    '''
    def __init__(self, start, end, start_type, end_type, points):
        self.start = start
        self.end = end
        self.start_type = start_type
        self.end_type = end_type
        self.points = points

    def inverse(self):
        '''Returns the equivalent line going in the opposite direction and
        with start/end reversed.

        '''
        inv = NeumannLine(self.end, self.start,
                          self.end_type, self.start_type,
                          self.points[::-1])
        return inv

    def jumps(self):
        '''Returns a tuple (x, y), each True/False depending on whether the
        line passes through a periodic boundary.

        '''
        x = False
        y = False
        for i in range(len(self)-1):
            cur = self[i]
            nex = self[i+1]
            if n.abs(nex[0]-cur[0]) > 5:
                x = True
            if n.abs(nex[1]-cur[1]) > 5:
                y = True
        return x, y

    def __str__(self):
        return '<Neumann line: {0} at {1} -> {2} at {3}>'.format(
            self.start_type, self.start, self.end_type, self.end)

    def __repr__(self):
        return self.__str__()

    def __getitem__(self,*args):
        return self.points.__getitem__(*args)

    def __setitem__(self,*args):
        return self.points.__setitem__(*args)

    def __len__(self,*args):
        return len(self.points)


class NeumannTracer(object):
    '''Stores information about a function and appropriate x/y
    information to store information about it in an array.

    Args:

    * xnum: Number of x pixels to sample
    * ynum: Number of y pixels to sample
    * dx: Step length dx for each pixel
    * dy: Step length dy for each pixel
    * func: A function to look for critical points in
    * start_point: The (x, y) tuple to take as (0, 0) in the function.
    * to_edges: May be False (ignore edges), or 'periodic' (use periodic
      boundary conditions)
    * verbose: Whether to print information about ongoing calculations
    * func_params: This may be used to pass extra information about the function,
                   which may let cython be quicker. Right now supports only a tuple of
                   the form "('rwm', wvs, amplitudes, phases)"
    * eigenvalue: the energy of the current eigenvalue, used to normalise rho

    '''
    def __init__(self, xnum, ynum, dx, dy, func,
                 start_point=(0.00123, 0.00123), to_edges=False,
                 upsample=5, verbose=True, isolate_gradients=20,
                 func_params=(), eigenvalue=1.0,
                 area_constraint=None):
        self.arr = n.zeros((xnum, ynum), dtype=n.float64)
        self.hessian_arr = n.zeros((xnum, ynum), dtype=n.float64)
        self.hessian_angle_arr = n.zeros(
            (xnum, ynum), dtype=n.float64)

        self.xnum = xnum
        self.ynum = ynum
        self.shape = (xnum, ynum)
        self.dx = dx
        self.dy = dy
        self.dr = (dx, dy)
        self.func = func
        self.start_point = n.array(start_point).astype(n.double)  # type for cython
        self.sx = start_point[0]
        self.sy = start_point[1]
        self.upsample = upsample
        self.isolate_gradients = isolate_gradients  # precision for max grad change check
        self.eigenvalue = eigenvalue

        self.func_params = func_params

        self.area_constraint = area_constraint

        self.to_edges = to_edges

        self.verbose = verbose

        self.isolated_saddles = []

        self.arr_filled = False
        self.found_crits = False
        self.traced_lines = False
        self.hessian_filled = False
        self.hessian_angle_filled = False
        self.got_hess_dirs = False
        self.graph_built = False
        self.found_domains = False
        self.upsampled_crits = False
        self.upsampling_canon = False  # True if crits have been replaced
                                       # through upsamplint

        # Containers for upsampled point information
        self.upsample_crits = [(), (), (), ()]
        self.upsample_crits_dict = {}

        self.graph = CriticalGraph(xnum, ynum, self.dx)
        self.domains = []
        self.igraph = None

        self.lines = []
        self.noncanon_lines = []
        self.start_points = []
        self.end_points = []
        self.extra_lines = []

        self.crits = [[], [], [], []]

        self.fill_arr()

        self.figax = (None, None)

    def vprint(self, text='', newline=True, force=False):
        '''Prints conditionally following self.verbose.'''
        if self.verbose or force:
            sys.stdout.write(text)
            if newline:
                sys.stdout.write('\n')
            sys.stdout.flush()

    def contains(self, x, y, border=0.):
        '''Return True if the x, y point is within self (sx -> sx + xnum*dx etc.),
        else False.'''
        sx, sy = self.start_point
        dx, dy = self.dr
        xnum, ynum = self.shape
        if ((sx + border < x < sx + xnum*dx - border) and
            (sy + border < y < sy + ynum*dy - border)):
            return True
        return False

    def func_at_coord(self, x, y):
        '''Return the value of self.func at abstract coordinates x, y
        which will be translated via dx and the correct start point.

        '''
        sx, sy = self.start_point
        dx, dy = self.dr
        return self.func(sx + x*dx, sy + y*dy)

    def fill_arr(self, compiled=True):
        '''
        Sample the function on a (self.xnum, self.ynum) array.

        Result stored in self.arr.
        '''
        self.vprint('Populating function sample array...')
        arr = self.arr
        sx, sy = self.start_point
        dx, dy = self.dr
        xnum, ynum = self.shape
        if compiled and cneu is not None:
            cneu.fill_arr(self.func, sx, sy, dx, dy, xnum, ynum, arr)
        else:
            for x in range(xnum):
                self.vprint('\r\tx = {0} / {1}'.format(x, xnum), False)
                for y in range(ynum):
                    arr[x, y] = self.func(sx + x*dx, sy + y*dy)
        self.vprint()
        self.arr_filled = True

    def find_critical_points(self, compiled=True):
        '''
        Iterate over self.arr walking about each pixel and checking the
        number of sign changes. Bins the result appropriately as a
        maximum, minimum, saddle or regular point, and store in
        self.crits.
        '''
        if not self.arr_filled:
            self.fill_arr()
        self.vprint('Finding critical points...')

        if compiled and cneu is not None:
            maxima, minima, saddles, degenerate = cneu.get_critical_points(
                self.arr, self.to_edges, self.verbose)
        else:
            maxima, minima, saddles, degenerate = get_critical_points(
                self.arr, self.to_edges, self.verbose)
            
        self.crits = (maxima, minima, saddles, degenerate)
        #self.prune_critical_points()
        self.maxima = self.crits[0]
        self.minima = self.crits[1]
        self.saddles = self.crits[2]
        self.degenerate = self.crits[3]

        self.crits_dict = critical_points_to_index_dict(self.crits)

        self.found_crits = True

    def upsample_critical_points(self, span=10, compiled=True):
        '''Go through all the critical points of self, upsampling them (or using an
        existing upsampling region if possible).
        '''
        if not self.found_crits:
            self.find_critical_points()

        if compiled and cneu is not None:
            UpsampleTracerClass = cneu.UpsampleTracer  # This isn't really any
                                                       # faster...not a success!
        else:
            UpsampleTracerClass = NeumannTracer

        maxima, minima, saddles, degenerate = self.crits

        upsample = self.upsample
        if upsample <= 1:
            return

        dx, dy = self.dr
        upsample_tracers = []
        self.upsample_tracers = upsample_tracers

        new_dx = dx / upsample
        new_dy = dy / upsample

        for i, crit in enumerate(maxima + minima + saddles):
            if i % 20 == 0:
                self.vprint('\r\tUpsampling critical point {} / {}'.format(
                    i, len(maxima) + len(minima) + len(saddles)), False)
            mx, my = crit
            rx = self.start_point[0] + mx*dx
            ry = self.start_point[1] + my*dy
            if not any([tracer.contains(rx, ry, max(tracer.dx, tracer.dy)) for
                        tracer in upsample_tracers]):
                new_tracer = UpsampleTracerClass(  # Class defined above
                    span*upsample, span*upsample, new_dx, new_dy, self.func,
                    start_point=(self.sx + (mx - span/2.00523)*dx,
                                 self.sy + (my - span/2.00523)*dy),
                    verbose=False)
                new_tracer.find_critical_points()
                upsample_tracers.append(new_tracer)

        self.vprint()
        self.vprint('Replacing critical points...')
        ups_max, ups_min, ups_sad, ups_deg = [], [], [], []
        for i, tracer in enumerate(upsample_tracers):
            new_max, new_min, new_sad, new_deg = tracer.crits
            new_crits_dict = tracer.crits_dict
            for crit, crit_type in new_crits_dict.items():
                mx, my = crit
                rx = tracer.sx + mx*tracer.dx
                ry = tracer.sy + my*tracer.dy

                for earlier_tracer in upsample_tracers[:i]:
                    if earlier_tracer.contains(rx, ry, max(tracer.dx, tracer.dy)):
                        break
                else:
                    # The critical point isn't in the purview of another tracer,
                    # so let's find where it is and store it

                    coord_x = (rx - self.start_point[0]) / self.dx
                    coord_y = (ry - self.start_point[1]) / self.dy

                    if crit_type == 'maximum':
                        ups_max.append((coord_x, coord_y))
                    elif crit_type == 'minimum':
                        ups_min.append((coord_x, coord_y))
                    elif crit_type == 'saddle':
                        ups_sad.append((coord_x, coord_y))
                    elif crit_type == 'degenerate':
                        ups_deg.append((coord_x, coord_y))

        self.upsample_crits = [ups_max, ups_min, ups_sad, ups_deg]
        self.upsample_crits_dict = critical_points_to_index_dict(
            self.upsample_crits)
        self.upsampled_crits = True
            
    def make_upsampling_canon(self):
        if not self.upsampled_crits:
            self.upsample_critical_points()
        ups_max, ups_min, ups_sad, ups_deg = self.upsample_crits
        self.maxima = ups_max
        self.minima = ups_min
        self.saddles = ups_sad
        self.degenerate = ups_deg
        self.crits = self.upsample_crits
        self.crits_dict =self.upsample_crits_dict
        self.upsampling_canon = True
        self.noncanon_lines = self.lines
        self.lines = []
        self.traced_lines = False

    def get_critical_points(self):
        '''Find the critical points, and return (minima, maxima).'''
        if not self.found_crits:
            self.find_critical_points()
        return (self.minima, self.maxima)

    def get_nearest_neighbours(self):
        '''Return a list of distances from maxima to the nearest minimum.

        It does *not* respect periodic boundary conditions, so the
        distances could be lower, but the effect probably isn't
        significant if there are many domains.

        '''
        minima, maxima = self.get_critical_points()
        minima = n.array(minima)
        maxima = n.array(maxima)
        minima_done = n.ones(len(minima), dtype=bool)
        distances = []

        for i in range(len(maxima)):
            dist = 1000.
            match_index = 0
            for j in range(len(minima)):
                if minima_done[j]:
                    cur_dist = mag(minima[j] - maxima[i])
                    if cur_dist < dist:
                        dist = cur_dist
                        match_index = j
            minima_done[j] = False
            distances.append(dist)
        return distances
    def critical_point_proximity_alert(self):
        pass

    def prune_critical_points(self):
        '''Detects and removes points where two extrema border a single
        saddle - it seems that these are (almost?) always numerical
        mistakes.

        This approach is a little crude, a later replacement might
        involve resampling.

        '''
        self.vprint('Pruning critical points')
        maxima, minima, saddles, degenerate = self.crits

        tupmaxima = [tuple(j) for j in maxima]
        tupminima = [tuple(j) for j in minima]

        realmaxima = n.ones(len(maxima), dtype=bool)
        realminima = n.ones(len(minima), dtype=bool)
        realsaddles = n.ones(len(saddles), dtype=bool)

        for i in range(len(saddles)):
            self.vprint('\r\tChecking saddle {0} / {1}'.format(i,
                                                              len(saddles)),
                       False)
            saddle = saddles[i]
            adj = all_adj_indices.copy()
            adj[:, 0] += saddle[0]
            adj[:, 1] += saddle[1]
            adj[:, 0] = adj[:, 0] % self.xnum
            adj[:, 1] = adj[:, 1] % self.ynum
            adjmax = []
            adjmin = []
            for coords in adj:
                coords = tuple(coords)
                if tuple(coords) in tupmaxima:
                    adjmax.append(coords)
                if tuple(coords) in tupminima:
                    adjmin.append(coords)
            if len(adjmax) > 1:
                heights = map(lambda j: self.arr[j[0], j[1]], adjmax)
                fakemax = n.argmin(heights)
                realmaxima[tupmaxima.index(adjmax[fakemax])] = False
                realsaddles[i] = False
            elif len(adjmin) > 1:
                heights = map(lambda j: self.arr[j[0], j[1]], adjmin)
                fakemin = n.argmax(heights)
                realminima[tupminima.index(adjmin[fakemin])] = False
                realsaddles[i] = False
        self.vprint()
        maxima = n.array(maxima)[realmaxima]
        maxima = [tuple(c) for c in maxima]
        minima = n.array(minima)[realminima]
        minima = [tuple(c) for c in minima]
        saddles = n.array(saddles)[realsaddles]
        saddles = [tuple(c) for c in saddles]
        self.crits = maxima, minima, saddles, degenerate

    def gradient_directions_around_saddle(self, saddle, jump, samples=10):
        '''Returns the 4 directions of maximal gradient from the given saddle
        point, from a sampling at the given jump distance.
        '''
        sx, sy = self.start_point
        cx, cy = saddle
        dx, dy = self.dr
        heights = n.zeros(samples)
        func = self.func
        angles = n.linspace(0, 2*n.pi, samples+1)[:-1]
        for index, angle in enumerate(angles):
            jx = sx + cx*dx + jump * n.cos(angle)
            jy = sy + cy*dy + jump * n.sin(angle)
            height = func(jx, jy)
            heights[index] = height

        heights -= func(sx + cx*dx, sy + cy*dy)

        # check sign changes 4 times
        changes = 0
        for index, height in enumerate(heights):
            cur = heights[index]
            nex = heights[(index+1) % samples]
            if n.sign(cur) != n.sign(nex):
                changes += 1

        if changes != 4:
            return None

        cut = int(0.25*samples)

        imax = n.argmax(heights)
        ibefore = imax - cut
        iafter = imax + cut
        # iothermax = (iafter % samples) + n.argmax(n.hstack(
        #     [heights[(iafter % samples):], heights[(ibefore % samples):]]))
        if iafter >= samples or ibefore < 0:
            iothermax = (iafter % samples) + n.argmax(
                heights[(iafter % samples):(ibefore % samples)])
        else:
            iothermax = (iafter % samples) + n.argmax(n.hstack(
                [heights[iafter:], heights[:ibefore]]))
        iothermax %= samples

        imin = n.argmin(heights)
        ibefore = imin - cut
        iafter = imin + cut
        if iafter >= samples or ibefore < 0:
            iothermin = (iafter % samples) + n.argmin(
                heights[(iafter % samples):(ibefore % samples)])
        else:
            iothermin = (iafter % samples) + n.argmin(n.hstack(
                [heights[iafter:], heights[:ibefore]]))
        iothermin %= samples

        return sorted(((1.0, angles[imax]), (-1.0, angles[imin]),
                      (1.0, angles[iothermax]), (-1.0, angles[iothermin])),
                      key=lambda j: j[1])
                                           
        

    def trace_neumann_lines(self, compiled=True):
        '''For every saddle in self.crits, drop 4 Neumann lines at the points
        of adjacent sign change, each gradient ascending/descending
        appropriately until they hit another critical point or appear
        to have stopped. The resulting lines are stored in self.lines.

        If isolate_gradients > 4, the function samples extra points around
        the saddle to find the optimal places to start tracing the gradient.
        This is also performed (by necessity) if upsampling is canon.
        '''
        if not self.arr_filled:
            self.fill_arr()
        if not self.found_crits:
            self.find_critical_points()

        if compiled and cneu is not None:
            gradient_trace_func = cneu.trace_gradient_line
        else:
            gradient_trace_func = trace_gradient_line

        self.vprint('Tracing Neumann lines...')

        arr = self.arr

        isolate_gradients = self.isolate_gradients
        isolated_saddles = []

        integer_saddles = not self.upsampling_canon

        curs = 0
        for saddle in self.saddles:
            if curs % 50 == 0:
                self.vprint('\r\tCurrent saddle {0} / {1}'.format(
                    curs, len(self.saddles)), False)
            curs += 1

            saddlex, saddley = saddle
            if saddlex % 2 == 0:
                ais = even_adj_indices.copy()
            else:
                ais = odd_adj_indices.copy()

            ais[:, 0] += saddlex
            ais[:, 1] += saddley

            val = arr[saddlex, saddley]
            adjs = n.zeros(6, dtype=n.float64)
            for i in range(6):
                adjs[i] = arr[ais[i][0] % self.xnum, ais[i][1] % self.ynum]
                # Always checks periodicity, as saddles on the
                # boundary are already ignored by get_critical_points

            adjs = adjs - val

            if self.upsampling_canon:
                nearby_distance = 1.5 / self.upsample
            else:
                nearby_distance = 1.

            # Find tracing start points
            tracing_start_points = []
            if isolate_gradients or self.upsampling_canon:
                jump_frac = 1.1
                if nearly_integer(saddlex) and nearly_integer(saddley):
                    # Probably could do a real int check but lets be safe
                    jump = jump_frac * self.dx
                else:
                    jump = jump_frac*self.dx/(
                        float(self.upsample) if self.upsampling_canon else 1.)
                tracing_start_points = self.gradient_directions_around_saddle(
                    (saddlex, saddley), jump,
                    isolate_gradients)
                if tracing_start_points is not None:
                    # tracing_start_points = sorted(tracing_start_points,
                    #                               key=lambda j: j[1])
                    tracing_start_points = [
                        (sign, (saddle[0] + jump_frac*n.cos(angle),
                                saddle[1] + jump_frac*n.sin(angle))) for
                        sign, angle in tracing_start_points]
            if not tracing_start_points:
                tracing_start_points = []
                current_region_angles = []
                for i in range(6):
                    cur_adj = adjs[i]
                    next_adj = adjs[(i+1) % 6]
                    current_region_angles.append(n.arctan2(ais[i][1]-saddley,
                                                           ais[i][0]-saddlex))
                    if n.sign(next_adj) != n.sign(cur_adj):
                        sign = n.sign(cur_adj)
                        tracing_start_points.append([sign, ais[i]])
            else:
                isolated_saddles.append(saddle)

            for coords in tracing_start_points:
                sign, coords = coords
                if sign == 1.0:
                    direction = 'down'
                else:
                    direction = 'up'
                diff = [coords[0]-saddlex, coords[1]-saddley]
                self.extra_lines.append(n.array([[saddlex, saddley],
                                                 [saddlex+2*diff[0],
                                                  saddley+2*diff[1]]]))
                points, endcoord = gradient_trace_func(
                    coords[0] + 0.1*diff[0],
                    coords[1] + 0.1*diff[1],
                    self.dx, self.dy,
                    self.xnum, self.ynum,
                    self.func, self.crits_dict,
                    self.start_point,
                    direction, self.to_edges,
                    func_params=self.func_params,
                    area_constraint=self.area_constraint,
                    integer_saddles=integer_saddles,
                    nearby_distance=nearby_distance)

                # Check here whether the line has gone in totally the
                # wrong direction. Try a simple retest if so.
                if endcoord is not None and len(points) > 5:
                    first_vec = n.array(points[1]) - n.array(points[0])
                    last_vec = n.array(points[-1]) - n.array(points[-2])

                    ang_diff = n.arccos(first_vec.dot(last_vec) /
                                        (mag(first_vec) * mag(last_vec)))
                    ang_diff = n.abs(ang_diff)

                    if ang_diff > 0.6*n.pi:
                        print 'trying again'
                        points, endcoord = gradient_trace_func(
                            coords[0] + 0.3*diff[0],
                            coords[1] + 0.3*diff[1],
                            self.dx, self.dy,
                            self.xnum, self.ynum,
                            self.func, self.crits_dict,
                            self.start_point,
                            direction, self.to_edges,
                            func_params=self.func_params,
                            area_constraint=self.area_constraint,
                            integer_saddles=integer_saddles,
                            nearby_distance=nearby_distance)


                if len(points) > 4:
                    points = points[4:]
                elif len(points) < 4:
                    points = points[-1:]
                points = [saddle] + points

                self.start_points.append(tuple(saddle))
                if endcoord is not None:
                    self.end_points.append(tuple(endcoord))
                else:
                    self.end_points.append(None)
                self.lines.append(n.array(points))

                if (endcoord is not None and
                    tuple(endcoord) not in self.minima and
                    tuple(endcoord) not in self.maxima):
                    print 'adding new endcoord', tuple(endcoord)
                    fval = self.func(self.sx + endcoord[0]*self.dx,
                                     self.sy + endcoord[1]*self.dy)
                    if fval > 0:
                        self.maxima.append(tuple(endcoord))
                    else:
                        self.minima.append(tuple(endcoord))

        self.vprint()
        self.crits_dict = critical_points_to_index_dict(self.crits)
        self.isolated_saddles = isolated_saddles
        self.traced_lines = True

    def print_critical_heights(self):
        '''
        Print the height of every critical point in self.crits.
        '''
        if not self.arr_filled:
            self.fill_arr()
        if not self.found_crits:
            self.find_critical_points()

            
        print 'maxima'
        for entry in self.maxima:
            print (entry, self.func(self.sx+entry[0]*self.dx,
                                    self.sy+entry[1]*self.dy),
                   self.arr[entry[0], entry[1]])
        print 'minima'
        for entry in self.minima:
            print (entry, self.func(self.sx+entry[0]*self.dx,
                                    self.sy+entry[1]*self.dy),
                   self.arr[entry[0], entry[1]])
        print 'saddles'
        for entry in self.saddles:
            print (entry, self.func(self.sx+entry[0]*self.dx,
                                    self.sy+entry[1]*self.dy),
                   self.arr[entry[0], entry[1]])
        print 'degenerate'
        for entry in self.degenerate:
            print (entry, self.func(self.sx+entry[0]*self.dx,
                                    self.sy+entry[1]*self.dy),
                   self.arr[entry[0], entry[1]])

    def calculate_hessian_directions(self):
        '''Associate hessian direction vectors with each critical point.'''
        if not self.found_crits:
            self.find_critical_points()
        maxima, minima, saddles, degenerate = self.crits
        hess_dirs = {}
        for point in maxima + minima + saddles:
            hess = hessian(self.func, point[0], point[1],
                           self.dx, self.dy)
            eigs, eigvs = n.linalg.eig(hess)
            ev1, ev2 = eigvs
            ev1 *= 0.01*self.xnum / mag(ev1) 
            ev2 *= 0.01*self.ynum / mag(ev2) 
            vpoint = n.array(point)
            v1 = n.array([vpoint - 2*ev1, vpoint + 2*ev1])
            v2 = n.array([vpoint - 2*ev2, vpoint + 2*ev2])
            hess_dirs[point] = (v1, v2)

        self.hess_dirs = hess_dirs
        self.got_hess_dirs = True
        return hess_dirs

    def make_hessian_array(self, compiled=True):
        '''
        Calculate the hessian at every point of self.arr.
        '''
        self.vprint('Filling Hessian domain array...')
        arr = self.hessian_arr
        sx, sy = self.start_point
        dx, dy = self.dr
        xnum, ynum = self.shape
        if compiled and cneu is not None:
            hessian_det_func = cneu.hessian_det
        else:
            hessian_det_func = hessian_det
        for x in range(xnum):
            self.vprint('\r\tx = {0} / {1}'.format(x, xnum), False)
            for y in range(ynum):
                arr[x, y] = hessian_det_func(self.func, sx + x*dx,
                                             sy + y*dy, dx, dy)

        self.hessian_filled = True
        self.vprint()

    def make_hessian_angle_array(self, compiled=True):
        '''Calculate the angle of the Hessian major axis
        at every point of self.arr.'''
        self.vprint('Filling Hessian angle array...')
        arr = self.hessian_angle_arr
        sx, sy = self.start_point
        dx, dy = self.dr
        xnum, ynum = self.shape
        if compiled and cneu is not None:
            hessian_func = cneu.hessian
        else:
            hessian_func = hessian
        for x in range(xnum):
            self.vprint('\r\tx = {0} / {1}'.format(x, xnum), False)
            for y in range(ynum):
                hessian = hessian_func(self.func, sx + x*dx,
                                       sy + y*dy, dx, dy)
                eigs, eigvs = n.linalg.eig(hessian)
                vec = eigvs[n.argmin(eigs)]
                angle = n.arctan2(vec[1], vec[0])
                if angle < 0:
                    angle += n.pi
                arr[x, y] = angle

        self.hessian_angle_filled = True
        self.vprint()

    def build_graph(self):
        '''Build a CriticalGraph with all the critical points of self as
        nodes and all the lines of self.lines joining them.

        '''
        if not self.arr_filled:
            self.fill_arr()
        if not self.found_crits:
            self.find_critical_points()
        if not self.traced_lines:
            self.trace_neumann_lines()

        for i in range(len(self.lines)):
            start = self.start_points[i]
            end  = self.end_points[i]
            line = self.lines[i]
            start_crit_type = self.crits_dict[start]
            if end is not None:
                # print end
                # if end in self.maxima:
                #     print 'maximum'
                # elif end in self.minima:
                #     print 'minimum'
                # elif end in self.saddles:
                #     print 'saddle'
                # else:
                #     print 'none'
                #print self.crits_dict
                end_crit_type = self.crits_dict[end]
            else:
                end_crit_type = None
            neuline = NeumannLine(start, end, start_crit_type,
                                  end_crit_type, line)
            self.graph.add_or_edit_node(start, start_crit_type, neuline)
            self.graph.add_or_edit_node(end, end_crit_type, neuline.inverse())

        self.graph_built = True

    def get_igraph(self):
        if not self.graph_built:
            self.build_graph()
        crit_graph = self.graph

        nodes = crit_graph.keys()
        numbers_dict = {}
        for node, number in zip(nodes, range(len(nodes))):
            numbers_dict[node] = number

        g = ig.Graph()
        g.add_vertices(len(nodes))
        for node in crit_graph:
            crit_type, lines = crit_graph[node]
            if crit_type == 'saddle':
                for line in lines:
                    end = line.end
                    g.add_edge(numbers_dict[node], numbers_dict[end])
        self.igraph = g
        return g

    def get_igraph_dual(self):
        if not self.graph_built:
            self.build_graph()
        dual = self.graph.dual_graph()

        numbers_dict = {}
        nodes = dual.keys()
        for node, number in zip(nodes, range(len(nodes))):
            numbers_dict[node] = number

        g = ig.Graph()
        g.add_vertices(len(dual))
        for node, adjacent in dual.iteritems():
            for other_node in adjacent:
                g.add_edge(numbers_dict[node], numbers_dict[other_node])
        self.dual_graph = g
        return g

    def get_recognised_domains(self):
        '''Get the list of all closed Neumann domains in self.graph (created
        by self.build_graph)

        '''
        if not self.arr_filled:
            self.fill_arr()
        if not self.found_crits:
            self.find_critical_points()
        if not self.traced_lines:
            self.trace_neumann_lines()
        if not self.graph_built:
            self.build_graph()
        if self.found_domains:
            return self.domains
        self.found_domains = True
        domains = self.graph.get_closed_domains()
        self.domains = domains

        return domains

    def get_critical_degrees(self):
        ds = self.get_recognised_domains()

        maxima = self.maxima
        minima = self.minima

        maxdegrees = {}
        mindegrees = {}
        for node in maxima:
            if node in self.graph:
                lines = self.graph[node][1]
                maxdegrees[node] = len(lines)
            else:
                maxdegrees[node] = 0
        for node in minima:
            if node in self.graph:
                lines = self.graph[node][1]
                mindegrees[node] = len(lines)
            else:
                mindegrees[node] = 0

        maxdomaindegrees = {key: 0 for key in maxima}
        mindomaindegrees = {key: 0 for key in minima}
        for domain in ds:
            lines = domain.lines
            for line in lines:
                end = line.end
                if end in maxdomaindegrees:
                    maxdomaindegrees[end] += 1
                if end in mindomaindegrees:
                    mindomaindegrees[end] += 1

        maxrealdegrees = []
        minrealdegrees = []
        realmaxima = []
        realminima = []
        for node in maxima:
            #print maxdegrees[node], maxdomaindegrees[node]
            line_degree = maxdegrees[node]
            domain_degree = maxdomaindegrees[node]
            if line_degree == domain_degree:
                maxrealdegrees.append(line_degree)
                realmaxima.append(node)
        for node in minima:
            line_degree = mindegrees[node]
            domain_degree = mindomaindegrees[node]
            if line_degree == domain_degree:
                minrealdegrees.append(line_degree)
                realminima.append(node)
            #print mindegrees[node], mindomaindegrees[node]
        return maxrealdegrees, minrealdegrees

    def get_critical_degree_dists(self):
        maxdegs, mindegs = self.get_critical_degrees()
        maxdegs = n.array(maxdegs, dtype=n.float64)
        mindegs = n.array(mindegs, dtype=n.float64)

        max2 = n.sum(maxdegs == 2)
        max3 = n.sum(maxdegs == 3)
        max4 = n.sum(maxdegs == 4)
        max5 = n.sum(maxdegs == 5)
        max6 = n.sum(maxdegs == 6)
        max7 = n.sum(maxdegs == 7)
        max8 = n.sum(maxdegs == 8)
        max9 = n.sum(maxdegs == 9)
        min2 = n.sum(mindegs == 2)
        min3 = n.sum(mindegs == 3)
        min4 = n.sum(mindegs == 4)
        min5 = n.sum(mindegs == 5)
        min6 = n.sum(mindegs == 6)
        min7 = n.sum(mindegs == 7)
        min8 = n.sum(mindegs == 8)
        min9 = n.sum(mindegs == 9)

        totdegrees = float(max2 + max3 + max4 + max5 + max6 + max7 + max8 +
                           min2 + min3 + min4 + min5 + min6 + min7 + min8)

        if totdegrees == 0:
            return n.array(zip(range(2, 10), n.ones(8)))
        frac2 = float(max2 + min2) / totdegrees
        frac3 = float(max3 + min3) / totdegrees
        frac4 = float(max4 + min4) / totdegrees
        frac5 = float(max5 + min5) / totdegrees
        frac6 = float(max6 + min6) / totdegrees
        frac7 = float(max7 + min7) / totdegrees
        frac8 = float(max8 + min8) / totdegrees
        frac9 = float(max9 + min9) / totdegrees

        return n.array(zip(range(2, 10),
                           (frac2, frac3, frac4, frac5,
                            frac6, frac7, frac8, frac9)))

    def get_domain_areas(self):
        if not self.graph_built:
            self.build_graph()
        return self.graph.get_domain_areas()

    def get_domain_perimeters(self):
        if not self.graph_built:
            self.build_graph()
        return self.graph.get_domain_perimeters()

    def get_domain_diameters(self):
        if not self.graph_built:
            self.build_graph()
        return self.graph.get_domain_diameters()

    def get_domain_rhos(self):
        if not self.graph_built:
            self.build_graph()
        return self.graph.get_domain_rhos(self.eigenvalue)

    def get_domain_rhos_by_type(self):
        if not self.graph_built:
            self.build_graph()
        return self.graph.get_domain_rhos_by_type(self.eigenvalue)

    def get_neumann_heights(self):
        '''Returns a list of heights of all points on detected
        Neumann lines.'''
        if not self.traced_lines:
            self.trace_neumann_lines()
        lines = self.lines
        heights = []
        for line in lines:
            for point in line:
                height = self.func(*point)
                heights.append(height)

        return heights

    def get_neumann_angles(self):
        '''Returns a list of angles of the tangent at every Neumann
        line point.'''
        if not self.traced_lines:
            self.trace_neumann_lines()
        lines = self.lines
        angles = []
        for line in lines:
            ts = n.roll(line, -1, axis=0) - line
            for x, y in ts:
                angle = n.arctan2(y, x)
                angles.append(angle)
        return angles

    def get_extrema_heights(self, mod=True):
        '''Returns a list of heights at extrema, optionally normed
        (taking the abs value).'''
        maxima_heights = [self.func(*maximum) for maximum in self.maxima]
        minima_heights = [self.func(*minimum) for minimum in self.minima]
        if mod:
            maxima_heights = [n.abs(height) for height in maxima_heights]
            minima_heights = [n.abs(height) for height in minima_heights]
        return maxima_heights + minima_heights

    def get_saddle_heights(self, mod=True):
        '''Returns a list of heights at saddles, optionally normed
        (taking the abs value).'''
        saddle_heights = [self.func(*saddle) for saddle in self.saddles]
        if mod:
            saddle_heights = [n.abs(height) for height in saddle_heights]
        return saddle_heights

    def get_critical_heights(self, mod=True):
        '''Returns a list of heiths at all extrema.'''
        return self.get_extrema_heights(mod) + self.get_saddle_heights(mod)

    def get_neumann_gradients(self, along_line=True, mod=True):
        '''Returns (very approximate) height gradients along the lines.'''
        if not self.traced_lines:
            self.trace_neumann_lines()
        lines = self.lines

    def build_everything(self, including_hessian=False):
        '''
        Build all the arrays and trace all the lines via the various
        methods of self.

        Args:
        - including_hessian: Default False, whether to bother doing the
          Hessian array (this can be quite slow to make)
        '''
        if not self.arr_filled:
            self.fill_arr()
        if not self.found_crits:
            self.find_critical_points()
        if not self.traced_lines:
            self.trace_neumann_lines()
        if not self.graph_built:
            self.build_graph()
        if not self.found_domains:
            self.get_recognised_domains()
        if not self.hessian_filled and including_hessian:
            self.make_hessian_array()

    def get_gradients_array(self):
        self.vprint('Calculating gradients array...')
        grad_arr = self.arr.copy()
        sx, sy = self.start_point
        dx, dy = self.dr
        xnum, ynum = self.shape
        for x in range(xnum):
            self.vprint('\r\tx = {0} / {1}'.format(x, xnum), False)
            for y in range(ynum):
                gradient = grad(self.func, sx+x*dx, sy+y*dy, dx, dy)
                grad_arr[x, y] = n.arctan2(gradient[1], gradient[0])
        self.vprint()
        return grad_arr

    def get_nearby_gradients_array(self):
        self.vprint('Calculating nearby gradients array...')
        grad_arr = self.arr.copy()
        sx, sy = self.start_point
        dx, dy = self.dr
        xnum, ynum = self.shape
        for x in range(xnum):
            self.vprint('\r\tx = {0} / {1}'.format(x, xnum), False)
            for y in range(ynum):
                gradient = grad(self.func, sx+x*dx, sy+y*dy, dx, dy)
                orthdir = rotation_matrix(n.pi/2).dot(gradient)
                orthdir /= mag(orthdir)
                orthdir /= 20.

                orth_grad_1 = grad(self.func, sx+x*dx + orthdir[0],
                                   sy+y*dy + orthdir[1],
                                   dx, dy)
                orth_grad_2 = grad(self.func, sx+x*dx - orthdir[0],
                                   sy+y*dy - orthdir[1],
                                   dx, dy)

                angle1 = n.arctan2(orth_grad_1[1], orth_grad_1[0])
                anglegrad = n.arctan2(gradient[1], gradient[0])
                angle2 = n.arctan2(orth_grad_2[1], orth_grad_2[0])

                corr_1 = n.abs(anglegrad - angle1)
                if corr_1 > n.pi:
                    corr_1 -= n.pi
                corr_2 = n.abs(anglegrad - angle2)
                if corr_2 > n.pi:
                    corr_2 -= n.pi

                # rel_ang_1 = angle_between_vectors(orth_grad_1, gradient)
                # rel_ang_2 = angle_between_vectors(orth_grad_2, gradient)

                grad_arr[x, y] = n.exp(-1 * (corr_1 + corr_2)**2 / (n.pi/4)**2)
        self.vprint()
        return grad_arr

    def get_correlation_functions(self):
        if not self.found_crits:
            self.find_critical_points()
        maxima, minima, saddles, degenerate = self.crits
        xnum, ynum = self.shape

        extremum_extremum = []
        extremum_saddle = []
        crit_crit = []
        max_min = []
        saddle_saddle = []

        for index, maximum in enumerate(maxima):
            self.vprint('\r\tmax index {} / {}'.format(index,
                                                       len(maxima)), False)
            pos = n.array(maximum)
            for other_index, other_maximum in enumerate(
                    maxima[(index+1):]):
                other_pos = n.array(other_maximum)
                distance = mag(other_pos - pos)
                if distance > 0.5*self.xnum:
                    distance = reduce_distance(maximum, other_maximum,
                                               xnum, ynum)
                extremum_extremum.append(distance)
                crit_crit.append(distance)
            for other_index, minimum in enumerate(minima):
                other_pos = n.array(minimum)
                distance = mag(other_pos - pos)
                if distance > 0.5*self.xnum:
                    distance = reduce_distance(maximum, minimum,
                                               xnum, ynum)
                max_min.append(distance)
                extremum_extremum.append(distance)
                crit_crit.append(distance)
                max_min.append(distance)
            for other_index, saddle in enumerate(saddles):
                other_pos = n.array(saddle)
                distance = mag(other_pos - pos)
                if distance > 0.5*self.xnum:
                    distance = reduce_distance(maximum, saddle,
                                               xnum, ynum)
                extremum_saddle.append(distance)
                crit_crit.append(distance)
                

        for index, minimum in enumerate(minima):
            self.vprint('\r\tmin index {} / {}'.format(index,
                                                       len(minima)), False)
            pos = n.array(minimum)
            for other_index, other_minimum in enumerate(
                    maxima[(index+1):]):
                other_pos = n.array(other_minimum)
                distance = mag(other_pos - pos)
                if distance > 0.5*self.xnum:
                    distance = reduce_distance(maximum, other_minimum,
                                               xnum, ynum)
                extremum_extremum.append(distance)
                crit_crit.append(distance)
            for other_index, saddle in enumerate(saddles):
                other_pos = n.array(saddle)
                distance = mag(other_pos - pos)
                if distance > 0.5*self.xnum:
                    distance = reduce_distance(minimum, saddle,
                                               xnum, ynum)
                extremum_saddle.append(distance)
                crit_crit.append(distance)

        for index, saddle in enumerate(saddles):
            self.vprint('\r\tsaddle index {} / {}'.format(index,
                                                       len(saddles)), False)
            pos = n.array(saddle)
            for other_index, saddle in enumerate(saddles):
                other_pos = n.array(saddle)
                distance = mag(other_pos - pos)
                if distance > 0.5*self.xnum:
                    distance = reduce_distance(minimum, saddle,
                                               xnum, ynum)
                saddle_saddle.append(distance)
        return (extremum_extremum, extremum_saddle, crit_crit,
                max_min, saddle_saddle)
                
            

    def get_voronoi_diagram(self):
        if not self.found_crits:
            self.find_critical_points()

        maxima, minima, saddles, degenerate = self.crits
        points = maxima + minima + saddles
        v = Voronoi(n.array(points))
        return v

    def get_delaunay_diagram(self):
        if not self.found_crits:
            self.find_critical_points()

        maxima, minima, saddles, degenerate = self.crits
        points = maxima + minima + saddles
        v = Delaunay(n.array(points))
        return v

    def get_delaunay_graph(self):
        if not self.found_crits:
            self.find_critical_points()

        maxima, minima, saddles, degenerate = self.crits
        points = maxima + minima + saddles
        return DelaunayGraph(points, self.crits_dict)

    def plot(self, show_critical_points=True,
             trace_lines=True, plot_hessian=False,
             plot_hessian_angles=False,
             show_saddle_directions=False,
             show_domain_patches=False,
             constrain_patch_rhos=None,
             colour_patches_by_rho=False,
             colour_patches_by_domain_type=False,
             print_patch_areas=False,
             print_patch_rhos=False,
             figsize=None,
             show_sample_directions=False,
             show_hessian_directions=False,
             plot_gradients=False,
             plot_nearby_gradients=False,
             plot_delaunay=False,
             plot_voronoi=False,
             plot_reduced_delaunay=False,
             highlight_isolated_saddles=False,
             show_noncanon_lines=True,
             save=False, figax=None,
             maxima_style=maxima_style_old,  # Defined near top of file
             minima_style=minima_style_old,
             saddle_style = saddle_style_rami,
             interpolation='none',
             symmetrise_colours=False,
             cmap='RdYlBu_r'):
        '''
        Plot and return a graph showing (optionally):
        - Neumann lines
        - Hessian domains
        - Hessian eigenvectors at detected saddles
        - Coloured patches representing each closed Neumann domain

        Options include:
        - figsize: set the size of the output figure in inches
        - show_sample_directions: plot (in red) markers showing the
                                  directions in which Neumann lines were
                                  dropped
        - plot_gradients: plot (in hsv) the gradient at each point
        - plot_nearby_gradients: same, but nearby?
        - save: save at the given filename
        - figax: plot to the given (fig, ax) tuple, if any
        - interpolation: the interpolation of the main image array
                         (try 'bicubic' to blend)
        - cmap: the name of the colormap, defaults to 'Spectral_r', used
                to be 'RdBuYl'.
        '''
        if save:
            mpl_interactive(False)
        else:
            mpl_interactive(True)

        if not self.arr_filled:
            self.fill_arr()
        if not self.found_crits and show_critical_points:
            self.find_critical_points()
        if not self.traced_lines and trace_lines:
            self.trace_neumann_lines()
        if not self.hessian_filled and plot_hessian:
            self.make_hessian_array()
        if not self.found_domains and show_domain_patches:
            self.get_recognised_domains()
        if not self.got_hess_dirs and show_hessian_directions:
            self.calculate_hessian_directions()
        if not self.hessian_angle_filled and plot_hessian_angles:
            self.make_hessian_angle_array()

        plotarr = n.rot90(self.arr[::-1], 3).copy()
        if symmetrise_colours:
            absmax = n.max((n.max(plotarr), n.abs(n.min(plotarr))))
            plotarr[0, 0] = absmax
            plotarr[0, 1] = -1*absmax

        #fig, ax = plot_arr_with_crits(plotarr, self.crits)
        maxima, minima, saddles, degenerate = self.crits
        maxima = n.array(maxima)
        minima = n.array(minima)
        saddles = n.array(saddles)
        degenerate = n.array(degenerate)

        # Upsampled critical points (if they exist)
        (upsampled_maxima, upsampled_minima, upsampled_saddles,
         upsampled_degenerate) = self.upsample_crits
        upsampled_maxima = n.array(upsampled_maxima)
        upsampled_minima = n.array(upsampled_minima)
        upsampled_saddles = n.array(upsampled_saddles)
        upsampled_degenerate = n.array(upsampled_degenerate)

        fig, ax = plt.subplots()
        # if figax is None:
        #     if self.figax[0] is not None and self.figax[1] is not None:
        #         fig, ax = self.figax
        #         ax.clear()
        #     else:
        #         fig, ax = plt.subplots()
        # else:
        #     fig, ax = figax
        #     ax.clear()
        # fig.show()

        if not show_domain_patches:
            ax.imshow(plotarr, cmap=cmap, interpolation='none',
                      alpha=0.6)
        else:
            ax.imshow(plotarr, cmap=cmap, interpolation='none',
                      alpha=0.4)

        if plot_gradients:
            grad_arr = self.get_gradients_array()
            grad_arr = n.rot90(grad_arr[::-1], 3)
            ax.imshow(grad_arr, cmap='hsv', interpolation='none')
        if plot_nearby_gradients:
            grad_arr = self.get_nearby_gradients_array()
            grad_arr = n.rot90(grad_arr[::-1], 3)
            ax.imshow(grad_arr, cmap='pink', interpolation='none')

        ax.set_xlim(0, plotarr.shape[0])
        ax.set_ylim(0, plotarr.shape[1])

        ax.set_xticks([])
        ax.set_yticks([])

        if show_domain_patches:
            for domain in self.domains:

                if colour_patches_by_rho:
                    cfunc = (hsv if colour_patches_by_rho == 'hsv'
                            else jet)
                    rho = domain.rho() * n.sqrt(self.eigenvalue)
                    if colour_patches_by_rho == 'hsv':
                        colour = cfunc((rho-0.3) / 0.6)[:3]
                    else:
                        colour = cfunc(rho / 1.3)[:3]
                    alpha = 0.95
                elif colour_patches_by_domain_type:
                    domain_type = domain.guess_type()
                    if domain_type == 'star':
                        colour = (0., 0., 1.)
                    elif domain_type == 'lens':
                        colour = (1., 0., 0.)
                    elif domain_type == 'wedge':
                        colour = (0., 1., 0.)
                    else:
                        colour = (0.3, 0.3, 0.3)
                    alpha = 0.7
                else:
                    colour = hsv_to_rgb(n.random.random(), 1., 1.)
                    alpha = 0.7
                ps = domain.as_closed_curves()
                cpr = constrain_patch_rhos
                if (cpr is not None and not
                    cpr[0] < domain.rho() * n.sqrt(self.eigenvalue) < cpr[1]):
                    continue
                for p in ps:
                    if len(ps) > 1:
                        patch = Polygon(p, alpha=alpha,
                                        closed=True,
                                        color=colour,
                                        # linestyle=n.random.choice(
                                        #     patch_linestyles),
                                        linewidth=2,
                                        )
                    else:
                        patch = Polygon(p, alpha=alpha,
                                        closed=True,
                                        color=colour,
                                        # linestyle=n.random.choice(
                                        #     patch_linestyles),
                                        linewidth=2,
                                        )
                    ax.add_patch(patch)
                if print_patch_areas:
                    area = domain.area()
                    patch_areas = map(area_from_border, ps)
                    pos = n.average(ps[n.argmax(patch_areas)], axis=0)
                    ax.text(pos[0], pos[1],'{:.1f}'.format(area))
                if print_patch_rhos:
                    rho = domain.rho() * n.sqrt(self.eigenvalue)
                    patch_areas = map(area_from_border, ps)
                    pos = n.average(ps[n.argmax(patch_areas)], axis=0)
                    ax.text(pos[0], pos[1],'{:.3f}'.format(rho))

        legend_entries = []
        if self.upsample_crits_dict and not self.upsampling_canon:
            initial_crit_alpha=0.3
        else:
            initial_crit_alpha=1.0
        if show_critical_points:
            if len(maxima) > 0:
                ax.scatter(maxima[:, 0], maxima[:, 1], 60,
                           alpha=initial_crit_alpha,
                           **maxima_style)
                legend_entries.append('maxima')
            if len(minima) > 0:
                ax.scatter(minima[:, 0], minima[:, 1], 60,
                           alpha=initial_crit_alpha,
                           **minima_style)
                legend_entries.append('minima')
            if len(saddles) > 0:
                ax.scatter(saddles[:, 0], saddles[:, 1], 60, 
                           alpha=initial_crit_alpha,
                           **saddle_style)
                legend_entries.append('saddles')
            if len(degenerate) > 0:
                ax.scatter(degenerate[:, 0], degenerate[:, 1], c='orange',
                           alpha=initial_crit_alpha)
                legend_entries.append('degenerate')

            # And the same for upsampled if they exist
            if not self.upsampling_canon:
                if len(upsampled_maxima) > 0:
                    print upsampled_maxima[:, 0], upsampled_maxima[:, 1]
                    ax.scatter(upsampled_maxima[:, 0], upsampled_maxima[:, 1],
                               60, c='r')
                    legend_entries.append('upsampled maxima')
                if len(upsampled_minima) > 0:
                    ax.scatter(upsampled_minima[:, 0], upsampled_minima[:, 1],
                               60, c='b')
                    legend_entries.append('upsampled minima')
                if len(upsampled_saddles) > 0:
                    ax.scatter(upsampled_saddles[:, 0], upsampled_saddles[:, 1],
                               60, color='yellow')
                    legend_entries.append('saddles')
                if len(upsampled_degenerate) > 0:
                    ax.scatter(upsampled_degenerate[:, 0],upsampled_degenerate[:, 1],
                               c='orange')
                    legend_entries.append('upsampled degenerate')

        if self.isolated_saddles and highlight_isolated_saddles:
            isols = n.array(self.isolated_saddles)
            ax.scatter(isols[:, 0], isols[:, 1], 200, color='black', alpha=0.2)
                

        if show_domain_patches:
            ax.contour(plotarr, levels=[0], alpha=0.3, linestyles=['--'])
        else:
            ax.contour(plotarr, levels=[0], alpha=0.2)

        if show_noncanon_lines:
            for line in self.noncanon_lines:
                segs = sanitise_line(line)
                for seg in segs:
                    ax.plot(seg[:, 0], seg[:, 1],'-', color='green', alpha=0.5)
        if trace_lines:
            for line in self.lines:
                segs = sanitise_line(line)
                for seg in segs:
                    ax.plot(seg[:, 0], seg[:, 1],'-', color='purple')


        if show_sample_directions:
            for line in self.extra_lines:
                if (n.abs(line[0, 0]-line[1, 0] < 5) and
                    n.abs(line[0, 1] - line[1, 1] < 5)):
                    ax.plot(line[:, 0], line[:, 1],'-', color='red')

        if plot_hessian:
            hessian_arr = n.rot90(self.hessian_arr[::-1], 3)
            ax.imshow(n.sign(hessian_arr), cmap='binary',
                      interpolation='none', alpha=0.5)
            ax.contour(hessian_arr, levels=[0],
                       linewidths=2, alpha=0.6, color='cyan')

        if plot_hessian_angles:
            angle_arr = n.rot90(self.hessian_angle_arr[::-1], 3)
            ax.imshow(angle_arr, cmap='hsv',
                      interpolation='none', alpha=0.75)

        if show_hessian_directions:
            hess_dirs = self.hess_dirs
            for point, vs in hess_dirs.items():
                v1, v2 = vs
                ax.plot(v1[:, 0], v1[:, 1],
                        color=(0, 0.5, 0),
                        linewidth=1.5)
                ax.plot(v2[:, 0], v2[:, 1],
                        color=(0, 0.5, 0),
                        linewidth=1.5)

        self.figax = (fig, ax)

        if plot_delaunay:
            d = self.get_delaunay_diagram()
            from matplotlib import tri
            tri.triplot(ax, d.points[:, 0], d.points[:, 1], color='green',
                        linewidth=1.5)

        if plot_voronoi:
            print 'WARNING: Voronoi plot doesn\'t, will just plot Delaunay.'''
            # Fix this!
            d = self.get_voronoi_diagram()
            from matplotlib import tri
            tri.triplot(ax, d.points[:, 0], d.points[:, 1], color='red',
                        linewidth=1.5)

        if plot_reduced_delaunay:
            d = self.get_delaunay_diagram()
            points = d.points
            vertices = d.simplices
            triangles = points[vertices]
            cd = self.crits_dict
            for triangle in triangles:
                triangle = triangle.astype(n.int64)
                types = [cd[tuple(triangle[i])] for i in range(3)]
                for i in range(3):
                    current_points = n.roll(triangle, i, axis=0)[:2]
                    current_types = n.roll(types, i)[:2]
                    if ((current_types[0] == 'saddle' and
                         current_types[1] == 'saddle') or
                        ('maximum' in current_types and
                         'minimum' in current_types)):
                        pass
                    else:
                        print [tuple(row) for row in current_points]
                        ax.plot(current_points[:, 0], current_points[:, 1],
                                color='green')


        if show_saddle_directions:
            saddles = self.saddles
            sx, sy = self.start_point
            dx, dy = self.dr
            for saddle in saddles:
                x, y = saddle
                hess = hessian(self.func, sx+x*dx, sy+y*dy, dx, dy)
                eigs, eigvs = n.linalg.eig(hess)
                dir1 = eigvs[0] / mag(eigvs[0])
                dir2 = eigvs[1] / mag(eigvs[1])
                xs1 = [x-3*dir1[0], x, x+3*dir1[0]]
                ys1 = [y-3*dir1[1], y, y+3*dir1[1]]
                xs2 = [x-3*dir2[0], x, x+3*dir2[0]]
                ys2 = [y-3*dir2[1], y, y+3*dir2[1]]
                ax.plot(xs1, ys1, color='black', linewidth=2)
                ax.plot(xs2, ys2, color='black', linewidth=2)

        if figsize is not None:
            fig.set_size_inches(figsize[0], figsize[1])

        if save:
            filen = save
            fig.savefig(filen)

        return fig, ax

    def plot_bokeh(self, trace_lines=True, plot_hessian=False,
                   show_saddle_directions=False,
                   show_domain_patches=False,
                   print_patch_areas=False,
                   figsize=None,
                   show_sample_directions=False,
                   save=False, figax=None):
        '''Duplicate plot method using bokeh rather than matplotlib.
        
        Plot and return a graph showing (optionally):
        - Neumann lines
        - Hessian domains
        - Hessian eigenvectors at detected saddles
        - Coloured patches representing each closed Neumann domain

        This is a WIP (and doesn't really work right now)
        '''

        try:
            import bokeh.plotting as plotting
        except ImportError:
            print 'Failed to import bokeh. Cancelling plot.'
            return

        plotting.output_file("bokehtest.html", title="bokeh test?")
        
        if not self.arr_filled:
            self.fill_arr()
        if not self.found_crits:
            self.find_critical_points()
        if not self.traced_lines and trace_lines:
            self.trace_neumann_lines()
        if not self.hessian_filled and plot_hessian:
            self.make_hessian_array()
        if not self.found_domains and show_domain_patches:
            self.get_recognised_domains()

        plotarr = n.rot90(self.arr[::-1], 3)

        #fig, ax = plot_arr_with_crits(plotarr, self.crits)
        maxima, minima, saddles, degenerate = self.crits
        maxima = n.array(maxima)
        minima = n.array(minima)
        saddles = n.array(saddles)
        degenerate = n.array(degenerate)

        # Show heightmap image
        # ax.imshow(plotarr, cmap='Spectral_r', interpolation='none', alpha=0.6)

        # Set axis scales
        # Remove ticks

        if show_domain_patches:
            patch_colours = []
            patch_xs = []
            patch_ys = []
            for domain in self.domains:
                colour = hsv_to_rgb(n.random.random(), 1., 1.)
                ps = domain.as_closed_curves()
                for p in ps:
                    patch_xs.append(p[:,0])
                    patch_ys.append(p[:,1])
                    patch_colours.append(colour)
            plotting.patches(patch_xs, patch_ys, fill_color=patch_colours,
                             fill_alpha=0.7, line_color='black',
                             line_width=0.5)
                # Print patch areas

        legend_entries = []
        if len(maxima) > 0:
            plotting.scatter(maxima[:,0], maxima[:,1], name='maxima',
                             color='#ff0000')
            legend_entries.append('maxima')
        if len(minima) > 0:
            plotting.scatter(minima[:,0], minima[:,1], name='minima',
                             color='#0000ff')
            legend_entries.append('minima')
        if len(saddles) > 0:
            plotting.scatter(saddles[:,0], saddles[:,1], name='saddles',
                             color='#ffff00')
            legend_entries.append('saddles')
        if len(degenerate) > 0:
            plotting.scatter(degenerate[:,0], degenerate[:,1],
                             name='degenerate', color='#ffffff')
            legend_entries.append('degenerate')

        # Contours of plotarr

        if trace_lines:
            for line in self.lines:
                segs = sanitise_line(line)
                for seg in segs: 
                    plotting.line(seg[:, 0], seg[:,1], color='#DD00DD')

        # Plot hessian imarr
            # hessian_arr = n.rot90(self.hessian_arr[::-1], 3)
            # ax.imshow(n.sign(hessian_arr), cmap='binary',
            #           interpolation='none', alpha=0.5)
            # ax.contour(hessian_arr, levels=[0],
            #            linewidths=2, alpha=0.6, color='cyan')

        plotting.show()

    def plot_torus(self, clf=True, cmap='RdYlBu',
                   trace_lines=True,
                   amplitude_modulation=0.2,
                   extra_theta=0.,
                   bigr=1.,
                   littler=0.3):

        if trace_lines and not self.traced_lines:
            self.trace_neumann_lines()
        
        import mayavi.mlab as may
        if clf:
            may.clf()

        plotarr = self.arr  # -1 factor corrects mayavi colours
        xnum, ynum = self.shape
        absmax = n.max((n.max(plotarr), n.abs(n.min(plotarr))))
        modulated_plotarr = amplitude_modulation * plotarr / absmax

        r = littler
        pi = n.pi
        cos = n.cos
        sin = n.sin
        phi, theta = n.mgrid[0:2*pi:xnum*1j, 0:2*pi:ynum*1j]
        theta += extra_theta

        
        x = (bigr + (r+modulated_plotarr)*sin(theta)) * sin(phi)
        y = (bigr + (r+modulated_plotarr)*sin(theta)) * cos(phi)
        z = (r+modulated_plotarr)*cos(theta)

        mesh = may.mesh(x, y, z, scalars=plotarr, colormap=cmap)

        # Invert colour scheme - flipping array doesn't seem to work
        lut1 = mesh.module_manager.scalar_lut_manager.lut.table.to_array()
        lut1 = lut1[::-1]
        mesh.module_manager.scalar_lut_manager.lut.table = lut1

        if trace_lines:
            lines = self.lines
            for index, line in enumerate(lines):
                self.vprint('\r\tPlotting line {} / {} on torus'.format(index, len(lines)), False)
                segs = sanitise_line(line)
                for seg in segs:
                    xs = seg[:, 0] / xnum * 2*n.pi
                    ys = seg[:, 1] / ynum * 2*n.pi
                    local
                    xpoints = (bigr + r*sin(xs)) * sin(ys)
                    ypoints = (bigr + r*sin(xs)) * cos(ys)
                    zpoints = r * cos(xs)
                    may.plot3d(xpoints, ypoints, zpoints, color=(1, 0, 1),
                               tube_radius=0.01)
            self.vprint()
             

    def plot3d(self, clf=True, save=''):
        import mayavi.mlab as may
        if clf:
            may.clf()

        if not self.arr_filled:
            self.fill_arr()
        if not self.found_crits:
            self.find_critical_points()
        # if not self.traced_lines and trace_lines:
        #     self.trace_neumann_lines()
        # if not self.hessian_filled and plot_hessian:
        #     self.make_hessian_array()
        # if not self.found_domains and show_domain_patches:
        #     self.get_recognised_domains()

        plotarr = n.rot90(self.arr[::-1], 3)
        plotarr = self.arr

        #fig, ax = plot_arr_with_crits(plotarr, self.crits)
        maxima, minima, saddles, degenerate = self.crits
        maxima = n.array(maxima)
        minima = n.array(minima)
        saddles = n.array(saddles)
        degenerate = n.array(degenerate)

        may.surf(plotarr, colormap='RdYlBu')
        may.contour_surf(plotarr, contours=[0])

        hmin = n.min(plotarr)
        hmax = n.max(plotarr)
        dh = (hmin-hmax)/20

        saddlezs = []
        if len(saddles) > 0:
            for saddle in saddles:
                x, y = saddle
                saddlezs.append(self.arr[x, y])
        saddlezs = n.array(saddlezs)

        maximumzs = []
        if len(maxima) > 0:
            for maximum in maxima:
                x, y = maximum
                maximumzs.append(self.arr[x, y])
        maximumzs = n.array(maximumzs)

        minimumzs = []
        if len(minima) > 0:
            for minimum in minima:
                x, y = minimum
                minimumzs.append(self.arr[x, y])
        minimumzs = n.array(minimumzs)

        extlen = self.arr.shape[0]
        saddleextent = [-extlen/2.0, extlen/2.0, -extlen/2.0,
                        extlen/2.0, n.min(saddlezs), n.max(saddlezs)]
        minimumextent = [-extlen/2.0, extlen/2.0, -extlen/2.0,
                         extlen/2.0, n.min(minimumzs), n.max(minimumzs)]
        maximumextent = [-extlen/2.0, extlen/2.0, -extlen/2.0,
                         extlen/2.0, n.min(maximumzs), n.max(maximumzs)]

        if len(saddles) > 0:
            may.points3d(saddles[:, 0], saddles[:, 1], saddlezs,
                         color=(1, 1, 0), scale_factor=1.5, extent=saddleextent)
        if len(minima) > 0:
            may.points3d(minima[:, 0], minima[:, 1], minimumzs,
                         color=(1, 0, 0), scale_factor=1.5,
                         extent=minimumextent)
        if len(saddles) > 0:
            may.points3d(maxima[:, 0], maxima[:, 1], maximumzs,
                         color=(0, 0, 1), scale_factor=1.5,
                         extent=maximumextent)

        if save:
            may.savefig(save)

def get_filled_array(xnum, ynum, dx, dy, func, start_point=(0.0, 0.0)):
    arr = n.zeros((xnum, ynum), dtype=n.float64)
    sx, sy = start_point
    for x in range(xnum):
        for y in range(ynum):
            arr[x, y] = func(sx + x*dx, sy + y*dy)
    return arr

def trace_gradient_line(sx, sy, dx, dy, xnum, ynum, func,
                        critdict, start_point, direction, to_edges,
                        func_params=(), area_constraint=None,
                        integer_saddles=True, nearby_distance=1.):
    '''Trace gradient (Neumann) line at point until reaching a critical point.
    func_params is ignored by this python implementation (but may be used by cython).
    '''
    cx, cy = sx, sy
    startx, starty = start_point
    if direction == 'down':
        dirfac = -1
    else:
        dirfac = 1

    points = [[cx, cy]]

    while True:
        gradient = grad(func, startx+cx*dx, starty+cy*dy, dx, dy) * dirfac
        angle = n.arctan2(gradient[1], gradient[0])

        cx += 0.25*n.cos(angle)
        cy += 0.25*n.sin(angle)

        if len(points)>20:
            if mag(n.array([cx, cy])-n.array(points[-20])) < 0.75:
                #print 'new points', [int(n.round(cx)), int(n.round(cy))]
                return (points, [int(n.round(cx)), int(n.round(cy))])

        # if len(points) == 3:
        #     print 'extra test'
        #     p = n.array(points[-3])
        #     dp = n.array(points[-2]) - n.array(points[-3])
        #     q = n.array(points[-1])
        #     dq = n.array([cx, cy])-n.array(q)
        #     print p, dp, q, dq
        #     cross, wx, wy = do_vectors_cross(p, dp, q, dq)
        #     print cross, wx, wy
        #     if cross:
        #         cx -= 0.5*n.cos(angle)
        #         cy -= 0.5*n.sin(angle)
        points.append([cx, cy])

        if (area_constraint is not None and not area_constraint(cx, cy)):
            return (points, None)
        elif cx < 0 or cx > xnum or cy < 0 or cy > ynum:
            if to_edges in ['periodic']:
                cx %= xnum
                cy %= ynum
            else:
                return (points, None)

        nearx, neary = int(n.round(cx)), int(n.round(cy))
        if to_edges in ['periodic']:
            nearx %= xnum
            neary %= ynum
        if integer_saddles:
            if (nearx, neary) in critdict:
                crit_type = critdict[nearx, neary]
                if ((crit_type == 'maximum' and direction == 'down') or
                    (crit_type == 'minimum' and direction == 'up')):
                # if crit_type in ['maximum','minimum']:
                    #print (nearx, neary), crit_type, direction
                    points.append((nearx, neary))
                    return (points, (nearx, neary))
        else:
            keys = critdict.keys()
            for key in keys:
                distance = reduce_distance((cx, cy), key, xnum, ynum)
                if distance < nearby_distance:
                    print 'distance nearby!', distance
                    points.append(key)
                    return (points, key)

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

def grad(func, x, y, dx, dy):
    dfdx = (func(x, y)-func(x+0.015*dx, y))/(0.015*dx)
    dfdy = (func(x, y)-func(x, y+0.015*dy))/(0.015*dy)
    return n.array([dfdx, dfdy])

def hessian(func, x, y, dx, dy):
    dfdx, dfdy = grad(func, x, y, dx, dy)

    dfdxdx = (grad(func, x+0.05*dx, y, dx, dy)[0] - dfdx) / (0.05*dx)
    dfdydy = (grad(func, x, y+0.05*dy, dx, dy)[1] - dfdy) / (0.05*dy)
    dfdxdy = (grad(func, x+0.05*dx, y, dx, dy)[1] - dfdy) / (0.05*dx)
    dfdydx = (grad(func, x, y+0.05*dy, dx, dy)[0] - dfdx) / (0.05*dy)

    return n.array([[dfdxdx, dfdxdy], [dfdydx, dfdydy]])

def hessian_det(func, x, y, dx, dy):
    hess_mat = hessian(func, x, y, dx, dy)
    return n.linalg.det(hess_mat)

def hessian_sign(func, x, y, dx, dy):
    hess_mat = hessian(func, x, y, dx, dy)
    return n.sign(n.linalg.det(hess_mat))

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

even_adj_indices = n.array([[-1, -1], [-1, 0], [0, 1], [1, 0], [1, -1], [0, -1]])
odd_adj_indices =  n.array([[-1, 0], [-1, 1], [0, 1], [1, 1], [1, 0], [0, -1]])
all_adj_indices = n.array([[-1, -1], [-1, 0], [-1, 1], [0, 1],
                           [1, 1], [1, 0], [1, -1], [0, -1]])
safe_adj_indices = n.array([[-1, 0], [0, 1], [1, 0], [0, -1]])
def get_critical_points(arr, to_edges=False, verbose=True):
    lx, ly = arr.shape
    adjs = n.zeros(6, dtype=n.float64)

    maxima = []
    minima = []
    saddles = []
    degenerate = []

    border_mult = n.ones(6, dtype=n.float64)

    prevx = -1
    for x, y in product(xrange(lx), xrange(ly)):
        if x != prevx:
            prevx = x
            if verbose:
                lineprint('\r\tx = {0} / {1}'.format(x, lx), False)
        if to_edges or (x != 0 and y != 0 and x != (lx-1) and y != (ly-1)):
            val = arr[x, y]
            if x % 2 == 0:
                ais = even_adj_indices.copy()
            else:
                ais = odd_adj_indices.copy()

            ais[:, 0] += x
            ais[:, 1] += y

            if to_edges == 'periodic':
                ais[:, 0] = ais[:, 0] % lx
                ais[:, 1] = ais[:, 1] % ly


            for i in range(6):
                adjs[i] = arr[ais[i][0] % lx, ais[i][1] % ly]

            point_type = classify_point(adjs-val)

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
                print ('=> odd number of sign changes, perhaps the function'
                       'is symmetrical about this point.')


    if verbose:
        lineprint()

    return (maxima, minima, saddles, degenerate)
    #return (n.array(maxima), n.array(minima),
#                   n.array(saddles), n.array(degenerate))

def classify_point(ds):
    if n.all(ds > 0):
        return 'minimum'
    elif n.all(ds < 0):
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

def plot_arr_with_crits(arr, crits):
    maxima, minima, saddles, degenerate = crits
    maxima = n.array(maxima)
    minima = n.array(minima)
    saddles = n.array(saddles)
    degenerate = n.array(degenerate)

    fig, ax = plt.subplots()

    ax.imshow(arr, cmap='Spectral_r', interpolation='none', alpha=0.6)

    legend_entries = []

    if len(maxima) > 0:
        ax.scatter(maxima[:, 0], maxima[:, 1], 60, c='r')
        legend_entries.append('maxima')
    if len(minima) > 0:
        ax.scatter(minima[:, 0], minima[:, 1], 60, c='b')
        legend_entries.append('minima')
    if len(saddles) > 0:
        ax.scatter(saddles[:, 0], saddles[:, 1], 60, color='yellow')
        legend_entries.append('saddles')
    if len(degenerate) > 0:
        ax.scatter(degenerate[:, 0], degenerate[:, 1], c='orange')
        legend_entries.append('degenerate')

    return fig, ax

def sanitise_line(l):
    curcut = 0
    segs = []
    for i in range(len(l)-1):
        next = l[i+1]
        cur = l[i]
        #print n.abs(next[0]-cur[0])
        #print n.abs(next[0]-cur[0]), n.abs(next[1]-cur[1])
        if n.abs(next[0]-cur[0]) > 5 or n.abs(next[1]-cur[1]) > 5:
            #print 'bigchange', len(segs)
            segs.append(l[curcut:(i+1)])
            #print '->', len(segs)
            curcut = i+1
    if curcut < len(l):
        segs.append(l[curcut:])
    return segs


def random_wave_function(number=50, wvmag=5, seed=0, returnall=False):
    seed = random.randint(0,10000000)
    generator = n.random.RandomState()
    generator.seed(seed)

    amps = generator.normal(size=number)
    phases = 2*n.pi*generator.rand(number)
    wvs = n.zeros((number, 2), dtype=n.float64)
    for i in range(number):
        wv = generator.normal(size=2)
        wv /= mag(wv)
        wv *= wvmag
        wvs[i] = wv

    def func(x, y):
        res = 0.0
        xyarr = n.array([x, y], dtype=n.float64)
        interior = wvs.dot(xyarr) + phases
        exterior = amps*n.sin(interior)
        return n.sum(exterior)
    if returnall:
        return (func, (amps, wvs, phases))
    return func

def duofactors(k):
    outs = []
    for i in range(int(n.sqrt(k))+2):
        for j in range(i, int(n.sqrt(k))+2):
            if (i**2 + j**2) == k:
                outs.append((i, j))
    return outs

def range_factors(a, b=None, best=False):
    if b is None:
        b = a
        a = 0
    results = zip(range(a, b), map(duofactors, range(a, b)))
    if best:
        bestnum = max([len(j[1]) for j in results])
        results = filter(lambda j: len(j[1]) == bestnum, results)
    return results

sign_permutations = set(chain(permutations([1, 1]),
                              permutations([1, -1]),
                              permutations([-1, -1])))
def periodic_random_wave_function(energy, seed=0, returnall=False,
                                  compiled=True, spectrum='shell'):
    if seed == 0:
        seed = n.random.randint(100000000)

    generator = n.random.RandomState()
    generator.seed(seed)

    wvs = []
    amps = []
    phases = []

    if spectrum == 'shell':
        possible_wvs = duofactors(energy)
        if len(possible_wvs) == 0:
            raise Exception('Energy must be a sum of two squares.')
        for wv in possible_wvs:
            for permutation in permutations(wv):
                for signs in sign_permutations:
                    wvs.append(
                        [p*q for (p, q) in zip(signs, permutation)])
                    amps.append(generator.normal(0., 1.))
                    phases.append(generator.rand() * 2 * n.pi)

    elif spectrum == 'gaussian':
        possible_wvs = []
        for i in range(int(-3.5*energy)-1, int(3.5*energy)):
            lineprint(
                'Populating Gaussian spectrum, i = {} from {} -> {}'.format(
                i, int(-3.5*energy-1), int(3.5*energy)))
            for j in range(int(-3.5*energy), int(3.5*energy)):
                if ((i != 0 or j != 0) and
                    (i % 2 and j % 2)):
                    wvs.append([i, j])
                    amps.append(generator.normal(0., energy) *
                                1/n.sqrt(energy) * n.exp(
                                    -0.5*(i**2 + j**2)/energy))
                    phases.append(generator.rand() * 2 * n.pi)
        print
                               
                                               
                                

    wvs = n.array(wvs)
    amps = n.array(amps)
    phases = n.array(phases)

    if compiled and cneu is not None:
        func = partial(cneu.random_wave_function, wvs.astype(n.double),
                       amps.astype(n.double), phases.astype(n.double))

    else:
        def func(x, y):
            res = 0.0
            xyarr = n.array([x, y], dtype=n.float64)
            interior = wvs.dot(xyarr) + phases
            exterior = amps*n.sin(interior)
            return n.sum(exterior)
    if returnall:
        return (func, (amps, wvs, phases))
    return func

def mag(v):
    return n.sqrt(v.dot(v))

def lineprint(s='', newline=True):
    sys.stdout.write(s)
    if newline:
        sys.stdout.write('\n')
    sys.stdout.flush()


def get_saddle_directions(vals):
    changes = []
    for i in range(len(vals)):
        cur = vals[i]
        next = vals[(i+1) % len(vals)]
        if n.sign(next) != n.sign(cur):
            changes.append((i+1) % len(vals))

    returns = []
    if len(changes) == 4:
        region_1 = vals[changes[0]:changes[1]]
        region_2 = vals[changes[1]:changes[2]]
        region_3 = vals[changes[2]:changes[3]]
        region_4a = vals[changes[3]:]
        region_4b = vals[:changes[0]]
        if len(region_4a) > 0 and len(region_4b) > 0:
            region_4 = n.hstack((region_4a, region_4b))
        elif len(region_4a) > 0:
            region_4 = region_4a
        else:
            region_4 = region_4b

        if n.sign(region_1[0]) > 0:
            returns.append((1.0, n.argmax(region_1) + changes[0]))
        else:
            returns.append((-1.0, n.argmax(region_1) + changes[0]))
        if n.sign(region_2[0]) > 0:
            returns.append((1.0, n.argmax(region_2) + changes[1]))
        else:
            returns.append((-1.0, n.argmax(region_2) + changes[1]))
        if n.sign(region_3[0]) > 0:
            returns.append((1.0, n.argmax(region_3) + changes[2]))
        else:
            returns.append((-1.0, n.argmax(region_3) + changes[2]))
        if n.sign(region_4[0]) > 0:
            returns.append((1.0, (n.argmax(region_4) + changes[3]) % len(vals)))
        else:
            returns.append(
                (-1.0, (n.argmax(region_4) + changes[3]) % len(vals)))
    else:
        print 'Saddle doesn\'t have 4 sign changes?'
        print vals
        print changes
    return returns

def angle_index(line, lines):
    endseg = line[-2:][::-1]
    angle = n.arctan2(endseg[1, 1]-endseg[0, 1], endseg[1, 0]-endseg[0, 0])
    angles = n.zeros(len(lines), dtype=n.float64)
    for i in range(len(lines)):
        curline = lines[i]
        seg = curline[:2]
        angles[i] = n.arctan2(seg[1, 1]-seg[0, 1], seg[1, 0]-seg[0, 0])
    if n.any(angle > angles):
        return n.argmax(angle > angles)
    else:
        return 0

def area_from_border(lines):
    area = 0.0
    num = len(lines)
    for i in range(len(lines)):
        area += (lines[i, 0]*lines[(i+1) % num, 1] -
                 lines[i, 1]*lines[(i+1) % num, 0])
    return n.abs(area / 2.)

sanitised_domains = []
def sanitise_domain(ps, shape=None):
    '''Takes a domain that may pass through a boundary, and sanitises
    it to have the right shape.'''
    if shape is None:
        dx, dy = None, None
    else:
        dx, dy = shape
    ps = ps.copy()

    for i in range(len(ps)):
        cur = ps[i]
        nex = ps[(i+1) % len(ps)]
        if n.abs(nex[0]-cur[0]) > 5:
            diff = n.round(nex[0] - cur[0])
            ps[i+1:, 0] -= diff
        if n.abs(nex[1]-cur[1]) > 5:
            diff = n.round(nex[1] - cur[1])
            ps[i+1:, 1] -= diff
    sanitised_domains.append(ps)
    return ps


def crude_area_from_border(lines, numsteps=100):
    maxxs = map(lambda j: n.max(j[:, 0]), lines)
    minxs = map(lambda j: n.min(j[:, 0]), lines)
    maxys = map(lambda j: n.max(j[:, 1]), lines)
    minys = map(lambda j: n.min(j[:, 1]), lines)

    minx = n.min(minxs)
    maxx = n.max(maxxs)
    miny = n.min(minys)
    maxy = n.max(maxys)

    dx = (maxx - minx) / numsteps
    dy = (maxy - miny) / numsteps

    fullpath = n.vstack(map(lambda j: j[:-1], lines))

    area = 0.0

    for x in n.linspace(minx-0.005*dx, maxx+0.005*dx, numsteps):
        testline = n.array([[x, miny], [x, maxy]])
        #print 'testline', testline
        p = testline[0]
        r = testline[1] - testline[0]
        intersections = []
        for i in range(len(fullpath)):
            q = fullpath[i]
            s = fullpath[(i+1) % len(fullpath)] - q
            #print minx, maxx, miny, maxy,' -> ', p, r, q, s
            intersect, t, u = do_vectors_cross(p, r, q, s)
            if intersect:
                intersect_y = (p + u*r)[1]
                intersections.append(intersect_y)
        if len(intersections) % 2 != 0:
            print 'Number of intersections is not even!'
        elif len(intersections) == 2:
            bottom = n.min(intersections)
            top = n.max(intersections)
            area += dx*(top-bottom)
        elif len(intersections) == 4:
            intersections = n.sort(intersections)
            area += dx*(intersections[1]-intersections[0] +
                        intersections[3]-intersections[2])
        #print len(intersections),'intersections found', intersections
    return area

    return (minx, maxx, miny, maxy)


def img_interpolate(x, y, arr):
    roundx = int(n.floor(x))
    roundy = int(n.floor(y))
    dx = roundx - x
    dy = roundy - y

    return (arr[roundx, roundy]*(1-dx)*(1-dy) +
            arr[roundx+1, roundy]*dx*(1-dy) +
            arr[roundx, roundy+1]*(1-dx)*dy +
            arr[roundx+1, roundy+1]*dx*dy)


def do_vectors_cross(p, r, q, s):
    """
    Takes four vectors p, dp and q, dq, then tests whether they cross in
    the dp/dq region. Returns this boolean, and the point where the
    crossing actually occurs.
    """
    p = n.float64(p)
    r = n.float64(r)
    q = n.float64(q)
    s = n.float64(s)
    if n.abs(n.cross(r, s)) < 0.00001:
        return (False, 0., 0.)

    t = n.cross(q-p, s) / n.cross(r, s)

    if 0.0 < t < 1.0:
        u = n.cross(q-p, r) / n.cross(r, s)
        if 0.0 < u < 1.0:
            return (True, t, u)
        return (False, t, u)
    return (False, t, -1.0)


def animate(filen='test', path='onephase', number=200, function=None):
    if function is None:
        f, d2 = random_wave_function(50, 10, returnall=True)
    else:
        f, d2 = function

    amps, wvs, phases = d2
    phases -= 2 * (2*n.pi)/50

    for i in range(number):
        if path == 'allphase':
            phases += 2*n.pi/number
        elif path[:9] == 'manyphase':
            val = int(path[9:])
            phases[:val] += 1/number * 4*(2*n.pi/50)
        else:
            phases[0] += 2*n.pi/number

        def func(x, y):
            res = 0.0
            xyarr = n.array([x, y], dtype=n.float64)
            interior = wvs.dot(xyarr) + phases
            exterior = amps*n.sin(interior)
            return n.sum(exterior)
        a = NeumannTracer(100, 100, n.pi/190, n.pi/190, f,
                          start_point=(-30*n.pi/190, -30*n.pi/190))
        a.plot(save='{0}_{1:04}.png'.format(filen, i))
        del a  # For some reason, the reference count builds up?


def periodic_animate(scale=5, number=50, frames=200, downscale=2, func=None,
                     plot3d=False):
    length = int(100/float(downscale) * float(scale)/5.)
    periodicity = n.sqrt(scale/2.)

    if func is None:
        f, d2 = periodic_random_wave_function(number, scale, returnall=True)
    else:
        f, d2 = func
    amps, wvs, phases = d2

    phases = phases.copy()
    phases += 25 * 2*n.pi/frames

    for i in range(frames):
        #phases += 2*n.pi/frames
        phases += (4*2*n.pi/frames) / frames

        def func(x, y):
            res = 0.0
            xyarr = n.array([x, y], dtype=n.float64)
            interior = wvs.dot(xyarr) + phases
            exterior = amps*n.sin(interior)
            return n.sum(exterior)
        a = NeumannTracer(length, length,
                          periodicity/float(length), periodicity/float(length),
                          func, to_edges='periodic')
        a.plot(save='panim1_{0:04}.png'.format(i))
        if plot3d:
            a.plot3d(save='panim3d_{0:04}.png'.format(i))



def get_periodic_tracer(energy, gridsize=None, downscale=2,
                        returnall=False, upsample=5,
                        seed=0, compiled=True,
                        pass_func_params=True,
                        spectrum='shell'):
    '''Returns a :class:`NeumannTracer` covering the periodic domain of a
    torus with the given scale and number of wavevectors.

    The sample grid size may be set manually with *gridsize*, or
    otherwise chosen automatically and scaled by downscale.'''
    if not duofactors(energy):
        raise ValueError('Input energy has no degenerate components.')
    if gridsize is None:
        length = int(100/float(downscale) * float(n.sqrt(energy))/3.)
    else:
        length = gridsize
    print 'length is', length
    periodicity = 2*n.pi
    f, d2 = periodic_random_wave_function(energy, returnall=True,
                                          seed=seed, compiled=compiled,
                                          spectrum=spectrum)
    if pass_func_params:
        func_params = ('rwm', d2[1].astype(n.double), d2[0].astype(n.double),
                       d2[2].astype(n.double))
    else:
        func_params = ()
    tracer = NeumannTracer(length, length,
                           periodicity/(1.*length), periodicity/(1.*length),
                           f,
                           to_edges='periodic', upsample=upsample,
                           func_params=func_params,
                           eigenvalue=energy)
    if returnall:
        return tracer, f, d2, length
    return tracer

def stadium_constraint(square_width, circle_radius, x, y):
    '''Returns True if x, y in the desymmetrised stadium with given
    parameters, else False.'''
    print 'stadium', square_width, circle_radius, x, y
    if x < 0 or x > square_width + circle_radius:
        return False
    if 0 <= x <= square_width and 0 <= y <= circle_radius:
        return True
    if square_width <= x <= (square_width + circle_radius):
        return True if 0 <= y <= (
            n.sqrt(circle_radius**2 - (x - square_width)**2)) else False
    return False

def get_stadium_constraint(square_width, circle_radius):
    '''Returns a stadium constraint function for square with the given
    width and height from the circle radius.'''
    return partial(stadium_constraint, square_width, circle_radius)

def get_stadium_tracer(square_width, circle_radius, coefficients):
    func = get_stadium_mode(coefficients)
    tracer = NeumannTracer(200, 200, 0.01, 0.005, func)
    xnum, ynum = tracer.shape
    dx, dy = tracer.dr
    sx, sy = tracer.start_point
    new_square_width = (square_width - sx) / dx
    new_circle_radius = circle_radius / dx
    tracer.area_constraint = get_stadium_constraint(new_square_width,
                                                    new_circle_radius)
    extremum = n.max((n.max(tracer.arr), n.abs(n.min(tracer.arr))))
    tracer.arr[0][-1] = extremum
    tracer.arr[0][-2] = -1*extremum
    return tracer

def get_stadium_tracer_from_file(number):
    import json
    with open('/home/asandy/neumann/stadium_mode_amplitudes.json',
              'r') as fileh:
        coeffs = json.load(fileh)[number]
    return get_stadium_tracer(1., 1., coeffs)
    
def get_stadium_mode(coefficients, prefactor=1.):
    assert len(coefficients) == 50
    k = 1.
    coefficients = n.array(coefficients)
    def func(x, y):
        angles = n.pi / 50. * (coefficients - 0.5)
        components = (coefficients *
                      n.cos(k * n.cos(angles) * x) *
                      n.cos(k * n.sin(angles) * y))
        return prefactor * n.sum(components)
    return func

def hermite(n, x):
    '''Computes the hermite polynomial with given n, at x'''
    h0 = 1
    h1 = 2*x
    i = 1

    if n == 0:
        return h0
    elif n == 1:
        return h1
    else:
        cur = h1
        prev = h0
        while i < n:
            new = 2*x*cur - 2*i*prev
            prev = cur
            cur = new
            i += 1
    return new

def get_random_hermite_mode(energy, adjust=[], seed=0):
    '''Returns a function of x and y.'''
    if seed == 0:
        seed = n.random.randint(10000000)
    g = n.random.RandomState()
    g.seed(seed)

    if cneu is not None:
        hermite_func = cneu.hermite
    else:
        hermite_func = hermite

    polys = []
    everything_prefac = 1/n.sqrt(2**energy)
    for k in range(0, energy+1):
        prefac = 1 / n.sqrt(factorial(k)*factorial(energy-k))
        hx = partial(hermite_func, k)
        hy = partial(hermite_func, energy-k)
        polys.append([g.normal() * prefac, hx, hy])

    for index, value in adjust:
        polys[index][0] = polys[index][0] * value

    def sumfunc(x, y):
        ans = 0.
        radial_prefactor = n.exp(-0.5*(x**2 + y**2))
        for i in range(len(polys)):
            row = polys[i]
            ans += row[0] * row[1](x) * row[2](y)
        return ans * radial_prefactor
    return n.vectorize(sumfunc)

def do_domain_statistics_at_scales(scales, domains=1000, downscale=1, filen='neumann_statistics'):
    '''For every scale in scales, generates random torus eigenfunctions at
    that scale, counting the domains unitl it has found at least the
    number specified as an argument.

    downscale is a resolution parameter, it can usually be ignored.
    '''
    results = {}

    for scale in scales:
        domains_found = 0
        while domains_found < domains:
            lineprint('energy {}: found {} / {} domains'.format(
                scale, domains_found, domains))
            tracer, results = do_domain_statistics_at_scale(
                scale, downscale)
            domains_found += 2*results['count'][2]

            cur_filen = get_next_filen(filen + '_{}_'.format(scale))
            with open(cur_filen, 'wb') as fileh:
                cPickle.dump(results, fileh)
    print
                          

def do_domain_statistics_at_scale(scale, downscale=1):
    areas = []
    perimeters = []
    rhos = []
    degree_dists = n.zeros((8,2))
    diameters = []

    a = get_periodic_tracer(scale, downscale=downscale)
         
    if a is None:
        return None
    areas = a.get_domain_areas()
    perimeters = a.get_domain_perimeters()
    rhos = a.get_domain_rhos()
    rhos_by_type = a.get_domain_rhos_by_type()
    degree_dists = a.get_critical_degree_dists()
    diameters = a.get_domain_diameters()
    count = map(len, a.crits)

    results = {'areas': areas, 'perimeters': perimeters,
               'rhos': rhos,
               'rhos_by_type': rhos_by_type,
               'degree_dists': degree_dists,
               'diameters': diameters,
               'energy': scale, 'count': count}

    return a, results

def load_domain_statistics_from_filens(filens, flatten=True):
    '''For each filen, tries to load a domain statistics
    dictionary, and returns a list indexed by energy.'''
    results = {}
    for filen in filens:
        try:
            with open(filen, 'rb') as fileh:
                data = cPickle.load(fileh)
            energy = data['energy']
            if energy not in results:
                results[energy] = []
            results[energy].append(data)
        except cPickle.UnpicklingError:
            pass

    if flatten:
        results = flatten_results(results)

    return results

def flatten_results(results):
    '''Takes a results dictionary and merges all the constituent results.'''
    new_results = {}
    for scale, energy_results in results.items():
        areas = []
        perimeters = []
        rhos = []
        degree_dists = n.zeros((8, 2), dtype=n.float)
        diameters = []
        count = []
        rhos_by_type = {'lens': [], 'star': [], 'wedge': [], 'bad': []}
        tracer_number = 0
        for tracer_results in energy_results:
            areas.append(tracer_results['areas'])
            perimeters.append(tracer_results['perimeters'])
            rhos.append(tracer_results['rhos'])
            degree_dists += tracer_results['degree_dists']
            diameters.append(tracer_results['diameters'])
            count.append(tracer_results['count'])
            if 'rhos_by_type' in tracer_results:
                rhos_by_type['lens'].append(
                    tracer_results['rhos_by_type']['lens'])
                rhos_by_type['star'].append(
                    tracer_results['rhos_by_type']['star'])
                rhos_by_type['wedge'].append(
                    tracer_results['rhos_by_type']['wedge'])
                rhos_by_type['bad'].append(
                    tracer_results['rhos_by_type']['bad'])
                
            tracer_number += 1
        areas = n.hstack(areas)
        perimeters = n.hstack(perimeters)
        rhos = n.hstack(rhos)
        degree_dists /= float(tracer_number)
        diameters = n.hstack(diameters)
        count = n.sum(n.array(count), axis=0)
        new_results[scale] = {'areas': areas, 'perimeters': perimeters,
                              'rhos': rhos,
                              'rhos_by_type': rhos_by_type,
                              'degree_dists': degree_dists,
                              'diameters': diameters,
                              'energy': scale, 'count': count}
    return new_results
             
def make_plots_from_results(results, filen='results_plots', bins=30):
    energies = sorted(results.keys())
    fig_size = (5, 4)
    mpl_interactive(False)

    legend = []
    fig, ax = plt.subplots()
    for energy, areas in [(energy, results[energy]['areas']) for
                          energy in energies]:
        ys, xs = n.histogram(areas*energy, bins=n.linspace(0, 30, bins),
                             density=True)
        dx = 0.5*(xs[1] - xs[0])
        ax.plot(xs[:-1] + dx, ys)
        legend.append(str(energy))
    ax.legend(legend)
    ax.set_xlabel('area $\\times\lambda$')
    ax.set_ylabel('PDF')
    ax.set_xlim(0, 30)
    fig.set_size_inches(fig_size)
    fig.tight_layout()
    fig.savefig(filen + '_areas.png', dpi=200)

    legend = []
    fig, ax = plt.subplots()
    for energy, perimeters in [(energy, results[energy]['perimeters']) for
                          energy in energies]:
        ys, xs = n.histogram(perimeters*n.sqrt(energy),
                             bins=n.linspace(0, 30, bins), density=True)
        dx = 0.5*(xs[1] - xs[0])
        ax.plot(xs[:-1] + dx, ys)
        legend.append(str(energy))
    ax.legend(legend)
    ax.set_xlabel('perimeter $\\times\sqrt{\lambda}$')
    ax.set_ylabel('PDF')
    ax.set_xlim(0, 30)
    fig.set_size_inches(fig_size)
    fig.tight_layout()
    fig.savefig(filen + '_perimeters.png', dpi=200)

    legend = []
    fig, ax = plt.subplots()
    for energy, diameters in [(energy, results[energy]['diameters']) for
                          energy in energies]:
        ys, xs = n.histogram(diameters*n.sqrt(energy),
                             bins=n.linspace(0, 10, bins),
                             density=True)
        dx = 0.5*(xs[1] - xs[0])
        ax.plot(xs[:-1] + dx, ys)
        legend.append(str(energy))
    ax.legend(legend)
    ax.set_xlabel('diameter $\\times\sqrt{\lambda}$')
    ax.set_ylabel('PDF')
    ax.set_xlim(0, 15)
    fig.set_size_inches(fig_size)
    fig.tight_layout()
    fig.savefig(filen + '_diameters.png', dpi=200)
    
    legend = []
    fig, ax = plt.subplots()
    for energy, rhos in [(energy, results[energy]['rhos']) for
                          energy in energies]:
        ys, xs = n.histogram(rhos, bins=n.linspace(0, 1.5, bins),
                             density=True)
        dx = 0.5*(xs[1] - xs[0])
        ax.plot(xs[:-1] + dx, ys)
        legend.append(str(energy))
    ax.vlines([0.785, 0.92], 0, 3.0)
    ax.legend(legend)
    ax.set_xlabel('$\\rho$')
    ax.set_ylabel('PDF')
    fig.set_size_inches(fig_size)
    fig.tight_layout()
    fig.savefig(filen + '_rhos.png', dpi=200)

    legend = []
    fig, ax = plt.subplots()
    for energy, degree_dists in [(energy, results[energy]['degree_dists']) for
                          energy in energies]:
        ax.plot(degree_dists[:, 0], degree_dists[:, 1])
        legend.append(str(energy))
    ax.legend(legend)
    ax.set_xlabel('degree')
    ax.set_ylabel('PDF')
    fig.set_size_inches(fig_size)
    fig.tight_layout()
    fig.savefig(filen + '_degrees.png', dpi=200)

    mpl_interactive(True)

        

def nearly_integer(v, cutoff=0.0001):
    diff = n.abs(int(v) - v)
    return True if diff < cutoff else False

def get_degree_statistics(scale=65, downscale=7, runs=100):
    threes = []
    fours = []
    fives = []
    sixs = []
    sevens = []
    eights = []
    nines = []
    for i in range(runs):
        tracer = get_periodic_tracer(scale, downscale=downscale)
        degrees = tracer.get_critical_degrees()
        mindegs, maxdegs = degrees
        mindegs = n.array(mindegs)
        threes.append(n.sum(mindegs==3))
        fours.append(n.sum(mindegs==4))
        fives.append(n.sum(mindegs==5))
        sixs.append(n.sum(mindegs==6))
        sevens.append(n.sum(mindegs==7))
        eights.append(n.sum(mindegs==8))
        nines.append(n.sum(mindegs==9))

    return [threes, fours, fives, sixs, sevens, eights, nines]

def get_next_filen(filen):
    '''Returns next iteration of filen that doesn't exist.'''
    i = 1
    while os.path.exists('{}_{:05d}.pickle'.format(filen, i)):
        i += 1
    return '{}_{:05d}.pickle'.format(filen, i)


def reduce_distance(p1, p2, xnum, ynum):
    '''Takes two positions and finds the shortest distance between
    them by taking account of periodic boundary conditions.
    '''
    p1 = n.array(p1)
    p2 = n.array(p2)
    distance = mag(p2 - p1)
    x1, y1 = p1
    x2, y2 = p2
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            new_distance = mag((p2 + n.array([dx*xnum, dy*ynum])))
            if new_distance < distance:
                distance = new_distance
    return distance

def angle_of(v):
    '''Calculates the angle of the given vector.'''
    return n.arctan2(v[1], v[0])

def degree_1_critical_point(x, y):
    
    if n.isnan(x) or n.isnan(y):
        raise Exception('nan')
        
    epsi = .1;
    x0 = 0.7;  y0 = n.pi/4+.45;  
    rx =  .80*n.pi/4;   ry = .6*n.pi/4;
    amp_bump = 10**3;

    if (n.abs(x-x0) < rx) and (n.abs(y-y0) < ry):
        return (n.sin(2*y) + epsi*n.cos(2*x)) + amp_bump * n.exp(-1/(rx**2-(x-x0)**2)-1/(ry**2-(y-y0)**2)) 
    else:
        return (n.sin(2*y) + epsi*n.cos(2*x))
