import numpy as n
from itertools import product
from colorsys import hsv_to_rgb

from matplotlib import interactive as mpl_interactive
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

mpl_linestyles = ['', ' ', 'None', '--', '-.', '-', ':']
patch_linestyles = ['solid','dashed','dashdot','dotted']

class UpsampleRegion(object):
    def __init__(self, sx, sy, lx, ly, degree=1):
        self.sx = sx
        self.sy = sy
        self.start = n.array([x, y])
        self.lx = lx
        self.ly = ly
        self.degree = degree
    def shape_parameters(self):
        return (self.sx, self.sy, self.lx, self.ly)
    def contains(self, x, y):
        sx, sy, lx, ly = self.shape_parameters()
        if (sx < x < (sx + lx)) and (sy < y < (sy+ly)):
            return True
        return False

class CriticalGraph(dict):
    '''A dict subclass representing the graph of Neumann nodes (critical
    points) and their connections (Neumann gradient lines). Nodes are
    stored via their coordinates, with node information as [node type,
    [list of Neumann lines]].

    Nodes should be added with add_or_edit_node rather than by direct
    dict setting.

    '''
    def __init__(self, xnum=0, ynum=0, *args):
        '''Arguments xnum and ynum should be the array size in each direction,
        and are used when calculating angles between critical
        points.

        '''
        super(CriticalGraph, self).__init__(*args)
        self.xnum = xnum
        self.ynum = ynum
        self.nodes_aligned = False

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
            print angles
            order = n.argsort(angles)
            newlines = []
            for i in range(len(order)):
                newlines.append(lines[order[i]])
            self[key][1] = newlines
        self.nodes_aligned = True

    def get_closed_domains(self):
        '''
        Returns a list of all closed domains, recognised by walking
        clockwise or anticlockwise from every saddle point.

        Explicitly fails when a node passes through a cell boundary or
        includes a line that returns to its own node.
        '''
        if not self.nodes_aligned:
            self.align_nodes()
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
        return unique_domains

    def get_domain_areas(self):
        domains = self.get_closed_domains()
        areas = []
        i = 0
        for domain in domains:
            lineprint('\rGetting area of domain %d / %d' % (i, len(domains)),
                      False)
            i += 1
            area = domain.area()
            areas.append(area)
        print # Lineprint newline
        return areas

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
        return NeumannDomain(lines)

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

class NeumannDomain(object):
    '''Represents a Neumann domain by storing a list of boundary
    NeumannLines. Also stores the number of saddles, maxima and minima
    participating in the domain.

    '''
    def __init__(self, lines):
        self.lines = lines
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

    def as_closed_curve(self):
        '''Joins the component lines and returns a single 2d array of points
        making up the domain boundary.

        '''
        points = []
        for line in self.lines:
            points.append(line.points)
        return n.vstack(points)

    def as_sanitised_curves(self):
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
the cell crosses a boundary

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
        return area_from_border(self.as_sanitised_curve())

    def crude_area(self):
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
    - xnum: Number of x pixels to sample
    - ynum: Number of y pixels to sample
    - dx: Step length dx for each pixel
    - dy: Step length dy for each pixel
    - func: A function to look for critical points in
    - start_point: The (x, y) tuple to take as (0, 0) in the function.
    - to_edges: May be False (ignore edges),
                'periodic' (use periodic boundary conditions) or
                'fourier' (periodic with a -1 factor if passing
                           through the boundary).

    '''
    def __init__(self, xnum, ynum, dx, dy, func,
                 start_point=(0.0, 0.0), to_edges=False, upsample_at=1000000):
        self.arr = n.zeros((xnum, ynum), dtype=n.float64)
        self.hessian_arr = n.zeros((xnum, ynum), dtype=n.float64)
        self.xnum = xnum
        self.ynum = ynum
        self.shape = (xnum, ynum)
        self.dx = dx
        self.dy = dy
        self.dr = (dx, dy)
        self.func = func
        self.start_point = start_point
        self.sx = start_point[0]
        self.sy = start_point[1]
        self.upsample_at = upsample_at

        self.to_edges = to_edges

        self.arr_filled = False
        self.found_crits = False
        self.traced_lines = False
        self.hessian_filled = False
        self.graph_built = False
        self.found_domains = False

        self.graph = CriticalGraph(xnum, ynum)
        self.domains = []

        self.upsample_regions = []

        self.lines = []
        self.start_points = []
        self.end_points = []
        self.extra_lines = []

        self.fill_arr()

        self.figax = (None, None)

    def func_at_coord(self, x, y):
        '''Return the value of self.func at abstract coordinates x, y
        translated via dx and the correct start point.

        '''
        sx, sy = self.start_point
        dx, dy = self.dr
        return self.func(sx + x*dx, sy + y*dy)

    def fill_arr(self):
        '''
        Sample the function on a (self.xnum, self.ynum) array.

        Result stored in self.arr.
        '''
        print 'Populating function sample array...'
        arr = self.arr
        sx, sy = self.start_point
        dx, dy = self.dr
        xnum, ynum = self.shape
        for x in range(xnum):
            lineprint('\r\tx = {0} / {1}'.format(x, xnum), False)
            for y in range(ynum):
                arr[x, y] = self.func(sx + x*dx, sy + y*dy)
        print
        self.arr_filled = True

    def find_critical_points(self):
        '''
        Iterate over self.arr walking about each pixel and checking the
        number of sign changes. Bin the result appropriately as a
        maximum, minimum, saddle or regular point, and store in
        self.crits.
        '''
        if not self.arr_filled:
            self.fill_arr()
        print 'Finding critical points...'
        maxima, minima, saddles, degenerate = get_critical_points(self.arr,
                                                                  self.to_edges)
        self.crits = (maxima, minima, saddles, degenerate)
        #self.prune_critical_points()
        self.maxima = self.crits[0]
        self.minima = self.crits[1]
        self.saddles = self.crits[2]
        self.degenerate = self.crits[3]

        self.crits_dict = critical_points_to_index_dict(self.crits)

        self.found_crits = True

    def critical_point_proximity_alert(self):
        pass

    def prune_critical_points(self):
        '''Detects and removes points where two extrema border a single
        saddle - it seems that these are (almost?) always numerical
        mistakes.

        This approach is a little crude, a later replacement might
        involve resampling.

        '''
        print 'Pruning critical points'
        maxima, minima, saddles, degenerate = self.crits

        tupmaxima = [tuple(j) for j in maxima]
        tupminima = [tuple(j) for j in minima]

        realmaxima = n.ones(len(maxima), dtype=bool)
        realminima = n.ones(len(minima), dtype=bool)
        realsaddles = n.ones(len(saddles), dtype=bool)

        for i in range(len(saddles)):
            lineprint('\r\tChecking saddle {0} / {1}'.format(i,
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
                    #print saddle, coords, tupmaxima
                    adjmax.append(coords)
                if tuple(coords) in tupminima:
                    adjmin.append(coords)
                    #print saddle, coords, tupminima
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
        lineprint()
        maxima = n.array(maxima)[realmaxima]
        maxima = [tuple(c) for c in maxima]
        minima = n.array(minima)[realminima]
        minima = [tuple(c) for c in minima]
        saddles = n.array(saddles)[realsaddles]
        saddles = [tuple(c) for c in saddles]
        self.crits = maxima, minima, saddles, degenerate

    def trace_neumann_lines(self):
        '''For every saddle in self.crits, drop 4 Neumann lines at the points
        of adjacent sign change, each gradient ascending/descending
        appropriately until they hit another critical point or appear
        to have stopped. The resulting lines are stored in self.lines.

        '''
        if not self.arr_filled:
            self.fill_arr()
        if not self.found_crits:
            self.find_critical_points()

        print 'Tracing Neumann lines...'

        arr = self.arr

        curs = 0
        for saddle in self.saddles:
            lineprint('\r\tCurrent saddle {0} / {1}'.format(
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
                points, endcoord = trace_gradient_line(coords[0] + 0.1*diff[0],
                                                      coords[1] + 0.1*diff[1],
                                                      self.dx, self.dy,
                                                      self.xnum, self.ynum,
                                                      self.func, self.crits_dict,
                                                      self.start_point,
                                                      direction, self.to_edges)

                # Check here whether the line has gone in totally the
                # wrong direction. Try a simple retest if so.
                if endcoord is not None and False:
                    origangle = n.arctan2(coords[1]-saddley, coords[0]-saddlex)
                    endangle = n.arctan2(endcoord[1]-saddley,
                                         endcoord[0]-saddlex)
                    angle_diff = (endangle % (2*n.pi)) - (origangle % (2*n.pi))
                    if angle_diff > n.pi:
                        angle_diff -= n.pi
                    if angle_diff > n.pi/2:
                        print 'resampling line', origangle, endangle,'from', coords
                        points, endcoord = trace_gradient_line(
                            coords[0] + 0.5*diff[0], coords[1] + 0.5*diff[1],
                            self.dx, self.dy,
                            self.xnum, self.ynum, self.func,
                            self.crits_dict, self.start_point,
                            direction, self.to_edges)

                if len(points) > 4:
                    points = points[4:]
                else:
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
                    fval = self.func(self.sx + endcoord[0]*self.dx,
                                     self.sy + endcoord[1]*self.dy)
                    if fval > 0:
                        self.maxima.append(tuple(endcoord))
                    else:
                        self.minima.append(tuple(endcoord))

        lineprint()
        self.crits_dict = critical_points_to_index_dict(self.crits)
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

    def make_hessian_array(self):
        '''
        Calculate the hessian at every point of self.arr.
        '''
        lineprint('Filling Hessian domain array...')
        arr = self.hessian_arr
        sx, sy = self.start_point
        dx, dy = self.dr
        xnum, ynum = self.shape
        for x in range(xnum):
            lineprint('\r\tx = {0} / {1}'.format(x, xnum), False)
            for y in range(ynum):
                arr[x, y] = hessian_det(self.func, sx + x*dx, sy + y*dy, dx, dy)

        self.hessian_filled = True
        lineprint()

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
            neuline = NeumannLine(start, end, start_crit_type, end_crit_type, line)
            self.graph.add_or_edit_node(start, start_crit_type, neuline)
            self.graph.add_or_edit_node(end, end_crit_type, neuline.inverse())

        self.graph_built = True

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
        min2 = n.sum(mindegs == 2)
        min3 = n.sum(mindegs == 3)
        min4 = n.sum(mindegs == 4)
        min5 = n.sum(mindegs == 5)
        min6 = n.sum(mindegs == 6)
        min7 = n.sum(mindegs == 7)
        min8 = n.sum(mindegs == 8)

        totdegrees = float(max2 + max3 + max4 + max5 + max6 + max7 + max8 +
                           min2 + min3 + min4 + min5 + min6 + min7 + min8)

        frac2 = float(max2 + min2) / totdegrees
        frac3 = float(max3 + min3) / totdegrees
        frac4 = float(max4 + min4) / totdegrees
        frac5 = float(max5 + min5) / totdegrees
        frac6 = float(max6 + min6) / totdegrees
        frac7 = float(max7 + min7) / totdegrees
        frac8 = float(max8 + min8) / totdegrees

        return n.array(zip(range(2, 9),
                           (frac2, frac3, frac4, frac5, frac6, frac7, frac8)))

    def get_domain_areas(self):
        if not self.graph_built:
            self.build_graph()
        return self.graph.get_domain_areas()

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
    def plot(self, trace_lines=True, plot_hessian=False,
             show_saddle_directions=False,
             show_domain_patches=False,
             print_patch_areas=False,
             figsize=None,
             show_sample_directions=False,
             save=False):
        '''
        Plot and return a graph showing (optionally):
        - Neumann lines
        - Hessian domains
        - Hessian eigenvectors at detected saddles
        - Coloured patches representing each closed Neumann domain
        '''
        if save:
            mpl_interactive(False)
        else:
            mpl_interactive(True)

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

        fig, ax = plt.subplots()

        if not show_domain_patches:
            ax.imshow(plotarr, cmap='RdYlBu_r', interpolation='none', alpha=0.6)
        else:
            ax.imshow(plotarr, cmap='RdYlBu_r', interpolation='none', alpha=0.4)

        ax.set_xlim(0, plotarr.shape[0])
        ax.set_ylim(0, plotarr.shape[1])

        ax.set_xticks([])
        ax.set_yticks([])

        if show_domain_patches:
            for domain in self.domains:
                colour = hsv_to_rgb(n.random.random(), 1., 1.)
                ps = domain.as_closed_curves()
                for p in ps:
                    if len(ps) > 1:
                        patch = Polygon(p, alpha=1.0,
                                        closed=True,
                                        color=colour,
                                        linestyle=n.random.choice(
                                            patch_linestyles),
                                        linewidth=2,
                                        )
                    else:
                        patch = Polygon(p, alpha=0.7,
                                        closed=True,
                                        color=colour,
                                        linestyle=n.random.choice(
                                            patch_linestyles),
                                        linewidth=2,
                                        )
                    ax.add_patch(patch)
                if print_patch_areas:
                    area = domain.area()
                    pos = n.average(p, axis=0)
                    ax.text(pos[0], pos[1],'{:.1f}'.format(area))

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

        if show_domain_patches:
            ax.contour(plotarr, levels=[0], alpha=0.3, linestyles=['--'])
        else:
            ax.contour(plotarr, levels=[0], alpha=0.2)

        if trace_lines:
            for line in self.lines:
                segs = sanitise_line(line)
                # if len(segs) > 0:
                #     print len(segs), map(len, segs),'segs are', segs
                for seg in segs: #filter(lambda j: len(j)==1, segs):
                    ax.plot(seg[:, 0], seg[:, 1],'-', color='purple')
                #ax.plot(line[:, 0], line[:, 1],'-', color='purple')

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

        self.figax = (fig, ax)


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
    def plot3d(self, clf=True):
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

def get_filled_array(xnum, ynum, dx, dy, func, start_point=(0.0, 0.0)):
    arr = n.zeros((xnum, ynum), dtype=n.float64)
    sx, sy = start_point
    for x in range(xnum):
        for y in range(ynum):
            arr[x, y] = func(sx + x*dx, sy + y*dy)
    return arr

def trace_gradient_line(sx, sy, dx, dy, xnum, ynum, func,
                        critdict, start_point, direction, to_edges):
    #print 'Tracing gradient line'
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
        #     cross, wx, wy = doVectorsCross(p, dp, q, dq)
        #     print cross, wx, wy
        #     if cross:
        #         cx -= 0.5*n.cos(angle)
        #         cy -= 0.5*n.sin(angle)
        points.append([cx, cy])

        if cx < 0 or cx > xnum or cy < 0 or cy > ynum:
            if to_edges in ['periodic','fourier']:
                   # Need extra condition with fourier to get the signs right
                cx %= xnum
                cy %= ynum
            else:
                return (points, None)

        nearx, neary = int(n.round(cx)), int(n.round(cy))
        if to_edges in ['periodic','fourier']:
                   # Need extra condition with fourier to get the signs right
            nearx %= xnum
            neary %= ynum
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
#even_adj_indices = n.array([[-1, -1], [-1, 0], [-1, 1],
#                   [0, 1], [1, 1], [1, 0], [1, -1], [0, -1]])
#odd_adj_indices = even_adj_indices
def get_critical_points(arr, to_edges=False):
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
            elif to_edges == 'fourier':
                off_x = n.logical_or(ais[:, 0]<0, ais[:, 0]>=lx)
                off_y = n.logical_or(ais[:, 1]<0, ais[:, 1]>=ly)
                off_edge = n.logical_or(off_x, off_y)
                border_mult = ((-1.0*off_edge) + 0.5)*2.0
                ais[:, 0] = ais[:, 0] % lx
                ais[:, 1] = ais[:, 1] % ly


            for i in range(6):
                adjs[i] = arr[ais[i][0] % lx, ais[i][1] % ly]
            if to_edges == 'fourier':
                adjs *= border_mult

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
                print 'A failure occurred at', x, y

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
        return 'fail'

def plot_arr_with_crits(arr, crits):
    maxima, minima, saddles, degenerate = crits
    maxima = n.array(maxima)
    minima = n.array(minima)
    saddles = n.array(saddles)
    degenerate = n.array(degenerate)

    fig, ax = plt.subplots()

    ax.imshow(arr, cmap='RdYlBu_r', interpolation='none', alpha=0.6)

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
    if seed == 0:
        seed = n.random.randint(10000000000)
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

def periodic_random_wave_function(number=50, scale=5, seed=0, returnall=False):
    if seed == 0:
        seed = n.random.randint(10000000000)
    generator = n.random.RandomState()
    generator.seed(seed)

    possible_wvs = duofactors(scale)
    possible_signs = map(n.array, [[1, 1], [1, -1], [-1, 1], [-1, -1]])

    amps = generator.normal(size=number)
    phases = 2*n.pi*generator.rand(number)
    wvs = n.zeros((number, 2), dtype=n.float64)
    for i in range(number):
        wv = n.array(possible_wvs[generator.randint(len(possible_wvs))],
                     dtype=n.float64)
        generator.shuffle(wv)
        wv *= possible_signs[generator.randint(len(possible_signs))]
        wv *= 2*n.pi/n.sqrt(scale/2.)
        wvs[i] = wv

    def func(x, y):
        res = 0.0
        xyarr = n.array([x, y], dtype=n.float64)
        interior = wvs.dot(xyarr) + phases
        exterior = amps*n.sin(interior)
        return n.sum(exterior)
    if returnall:
        return func, (amps, wvs, phases)
    return func



def mag(v):
    return n.sqrt(v.dot(v))
    

import sys
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
            returns.append((-1.0, (n.argmax(region_4) + changes[3]) % len(vals)))
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
    return n.abs(area/2.)


def sanitise_domain(ps, shape=None):
    '''Takes a domain that may pass through a boundary, and sanitises
    it to have the right shape.'''
    if shape is None:
        dx, dy = None, None
    else:
        dx, dy = shape

    for i in range(len(ps)):
        cur = ps[i]
        nex = ps[(i+1) % len(ps)]
        if n.abs(nex[0]-cur[0]) > 5:
            diff = n.round(nex[0] - cur[0])
            ps[i+1:, 0] -= diff
        if n.abs(nex[1]-cur[1]) > 5:
            diff = n.round(nex[1] - cur[1])
            ps[i+1:, 0] -= diff
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
            intersect, t, u = doVectorsCross(p, r, q, s)
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


def doVectorsCross(p, r, q, s):
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


def animate(filen='test', path='onephase', number=200):
    f, d2 = random_wave_function(50, 10, returnall=True)

    amps, wvs, phases = d2

    for i in range(number):
        if path == 'allphase':
            phases += 2*n.pi/number
        elif path[:9] == 'manyphase':
            val = int(path[9:])
            phases[:val] += 2*n.pi/number
        else:
            phases[0] += 2*n.pi/number

        def func(x, y):
            res = 0.0
            xyarr = n.array([x, y], dtype=n.float64)
            interior = wvs.dot(xyarr) + phases
            exterior = amps*n.sin(interior)
            return n.sum(exterior)
        a = NeumannTracer(300, 300, n.pi/190, n.pi/190, f)
        a.plot(save='{0}_{1:04}.png'.format(filen, i))
        del a  # For some reason, the reference count builds up?


def periodic_animate(scale=5, number=50, frames=200, downscale=2):
    length = int(100/float(downscale) * float(scale)/5.)
    periodicity = n.sqrt(scale/2.)

    f, d2 = periodic_random_wave_function(number, scale, returnall=True)
    amps, wvs, phases = d2

    for i in range(frames):
        phases[0] += 2*n.pi/frames

        def func(x, y):
            res = 0.0
            xyarr = n.array([x, y], dtype=n.float64)
            interior = wvs.dot(xyarr) + phases
            exterior = amps*n.sin(interior)
            return n.sum(exterior)
        a = NeumannTracer(length, length,
                          periodicity/float(length), periodicity/float(length),
                          func, to_edges=True)
        a.plot(save='panim1_{0:04}.png'.format(i))


def get_periodic_tracer(scale=5, number=50, downscale=2, returnall=False):
    length = int(100/float(downscale) * float(scale)/5.)
    periodicity = n.sqrt(scale/2.)
    f, d2 = periodic_random_wave_function(number, scale, returnall=True)
    tracer = NeumannTracer(length, length,
                           periodicity/(1.*length), periodicity/(1.*length),
                           f,
                           to_edges='periodic')
    if returnall:
        return tracer, f, d2
    return tracer
