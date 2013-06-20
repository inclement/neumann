import numpy as n
from itertools import product

import matplotlib.pyplot as plt

class NeumannTracer(object):
    def __init__(self,xnum,ynum,dx,dy,func,start_point=(0.0,0.0),to_edges=False):
        self.arr = n.zeros((xnum,ynum),dtype=n.float64)
        self.xnum = xnum
        self.ynum = ynum
        self.shape = (xnum,ynum)
        self.dx = dx
        self.dy = dy
        self.dr = (dx,dy)
        self.func = func
        self.start_point = start_point
        self.sx = start_point[0]
        self.sy = start_point[1]

        self.to_edges = to_edges

        self.arr_filled = False
        self.found_crits = False
        self.traced_lines = False

        self.lines = []

        self.fill_arr()

        self.figax = (None,None)

    def fill_arr(self):
        print 'Populating function sample array...'
        arr = self.arr
        sx,sy = self.start_point
        dx,dy = self.dr
        xnum,ynum = self.shape
        for x in range(xnum):
            lineprint('\r\tx = {0} / {1}'.format(x,xnum),False)
            for y in range(ynum):
                arr[x,y] = self.func(sx + x*dx, sy + y*dy)
        self.arr_filled = True

    def find_critical_points(self):
        if not self.arr_filled:
            self.fill_arr()
        print 'Finding critical points...'
        maxima,minima,saddles,degenerate = get_critical_points(self.arr,self.to_edges)
        self.crits = (maxima,minima,saddles,degenerate)
        self.prune_critical_points()
        self.maxima = self.crits[0]
        self.minima = self.crits[1]
        self.saddles = self.crits[2]
        self.degenerate = self.crits[3]

        self.crits_dict = critical_points_to_index_dict(self.crits)

        self.found_crits = True

    def prune_critical_points(self):
        print 'Pruning critical points'
        maxima,minima,saddles,degenerate = self.crits

        tupmaxima = [tuple(j) for j in maxima]
        tupminima = [tuple(j) for j in minima]

        realmaxima = n.ones(len(maxima),dtype=bool)
        realminima = n.ones(len(minima),dtype=bool)
        realsaddles = n.ones(len(saddles),dtype=bool)

        for i in range(len(saddles)):
            saddle = saddles[i]
            adj = all_adj_indices.copy()
            adj[:,0] += saddle[0]
            adj[:,1] += saddle[1]
            adj[:,0] = adj[:,0] % self.xnum
            adj[:,1] = adj[:,1] % self.ynum
            adjmax = []
            adjmin = []
            for coords in adj:
                coords = tuple(coords)
                if tuple(coords) in tupmaxima:
                    #print saddle,coords,tupmaxima
                    adjmax.append(coords)
                if tuple(coords) in tupminima:
                    adjmin.append(coords)
                    #print saddle,coords,tupminima
            if len(adjmax) > 1:
                heights = map(lambda j: self.arr[j[0],j[1]],adjmax)
                fakemax = n.argmin(heights)
                realmaxima[tupmaxima.index(adjmax[fakemax])] = False
                realsaddles[i] = False
            elif len(adjmin) > 1:
                heights = map(lambda j: self.arr[j[0],j[1]],adjmin)
                fakemin = n.argmax(heights)
                realminima[tupminima.index(adjmin[fakemin])] = False
                realsaddles[i] = False
        maxima = maxima[realmaxima]
        minima = minima[realminima]
        saddles = saddles[realsaddles]
        self.crits = maxima,minima,saddles,degenerate
    def trace_neumann_lines(self):
        if not self.arr_filled:
            self.fill_arr()
        if not self.found_crits:
            self.find_critical_points()

        print 'Tracing Neumann lines...'

        arr = self.arr

        curs = 0
        for saddle in self.saddles:
            lineprint('\r\tCurrent saddle {0} / {1}'.format(curs,len(self.saddles)),False)
            curs += 1

            saddlex,saddley = saddle
            if saddlex % 2 == 0:
                ais = even_adj_indices.copy()
            else:
                ais = odd_adj_indices.copy()

            ais[:,0] += saddlex
            ais[:,1] += saddley

            val = arr[saddlex,saddley]
            adjs = n.zeros(6,dtype=n.float64)
            for i in range(6):
                adjs[i] = arr[ais[i][0] % self.xnum, ais[i][1] % self.ynum]

            adjs = adjs - val

            tracing_start_points = []
            for i in range(6):
                cur_adj = adjs[i]
                next_adj = adjs[(i+1) % 6]
                if n.sign(next_adj) != n.sign(cur_adj):
                    sign = n.sign(cur_adj)
                    tracing_start_points.append((sign,ais[i]))

            for coords in tracing_start_points:
                sign,coords = coords
                if sign == 1.0:
                    direction = 'down'
                else:
                    direction = 'up'
                points,endcoord = trace_gradient_line(coords[0],coords[1],self.dx,self.dy,self.xnum,self.ynum,self.func,self.crits_dict,self.start_point,direction,self.to_edges)
                points = [saddle] + points
                self.lines.append(n.array(points))

        lineprint()
        self.traced_lines = True

    def print_critical_heights(self):
        if not self.arr_filled:
            self.fill_arr()
        if not self.found_crits:
            self.find_critical_points()
        
        print 'maxima'
        for entry in self.maxima:
            print entry,self.func(self.sx+entry[0]*self.dx,self.sy+entry[1]*self.dy),self.arr[entry[0],entry[1]]
        print 'minima'
        for entry in self.minima:
            print entry,self.func(self.sx+entry[0]*self.dx,self.sy+entry[1]*self.dy),self.arr[entry[0],entry[1]]
        print 'saddles'
        for entry in self.saddles:
            print entry,self.func(self.sx+entry[0]*self.dx,self.sy+entry[1]*self.dy),self.arr[entry[0],entry[1]]
        print 'degenerate'
        for entry in self.degenerate:
            print entry,self.func(self.sx+entry[0]*self.dx,self.sy+entry[1]*self.dy),self.arr[entry[0],entry[1]]

    def plot(self,trace_lines=True):
        if not self.arr_filled:
            self.fill_arr()
        if not self.found_crits:
            self.find_critical_points()
        if not self.traced_lines and trace_lines:
            self.trace_neumann_lines()

        plotarr = n.rot90(self.arr[::-1],3)
        
        fig,ax = plot_arr_with_crits(plotarr,self.crits)

        for line in self.lines:
            segs = sanitise_line(line)
            for seg in segs:
                ax.plot(seg[:,0],seg[:,1],'-',color='purple')
            ax.plot(line[:,0],line[:,1],'-',color='purple')

        ax.contour(plotarr,levels=[0],alpha=0.2)
        
        self.figax = (fig,ax)
        return fig,ax

def get_filled_array(xnum,ynum,dx,dy,func,start_point=(0.0,0.0)):
    arr = n.zeros((xnum,ynum),dtype=n.float64)
    sx,sy = start_point
    for x in range(xnum):
        for y in range(ynum):
            arr[x,y] = func(sx + x*dx, sy + y*dy)
    return arr

def trace_gradient_line(sx,sy,dx,dy,xnum,ynum,func,critdict,start_point,direction,to_edges):
    #print 'Tracing gradient line'
    cx,cy = sx,sy
    startx,starty = start_point
    if direction == 'down':
        dirfac = -1
    else:
        dirfac = 1

    points = [[cx,cy]]

    while True:
        gradient = grad(func, startx+cx*dx,starty+cy*dy,dx,dy) * dirfac
        angle = n.arctan2(gradient[1],gradient[0])

        cx += 0.25*n.cos(angle)
        cy += 0.25*n.sin(angle)

        if len(points)>10:
            if mag(n.array([cx,cy])-n.array(points[-10])) < 1.0:
                return (points,None)
        
        points.append([cx,cy])

        if cx < 0 or cx > xnum or cy < 0 or cy > ynum:
            if to_edges in ['periodic','fourier']: # Need extra condition with fourier to get the signs right
                cx %= xnum
                cy %= ynum
            else:    
                return (points,None)

        nearx,neary = int(n.round(cx)),int(n.round(cy))
        nearx %= xnum
        neary %= ynum
        if (nearx,neary) in critdict:
            points.append((nearx,neary))
            return (points,(nearx,neary))

def grad(func,x,y,dx,dy):
    dfdx = (func(x,y)-func(x+0.05*dx,y))/(0.05*dx)
    dfdy = (func(x,y)-func(x,y+0.05*dy))/(0.05*dy)
    return n.array([dfdx,dfdy])

# def hessian(func,x,y,dx,dy):
#     dfdx,dfdy = grad(func,x,y,dx,dy)

#     dfdxdx = (grad(func,x+0.05*dx,y,dx,dy)[0] - dfdx) / (0.05*dx)
#     dfdydy = (grad(func,x2
    
def critical_points_to_index_dict(crits):
    maxima,minima,saddles,degenerate = crits
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

even_adj_indices = n.array([[-1,-1],[-1,0],[0,1],[1,0],[1,-1],[0,-1]])
odd_adj_indices =  n.array([[-1,0],[-1,1],[0,1],[1,1],[1,0],[0,-1]])
all_adj_indices = n.array([[-1,-1],[-1,0],[-1,1],[0,1],[1,1],[1,0],[1,-1],[0,-1]])
#even_adj_indices = n.array([[-1,-1],[-1,0],[-1,1],[0,1],[1,1],[1,0],[1,-1],[0,-1]])
#odd_adj_indices = even_adj_indices
def get_critical_points(arr,to_edges=False):
    lx,ly = arr.shape
    adjs = n.zeros(6,dtype=n.float64)

    maxima = []
    minima = []
    saddles = []
    degenerate = []

    border_mult = n.ones(6,dtype=n.float64)

    for x,y in product(xrange(lx),xrange(ly)):
        if to_edges or (x != 0 and y != 0 and x != (lx-1) and y != (ly-1)):
            val = arr[x,y]
            if x % 2 == 0:
                ais = even_adj_indices.copy()
            else:
                ais = odd_adj_indices.copy()

            ais[:,0] += x
            ais[:,1] += y

            if to_edges == 'periodic':
                ais[:,0] = ais[:,0] % lx
                ais[:,1] = ais[:,1] % ly
            elif to_edges == 'fourier':
                off_x = n.logical_or(ais[:,0]<0,ais[:,0]>=lx)
                off_y = n.logical_or(ais[:,1]<0,ais[:,1]>=ly)
                off_edge = n.logical_or(off_x,off_y)
                border_mult = ((-1.0*off_edge) + 0.5)*2.0
                ais[:,0] = ais[:,0] % lx
                ais[:,1] = ais[:,1] % ly


            for i in range(6):
                adjs[i] = arr[ais[i][0] % lx, ais[i][1] % ly]
            if to_edges == 'fourier':
                adjs *= border_mult

            point_type = classify_point(adjs-val)

            if point_type == 'maximum':
                maxima.append((x,y))
            elif point_type == 'minimum':
                minima.append((x,y))
            elif point_type == 'saddle':
                saddles.append((x,y))
            elif point_type == 'degenerate':
                degenerate.append((x,y))
            elif point_type == 'fail':
                print 'A failure occurred at',x,y


    return (n.array(maxima),n.array(minima),n.array(saddles),n.array(degenerate))

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

def plot_arr_with_crits(arr,crits):
    maxima,minima,saddles,degenerate = crits

    fig,ax = plt.subplots()

    ax.imshow(arr,cmap='RdYlBu',interpolation='none',alpha=0.6)

    legend_entries = []

    if len(maxima) > 0:
        ax.scatter(maxima[:,0],maxima[:,1],60,c='r')
        legend_entries.append('maxima')
    if len(minima) > 0:
        ax.scatter(minima[:,0],minima[:,1],60,c='b')
        legend_entries.append('minima')
    if len(saddles) > 0:
        ax.scatter(saddles[:,0],saddles[:,1],60,color='yellow')
        legend_entries.append('saddles')
    if len(degenerate) > 0:
        ax.scatter(degenerate[:,0],degenerate[:,1],c='orange')
        legend_entries.append('degenerate')

    return fig,ax
        

def sanitise_line(l):
    curcut = 0
    segs = []
    for i in range(len(l)-1):
        next = l[i+1]
        cur = l[i]
        #print n.abs(next[0]-cur[0])
        #print n.abs(next[0]-cur[0]),n.abs(next[1]-cur[1])
        if n.abs(next[0]-cur[0]) > 5 or n.abs(next[1]-cur[1]) > 5:
            #print 'bigchange',len(segs)
            segs.append(l[curcut:(i+1)])
            #print '->',len(segs)
            curcut = i+1
    if curcut < len(l):
        segs.append(l[curcut:])
    return segs
        

def random_wave_function(number=50,wvmag=5,seed=0):
    if seed == 0:
        seed = n.random.randint(10000000000)
    generator = n.random.RandomState()
    generator.seed(seed)

    amps = generator.normal(size=number)
    phases = generator.rand(number)
    wvs = n.zeros((number,2),dtype=n.float64)
    for i in range(number):
        wv = generator.normal(size=2)
        wv /= mag(wv)
        wv *= wvmag
        wvs[i] = wv

    def func(x,y):
        res = 0.0
        for i in range(number):
            res += amps[i] * n.sin(2*n.pi/wvmag * wvs[i].dot(n.array([x,y])) + phases[i])
        return res

    return func

def mag(v):
    return n.sqrt(v.dot(v))

import sys
def lineprint(s='',newline=True):
    sys.stdout.write(s)
    if newline:
        sys.stdout.write('\n')
    sys.stdout.flush()

