from neumann.ctwosphere import spherical_harmonic_2d
from numpy.random import RandomState
import numpy as n
from functools import partial

from scipy.special import sph_harm

def get_random_coefficients(l, seed=0):
    if seed == 0:
        seed = n.random.randint(10000000)
    generator = RandomState()
    generator.seed(seed)
    
    coeffs = []
    for m in range(-l, l+1):
        coeffs.append(generator.normal())
    return l, coeffs

def random_spherical_harmonic(l, coeffs, theta, phi):
    result = 0.0
    for index, m in enumerate(range(-l, l+1)):
        result += coeffs[index] * spherical_harmonic_2d(l, m, theta, phi)
    return result

def random_real_spherical_harmonic(l, coeffs, theta, phi):
    result = 0.0
    for index, m in enumerate(range(-l, l+1)):
        if m < 0:
            continue
        if m == 0:
            result += coeffs[index] * spherical_harmonic_2d(l, m, theta, phi).real
        else:
            result += coeffs[index] * spherical_harmonic_2d(l, m, theta, phi).real
            result += coeffs[index] * spherical_harmonic_2d(l, m, theta, phi).imag
    return result
        
def get_random_spherical_harmonic(l):
    l, coeffs = get_random_coefficients(l)
    func = partial(random_spherical_harmonic, l, coeffs)
    return func

def get_random_real_spherical_harmonic(l):
    l, coeffs = get_random_coefficients(l)
    func = partial(random_real_spherical_harmonic, l, coeffs)
    return func

def plot_func(func, modulation=0., shape=(101, 101), cmap='RdYlBu',
              emph_doms=False, clf=True, invert_modulation=True,
              stereographic=False, offset=(0, 0, 0)):
    '''Takes a function and plots it over a sphere.

    Modulation argument controls surface modulation
    '''

    import mayavi.mlab as may
    if clf:
        may.clf()

    if not stereographic:
        theta, phi = n.mgrid[0:n.pi:shape[0]*1j, 0:2*n.pi:shape[1]*1j]
    else:
        theta, phi = n.mgrid[0:0.93*n.pi:shape[0]*1j, 0:2*n.pi:shape[1]*1j]

    func = n.vectorize(func)
    s = func(theta, phi).real

    maxs = n.max(s)
    norms = s / maxs
    if invert_modulation:
        r = 1. - modulation*norms
    else:
        r = 1. + modulation * norms

    if not stereographic:
        x = r*n.sin(theta)*n.cos(phi)
        y = r*n.sin(theta)*n.sin(phi)
        z = r*n.cos(theta)
    else:
        x = 2*r*n.tan(theta/2.)*n.cos(phi)
        y = 2*r*n.tan(theta/2.)*n.sin(phi)
        z = n.zeros(theta.shape) + r
        
    x += offset[0]
    y += offset[1]
    z += offset[2]
    print('xmin is', n.min(x))

    surf = may.mesh(x, y, z, scalars=s, colormap=cmap)
    if emph_doms:
        surf = may.mesh(x, y, z, scalars=s, colormap=cmap)
        lut = surf.module_manager.scalar_lut_manager.lut.table.to_array()
        lut[128:, 0] = n.zeros(128)
        lut[128:, 1] = n.zeros(128)
        lut[128:, 2] = n.zeros(128)
        lut[:128, 0] = n.ones(128)*255
        lut[:128, 1] = n.ones(128)*255
        lut[:128, 2] = n.ones(128)*255
        print('lut is', lut)
        surf.module_manager.scalar_lut_manager.lut.table = lut
        may.draw()

def plot_func_vispy(func, modulation=0., shape=(101, 101), cmap='RdYlBu',
                    emph_doms=False, clf=True, invert_modulation=True,
                    stereographic=False, offset=(0, 0, 0)):
    '''Takes a function and plots it over a sphere.

    Modulation argument controls surface modulation
    '''

    import mayavi.mlab as may
    if clf:
        may.clf()

    if not stereographic:
        theta, phi = n.mgrid[0:n.pi:shape[0]*1j, 0:2*n.pi:shape[1]*1j]
    else:
        theta, phi = n.mgrid[0:0.93*n.pi:shape[0]*1j, 0:2*n.pi:shape[1]*1j]

    func = n.vectorize(func)
    s = func(theta, phi).real

    maxs = n.max(s)
    norms = s / maxs
    if invert_modulation:
        r = 1. - modulation*norms
    else:
        r = 1. + modulation * norms

    if not stereographic:
        x = r*n.sin(theta)*n.cos(phi)
        y = r*n.sin(theta)*n.sin(phi)
        z = r*n.cos(theta)
    else:
        x = 2*r*n.tan(theta/2.)*n.cos(phi)
        y = 2*r*n.tan(theta/2.)*n.sin(phi)
        z = n.zeros(theta.shape) + r
        
    x += offset[0]
    y += offset[1]
    z += offset[2]
    print('xmin is', n.min(x))

    surf = may.mesh(x, y, z, scalars=s, colormap=cmap)
    if emph_doms:
        surf = may.mesh(x, y, z, scalars=s, colormap=cmap)
        lut = surf.module_manager.scalar_lut_manager.lut.table.to_array()
        lut[128:, 0] = n.zeros(128)
        lut[128:, 1] = n.zeros(128)
        lut[128:, 2] = n.zeros(128)
        lut[:128, 0] = n.ones(128)*255
        lut[:128, 1] = n.ones(128)*255
        lut[:128, 2] = n.ones(128)*255
        print('lut is', lut)
        surf.module_manager.scalar_lut_manager.lut.table = lut
        may.draw()

    
def plot_comparison(l, m, shape=(101, 101)):

    import mayavi.mlab as may

    theta, phi = n.mgrid[0:n.pi:shape[0]*1j, 0:2*n.pi:shape[1]*1j]

    r = 1.
    x = r*n.sin(theta)*n.cos(phi)
    y = r*n.sin(theta)*n.sin(phi)
    z = r*n.cos(theta)

    func1 = n.vectorize(partial(spherical_harmonic_2d, l, m))
    s1 = func1(theta, phi).real

    func2 = n.vectorize(partial(sph_harm, m, l))
    s2 = func2(phi, theta).real

    may.mesh(x, y, z, scalars=s1, colormap='RdYlBu')
    may.mesh(x-3, y, z, scalars=s2, colormap='RdYlBu')
