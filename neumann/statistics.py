import neumann as neu
import numpy as n


from os.path import exists
def get_filen(base, extension):
    i = 1
    while exists('{}_{}.{}'.format(base, i, extension)):
        i += 1

    return '{}_{}.{}'.format(base, i, extension)


def save_domain_statistics_at(scales, domains=1000, downscale=3):
    '''Runs get_statistics_at with the same arguments, and saves the
    results in a file starting with filen.'''
    import cPickle
    import json
    results = get_domain_statistics_at(scales, domains, downscale)
    i = 1

    for key, value in results.iteritems():

        length_vactor = value['length_factor']
        for quantity in ['areas', 'perimeters', 'rhos', 'diameters',
                         'rhos_lens', 'rhos_wedge', 'rhos_star']:
            values = value[quantity]
            base_filen = 'neumann-energy{}-{}-lengthfactor{}-{}results'.format(
                key, quantity, value['length_factor'], len(values))
            filen = get_filen(base_filen, 'txt')
            with open(filen, 'w') as fileh:
                for v in values:
                    fileh.write('{}\n'.format(v))
    return results
    

def get_domain_statistics_at(scales, domains=1000, downscale=3):
    '''For every scale in scales, generates random torus eigenfunctions at
    that scale, counting the domains unitl it has found at least the
    number specified as an argument.

    downscale is a resolution parameter, it can usually be ignored.
    '''
    results = {}

    for scale in scales:
        print 'Getting statistics at scale', scale
        areas = []
        perimeters = []
        rhos = []
        rhos_lens = []
        rhos_wedge = []
        rhos_star = []
        diameters = []
        degree_dists = n.zeros((8, 2))

        while len(areas) < domains:
            print 'Currently done', len(areas), 'domains'
            a = neu.get_periodic_tracer(scale, downscale=downscale)
            if a is None:
                break  # None is returned if there are no compatible
                       # periodic wavevectors
            areas.extend(a.get_domain_areas())
            perimeters.extend(a.get_domain_perimeters())

            rhos_by_type = a.get_domain_rhos_by_type()
            rhos_lens.extend(rhos_by_type['lens'])
            rhos_wedge.extend(rhos_by_type['wedge'])
            rhos_star.extend(rhos_by_type['star'])
            
            rhos.extend(a.get_domain_rhos())
            diameters.extend(a.get_domain_diameters())
            degree_dists += a.get_critical_degree_dists()

        length = float(int(100/float(downscale) * float(scale)/5.))
        real_length = 2*n.pi
        length_factor = 1/length * real_length
        wavelength = 1/n.sqrt(2)

        results[scale] = {'areas': areas,
                          'perimeters': perimeters,
                          'rhos': rhos,
                          'rhos_lens': rhos_lens,
                          'rhos_wedge': rhos_wedge,
                          'rhos_star': rhos_star,
                          'diameters': diameters,
                          'degrees': degree_dists,
                          'length_factor': length_factor}

    for value in results.itervalues():
        value['degrees'] /= degree_dists[0][0]/2.

    return results


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description='Generate statistics of Neumann domains')

    parser.add_argument('energy', help='The eigenvalue to sample at', type=int)
    parser.add_argument('domains', help='The number of domains to sample',
                        type=int)
    parser.add_argument('--downscale', help='The resolution downsampling to use',
                        type=float, default=2)
    parser.add_argument('--num_per_file', help='The number of domains to find before saving',
                        type=int, default=100000)

    args = parser.parse_args()

    # stats = get_domain_statistics_at([parser.scales],
    #                                  parser.domains,
    #                                  parser.downscale)

    num_per_file = args.num_per_file
    if num_per_file > args.domains:
        num_per_file = args.domains
    
    for i in range(0, args.domains, num_per_file):
        save_domain_statistics_at([args.energy], num_per_file,
                                  args.downscale)
                              

