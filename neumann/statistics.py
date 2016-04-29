import neumann as neu

def save_domain_statistics_at(scales, domains=1000, downscale=3, filen='neumann_results'):
    '''Runs get_statistics_at with the same arguments, and saves the
    results in a file starting with filen.'''
    import cPickle
    import json
    results = get_statistics_at(scales, domains, downscale)
    i = 1

    while os.path.exists('{}_{}.pickle'.format(filen, i)):
        i += 1
    with open('{}_{}.pickle'.format(filen, i), 'wb') as fileh:
        cPickle.dump(results, fileh)

    results2 = {}
    for key, value in results.iteritems():
        results2[key] = (value[0], value[1],
                         value[2], map(list, list(value[3])), value[4])
    with open('{}_{}.json'.format(filen, i), 'w') as fileh:
        json.dump(results2, fileh)
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
        degree_dists = n.zeros((8,2))

        while len(areas) < domains:
            print 'Currently done', len(areas), 'domains'
            a = get_periodic_tracer(scale, downscale=downscale)
            if a is None:
                break  # None is returned if there are no compatible
                       # periodic wavevectors
            areas.extend(a.get_domain_areas())
            perimeters.extend(a.get_domain_perimeters())
            rhos.extend(a.get_domain_rhos())
            degree_dists += a.get_critical_degree_dists()

        length = float(int(100/float(downscale) * float(scale)/5.))
        real_length = 2*n.pi
        length_factor = 1/length * real_length
        wavelength = 1/n.sqrt(2)

        results[scale] = (areas, perimeters, rhos, degree_dists, length_factor)

    for areas, perimeters, rhos, degrees, length_factor in results.itervalues():
        degrees /= degrees[0][0]/2.

    return results
