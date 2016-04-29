from setuptools import setup, find_packages
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

ext_modules = [
        Extension("neumann.cneumann",["neumann/cneumann.pyx"],
                  libraries=["m"]),
        Extension("neumann.ctwosphere",["neumann/ctwosphere.pyx"],
                  libraries=["m"]),
        ]

setup(
    name = 'NeumannTracer',
    description=('Numerical routines for locating and analysing '
                 'Neumann domains in real scalar functions.'),
    author='Alexander Taylor',
    author_email='alexander.taylor@bristol.ac.uk',
    cmdclass = {'build_ext': build_ext},
    include_dirs = [numpy.get_include()],
    ext_modules = ext_modules
)
