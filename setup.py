from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
        Extension("cneumann",["cneumann.pyx"], libraries=["m"]),
#        Extension("generationclassesfail",["generationclassesfail.pyx"],libraries=["m"]),
        ]

setup(
  name = 'cneumann',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
