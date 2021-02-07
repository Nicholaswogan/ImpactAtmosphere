from numpy.distutils.core import setup, Extension

sources = ['src/diffusion.f90', 'src/EvolveAtmFort.f90', 'src/minpack/dpmpar.f', 'src/minpack/enorm.f',
           'src/minpack/fdjac2.f', 'src/minpack/lmdif.f', 'src/minpack/lmdif1.f',
           'src/minpack/lmdif2.f', 'src/minpack/lmpar.f', 'src/minpack/qrfac.f', 'src/minpack/qrsolv.f']

only = 'only: setup rhs rhs_verbose diffuse :'

extensions = [Extension(name="EvolveAtmFort",
                        sources=sources,
                        libraries=['lapack'],
                        f2py_options=only.split()),]

setup(name = 'ImpactAtmosphere',
      python_requires='>3.6.0',
      packages=['ImpactAtmosphere'],
      version='4.0.1',
      ext_modules=extensions)
