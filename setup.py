def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('',parent_package,top_path)
    config.add_data_dir(('ImpactAtmosphere/data','ImpactAtmosphere/data'))
    return config

from numpy.distutils.core import setup, Extension

sources = ['src/HCN_transport.f90', 'src/EvolveAtmFort.f90', 'src/minpack/dpmpar.f', 'src/minpack/enorm.f',
           'src/minpack/fdjac2.f', 'src/minpack/lmdif.f', 'src/minpack/lmdif1.f',
           'src/minpack/lmdif2.f', 'src/minpack/lmpar.f', 'src/minpack/qrfac.f', 'src/minpack/qrsolv.f']

only = 'only: setup rhs rhs_verbose hcn_transport :'

extensions = [Extension(name="EvolveAtmFort",
                        sources=sources,
                        libraries=['lapack'],
                        f2py_options=only.split()),]

setup(name = 'ImpactAtmosphere',
      python_requires='>3.6.0',
      packages=['ImpactAtmosphere'],
      version='4.1.2',
      ext_modules=extensions,
      configuration=configuration)
