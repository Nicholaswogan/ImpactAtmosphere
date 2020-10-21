from numpy.distutils.core import setup, Extension

extensions = [
    Extension(name="photochem",
              sources=["wrap_no_comments.f"]),
]

setup(name = 'EvolveAtm',
packages=['EvolveAtm'],
version='2.2',
ext_modules=extensions)
