from distutils.core import setup, Extension
from numarray.numarrayext import NumarrayExtension
import sys

if not hasattr(sys, 'version_info') or sys.version_info < (2,2,0,'alpha',0):
    raise SysError

setup(name = "raven",
    version = "0.1",
    description = "A tool for manipulating ENZO data",
    url = "http://www.stanford.edu/~mturk/raven.html",
    author="Matthew Turk",
    author_email="mturk@stanford.edu",
    packages=["raven"],
    license="GPL-2",
    package_dir={"":""},
    ext_modules=[NumarrayExtension("raven.RavenCombine",['raven/point_combine.c'],\
        include_dirs=["./"],
        library_dirs=["./"],
        libraries=['m'])],
    py_modules=["raven"])
