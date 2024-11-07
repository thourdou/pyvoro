# 
# setup.py : pyvoro python interface to voro++
# 
# this extension to voro++ is released under the original modified BSD license
# and constitutes an Extension to the original project.
#
# Copyright (c) Joe Jordan 2012
# contact: <joe.jordan@imperial.ac.uk> or <tehwalrus@h2j9k.org>
#

from distutils.core import setup, Extension
from Cython.Build import cythonize

# fall back to provided cpp file if Cython is not found
extensions = [
    Extension("pyvoro",
              sources=["pyvoro/pyvoro.pyx",
                       "pyvoro/wrapper.cpp",
                       "src/voro++.cc"],
              include_dirs=["src"],
              language="c++",
              )
]
extensions = cythonize(extensions, nthreads=8)

setup(
    name="pyvoro",
    version="2.0.0",
    description="3D Voronoi tessellations: a python entry point for the voro++ library.",
    author="Theophile Hourdou",
    author_email="joe.jordan@imperial.ac.uk",
    url="https://github.com/joe-jordan/pyvoro",
    download_url="https://github.com/joe-jordan/pyvoro/tarball/v1.3.2",
    packages=["pyvoro",],
    package_dir={"pyvoro": "pyvoro"},
    ext_modules=extensions,
    keywords=["geometry", "mathematics", "Voronoi"],
    classifiers=[],
)
