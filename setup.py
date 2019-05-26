# -*- coding: utf-8 -*-

"""Compile the Cython libraries of Pyslvs."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2019"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from os import listdir
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext as _build_ext
from platform import system

src_path = 'src/'
bfgs_path = src_path + 'bfgs_solver/'

"""
with open("__init__.py", "r") as f:
    for line in f:
        if line.startswith('__version__'):
            __version__ = line.split()[-1]
            break
"""
with open("README.md", "r") as f:
    long_description = f.read()

macros = [
    ('_hypot', 'hypot'),
]

compile_args = [
    '-O3',
    '-Wno-cpp',
]

if system() == 'Windows':
    # Avoid compile error with CYTHON_USE_PYLONG_INTERNALS.
    # https://github.com/cython/cython/issues/2670#issuecomment-432212671
    macros.append(('MS_WIN64', None))
    # Disable format warning.
    compile_args.append('-Wno-format')
    # Disable NumPy warning.
    macros.append(('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'))

ext_modules = [Extension(
    "bfgs",
    sources=[
        src_path + 'bfgs.pyx',
        bfgs_path + 'geometric_constraint.cpp',
        bfgs_path + 'derivatives.cpp',
        bfgs_path + 'solve.cpp',
        bfgs_path + 'calc.cpp',
    ],
    language="c++",
    include_dirs=[bfgs_path],
    define_macros=macros,
    extra_compile_args=compile_args
)]

for place in [src_path, "Adesign/src/"]:
    for source in listdir(place):
        source = place + source
        if not source.endswith('.pyx'):
            continue
        if source == "src/bfgs.pyx":
            continue
        ext_modules.append(Extension(
            source.split('/')[-1].split('.')[0],  # Base name
            sources=[source],
            language="c++",
            include_dirs=[],
            define_macros=macros,
            extra_compile_args=compile_args
        ))


class Build(_build_ext):
    def finalize_options(self):
        super(Build, self).finalize_options()
        # Prevent numpy from thinking it is still in its setup process
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())


setup(
    name='pyslvs',
    # version="18.5.0.dev0",
    author=__author__,
    author_email=__email__,
    description="Python library of Solvespace",
    long_description=long_description,
    url="https://github.com/solvespace/solvespace",
    packages=find_packages(),
    ext_modules=ext_modules,
    cmdclass={'build_ext': Build},
    setup_requires=[
        'setuptools>=18.0',
        'wheel',
        'numpy',
        'cython',
        'lark-parser',
        'pygments',
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Cython",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Operating System :: OS Independent",
    ]
)
