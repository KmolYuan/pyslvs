# -*- coding: utf-8 -*-

"""Compile the Cython libraries of Pyslvs."""

from distutils.core import setup, Extension
import os

from Cython.Distutils import build_ext
import numpy
np_include = numpy.get_include()


sources = []
for source in os.listdir("./src"):
    if source.split('.')[-1] == 'pyx':
        sources.append(source)

extra_compile_args = [
    # Compiler optimize.
    '-O3',
    # Disable NumPy warning only on Linux.
    '-Wno-cpp',
    # Avoid C++ math library.
    '-D_hypot=hypot',
    # Avoid compile error with CYTHON_USE_PYLONG_INTERNALS.
    # https://github.com/cython/cython/issues/2670#issuecomment-432212671
    '-DMS_WIN64',
]

ext_modules = [Extension(
    "bfgs",
    sources=[
        'src/' + 'bfgs.pyx',
        'src/bfgs_solver/' + 'geometric_constraint.cpp',
        'src/bfgs_solver/' + 'derivatives.cpp',
        'src/bfgs_solver/' + 'solve.cpp',
        'src/bfgs_solver/' + 'calc.cpp',
    ],
    language="c++",
    include_dirs=['src/bfgs_solver/', np_include],
    extra_compile_args=extra_compile_args,
)]

for source in sources:
    if source == "bfgs.pyx":
        continue
    ext_modules.append(Extension(
        # Base name
        source.split('.')[0],
        # path + file name
        sources=['src/' + source],
        language="c++",
        include_dirs=[np_include],
        extra_compile_args=extra_compile_args,
    ))

setup(ext_modules=ext_modules, cmdclass={'build_ext': build_ext})
