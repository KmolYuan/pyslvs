# -*- coding: utf-8 -*-

"""Compile the Cython libraries of Pyslvs."""

from os import listdir
from distutils.core import setup, Extension
from platform import system
from Cython.Distutils import build_ext
import numpy

np_include = numpy.get_include()


sources = []
for place in ["src/", "Adesign/src/"]:
    for source in listdir(place):
        if source.split('.')[-1] == 'pyx':
            sources.append(place + source)

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
        'src/' + 'bfgs.pyx',
        'src/bfgs_solver/' + 'geometric_constraint.cpp',
        'src/bfgs_solver/' + 'derivatives.cpp',
        'src/bfgs_solver/' + 'solve.cpp',
        'src/bfgs_solver/' + 'calc.cpp',
    ],
    language="c++",
    include_dirs=['src/bfgs_solver/', np_include],
    define_macros=macros,
    extra_compile_args=compile_args,
)]

for source in sources:
    if source == "src/bfgs.pyx":
        continue
    ext_modules.append(Extension(
        # Base name
        source.split('/')[-1].split('.')[0],
        # path + file name
        sources=[source],
        language="c++",
        include_dirs=[np_include],
        define_macros=macros,
        extra_compile_args=compile_args,
    ))

setup(ext_modules=ext_modules, cmdclass={'build_ext': build_ext})
