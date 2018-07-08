# -*- coding: utf-8 -*-

"""Compile the Cython libraries of Pyslvs."""

from distutils.core import setup, Extension
import os

from Cython.Distutils import build_ext
import numpy
numpy_include = [numpy.get_include()]


sources = [
    source for source in os.listdir("./src")
    if source.split('.')[-1] == 'pyx'
]

extra_compile_args = [
    #Compiler optimize.
    '-O3',
    #Disable NumPy warning only on Linux.
    '-Wno-cpp',
    #Windows format warning.
    '-Wno-format',
]

ext_modules = [Extension(
    "bfgs",
    sources = [
        'src/' + 'bfgs.pyx',
        'src/bfgs_solver/' + 'constraint_func.cpp',
        'src/bfgs_solver/' + 'derivatives.cpp',
        'src/bfgs_solver/' + 'solve.cpp',
    ],
    language = "c++",
    include_dirs = ['src/bfgs_solver/'] + numpy_include,
    
    extra_compile_args = extra_compile_args + ['-D_hypot=hypot'],
)]

for source in sources:
    if source == "bfgs.pyx":
        continue
    ext_modules.append(Extension(
        source.split('.')[0], #Base name
        sources = ['src/' + source], #path + file name
        
        include_dirs = numpy_include,
        extra_compile_args = extra_compile_args,
    ))


setup(
    ext_modules = ext_modules,
    cmdclass = {'build_ext': build_ext},
)
