# -*- coding: utf-8 -*-

from os import listdir
from os.path import sep, join as pth_join
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from platform import system

src_path = 'pyslvs'
graph_path = pth_join(src_path, 'graph')
bfgs_path = pth_join(src_path, 'bfgs_solver')
tinycadlib_path = pth_join(src_path, 'tinycadlib')
metaheuristics_path = pth_join(src_path, 'metaheuristics')
macros = [('_USE_MATH_DEFINES', None)]
compile_args_msvc = ['/O2', '/std:c++17']  # MSVC disabled OpenMP
compile_args = ['-Wno-cpp', '-std=c++17', '-fopenmp']
link_args = ['-fopenmp']
link_args_msvc = []
link_args_static = [
    '-static-libgcc',
    '-static-libstdc++',
    '-Wl,-Bstatic,--whole-archive',
    '-lwinpthread',
    '-lgomp',
    '-Wl,--no-whole-archive',
]
if system() == 'Windows':
    # Disable format warning
    compile_args.append('-Wno-format')
# Disable NumPy warning
macros.append(('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'))

ext_modules = [
    Extension(src_path.replace(sep, '.') + '.bfgs', [
        pth_join(src_path, 'bfgs.pyx'),
        pth_join(bfgs_path, 'constraints.cpp'),
        pth_join(bfgs_path, 'solve.cpp'),
        pth_join(bfgs_path, 'calc.cpp'),
    ], language="c++", include_dirs=[bfgs_path]),
    Extension(src_path.replace(sep, '.') + '.tinycadlib', [
        pth_join(src_path, 'tinycadlib.pyx'),
        pth_join(tinycadlib_path, 'solver.cpp'),
    ], language="c++", include_dirs=[tinycadlib_path]),
]
for place in [src_path, graph_path, metaheuristics_path]:
    for source in listdir(place):
        if not source.endswith('.pyx'):
            continue
        if place == src_path and source in {'bfgs.pyx', 'tinycadlib.pyx'}:
            continue
        ext_modules.append(Extension(
            place.replace(sep, '.') + '.' + source.split('.')[0],  # Base name
            [pth_join(place, source)],
            language="c++"
        ))
for ext in ext_modules:
    ext.cython_directives = {'binding': True, 'cdivision': True}


class Build(build_ext):

    def build_extensions(self):
        compiler = self.compiler.compiler_type
        if compiler in {'mingw32', 'unix'}:
            for e in self.extensions:
                e.define_macros = macros
                e.extra_compile_args = compile_args
                if compiler == 'mingw32':
                    e.extra_link_args = link_args_static
                else:
                    e.extra_link_args = link_args
        elif compiler == 'msvc':
            for e in self.extensions:
                e.define_macros = macros[:1]
                e.extra_compile_args = compile_args_msvc
                e.extra_link_args = link_args_msvc
        super(Build, self).build_extensions()


setup(ext_modules=ext_modules, cmdclass={'build_ext': Build})
