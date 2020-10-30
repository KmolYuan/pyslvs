# -*- coding: utf-8 -*-

"""Compile the Cython libraries of Pyslvs."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2020"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from re import MULTILINE, search
from os import listdir
from os.path import sep, join as pth_join
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from platform import system


def read(path: str):
    with open(path, 'r', encoding='utf-8') as f:
        return f.read()


def find_version(path: str):
    m = search(r"^__version__ = ['\"]([^'\"]*)['\"]", read(path), MULTILINE)
    if m:
        return m.group(1)
    raise RuntimeError("Unable to find version string.")


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
compiler_directives = {'binding': True, 'cdivision': True}

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
paths = [src_path, graph_path, metaheuristics_path]
for place in paths:
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
    ext.cython_directives = compiler_directives


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


setup(
    name='pyslvs',
    version=find_version(pth_join('pyslvs', '__init__.py')),
    author=__author__,
    author_email=__email__,
    license=__license__,
    description="Pyslvs library",
    long_description=read("README.md"),
    long_description_content_type='text/markdown',
    url="https://github.com/KmolYuan/pyslvs",
    packages=find_packages(exclude=('test',)),
    package_data={'': ['*.pyi', '*.pxd', '*.pyx'], 'pyslvs': ['py.typed']},
    ext_modules=ext_modules,
    cmdclass={'build_ext': Build},
    zip_safe=False,
    python_requires=">=3.7",
    install_requires=read('requirements.txt').splitlines(),
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Cython",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
        "Operating System :: OS Independent",
        "Typing :: Typed",
    ]
)
