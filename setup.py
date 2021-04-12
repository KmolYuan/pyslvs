# -*- coding: utf-8 -*-

from os import walk
from os.path import sep, join as pth_join
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from platform import system

root_path = 'pyslvs'
bfgs_path = pth_join(root_path, 'bfgs_solver')
tinycadlib_path = pth_join(root_path, 'tinycadlib')
macros = [('_USE_MATH_DEFINES', None)]
compile_args_msvc = ['/O2', '/std:c++17']
compile_args = ['-Wno-cpp', '-std=c++17']
link_args = []
link_args_msvc = []
link_args_static = [
    '-static-libgcc',
    '-static-libstdc++',
    '-Wl,-Bstatic,--whole-archive',
    '-lwinpthread',
    '-Wl,--no-whole-archive',
]
if system() == 'Windows':
    # Disable format warning
    compile_args.append('-Wno-format')
# Disable NumPy warning
macros.append(('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'))
# Special modules mixed with C++ (language must be unified)
ext_modules = [
    Extension(root_path.replace(sep, '.') + '.bfgs', [
        pth_join(root_path, 'bfgs.pyx'),
        pth_join(bfgs_path, 'constraints.cpp'),
        pth_join(bfgs_path, 'solve.cpp'),
        pth_join(bfgs_path, 'calc.cpp'),
    ], language="c++", include_dirs=[bfgs_path]),
    Extension(root_path.replace(sep, '.') + '.tinycadlib', [
        pth_join(root_path, 'tinycadlib.pyx'),
        pth_join(tinycadlib_path, 'solver.cpp'),
    ], language="c++", include_dirs=[tinycadlib_path]),
]
# To find other sources recursively
for root, _, files in walk(root_path):
    for source in files:
        if not source.endswith('.pyx'):
            continue
        if root == root_path and source in {'bfgs.pyx', 'tinycadlib.pyx'}:
            continue
        f_name = pth_join(root, source)
        ext_modules.append(Extension(
            f_name.replace(sep, '.').rsplit('.', maxsplit=1)[0],  # Name
            [f_name],
            language="c++", include_dirs=list({root_path, root})
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
