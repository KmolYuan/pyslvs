# -*- coding: utf-8 -*-

"""Compile the Cython libraries of Pyslvs."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2019"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from os import listdir
from os.path import sep, join as pth_join
import re
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from platform import system


def read(*parts):
    with open(pth_join(*parts), 'r', encoding='utf-8') as f:
        return f.read()


def find_version(*file_paths):
    m = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", read(*file_paths), re.M)
    if m:
        return m.group(1)
    raise RuntimeError("Unable to find version string.")


src_path = 'pyslvs'
bfgs_path = pth_join(src_path, 'bfgs_solver')
adesign_path = pth_join(src_path, 'Adesign')
macros = [('_hypot', 'hypot'), ('M_PI', 'PI')]
compile_args = ['-O3', '-Wno-cpp']
if system() == 'Windows':
    # Avoid compile error with CYTHON_USE_PYLONG_INTERNALS
    # https://github.com/cython/cython/issues/2670#issuecomment-432212671
    macros.append(('MS_WIN64', None))
    # Disable format warning
    compile_args.append('-Wno-format')
# Disable NumPy warning
macros.append(('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'))

ext_modules = [Extension(
    src_path.replace(sep, '.') + '.bfgs',
    [
        pth_join(src_path, 'bfgs.pyx'),
        pth_join(bfgs_path, 'geometric_constraint.cpp'),
        pth_join(bfgs_path, 'derivatives.cpp'),
        pth_join(bfgs_path, 'solve.cpp'),
        pth_join(bfgs_path, 'calc.cpp'),
    ],
    language="c++",
    include_dirs=[bfgs_path]
)]
for place in [src_path, adesign_path]:
    for source in listdir(place):
        if not source.endswith('.pyx'):
            continue
        if place == src_path and source == 'bfgs.pyx':
            continue
        ext_modules.append(Extension(
            place.replace(sep, '.') + '.' + source.split('.')[0],  # Base name
            [pth_join(place, source)],
            language="c++"
        ))


class Build(build_ext):
    def build_extensions(self):
        compiler = self.compiler.compiler_type
        if compiler in {'mingw32', 'unix'}:
            for e in self.extensions:
                e.define_macros = macros
                e.extra_compile_args = compile_args
        elif compiler == 'msvc':
            for e in self.extensions:
                e.define_macros = [('_USE_MATH_DEFINES', None)]
        super(Build, self).build_extensions()

    def finalize_options(self):
        super(Build, self).finalize_options()
        import numpy
        self.include_dirs.append(numpy.get_include())


setup(
    name='pyslvs',
    version=find_version('pyslvs', '__init__.py'),
    author=__author__,
    author_email=__email__,
    license=__license__,
    description="Pyslvs library",
    long_description=read("README.md"),
    long_description_content_type='text/markdown',
    url="https://github.com/KmolYuan/pyslvs",
    packages=find_packages(exclude=('tests',)),
    package_data={'': ["*.pyi", "*.pxd"], 'pyslvs': ['py.typed']},
    ext_modules=ext_modules,
    cmdclass={'build_ext': Build},
    zip_safe=False,
    python_requires=">=3.7",
    install_requires=read('requirements.txt').splitlines(),
    test_suite="tests",
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Cython",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
        "Operating System :: OS Independent",
    ]
)
