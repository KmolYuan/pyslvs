# -*- coding: utf-8 -*-

"""Compile the Cython libraries of Pyslvs."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2019"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

import os
import re
import codecs
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from platform import system

here = os.path.abspath(os.path.dirname(__file__))
src_path = 'pyslvs/'
bfgs_path = src_path + 'bfgs_solver/'
adesign_path = src_path + 'Adesign/'


def read(*parts):
    with codecs.open(os.path.join(here, *parts), 'r') as f:
        return f.read()


def find_version(*file_paths):
    m = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", read(*file_paths), re.M)
    if m:
        return m.group(1)
    raise RuntimeError("Unable to find version string.")


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
    src_path.replace('/', '.') + 'bfgs',
    [
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

git_modules = [src_path, adesign_path]
for place in git_modules:
    for source in os.listdir(place):
        if not source.endswith('.pyx'):
            continue
        if place == src_path and source == 'bfgs.pyx':
            continue
        ext_modules.append(Extension(
            place.replace('/', '.') + source.split('.')[0],  # Base name
            [place + source],
            language="c++",
            include_dirs=[],
            define_macros=macros,
            extra_compile_args=compile_args
        ))


class Build(build_ext):
    def finalize_options(self):
        super(Build, self).finalize_options()
        # Prevent numpy from thinking it is still in its setup process
        if isinstance(__builtins__, dict):
            __builtins__['__NUMPY_SETUP__'] = False
        else:
            setattr(__builtins__, '__NUMPY_SETUP__', False)
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
    package_data={'': ["*.pyi"]},
    ext_modules=ext_modules,
    cmdclass={'build_ext': Build},
    zip_safe=False,
    python_requires=">=3.6",
    setup_requires=[
        'setuptools',
        'wheel',
        'numpy',
        'cython',
    ],
    install_requires=[
        'numpy',
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
