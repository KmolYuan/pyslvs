[![Version](https://img.shields.io/badge/version-18.11.0-yellow.svg)](https://github.com/KmolYuan/pyslvs/releases/latest)
[![Build status](https://ci.appveyor.com/api/projects/status/6l1bh1197ncahd0q?svg=true)](https://ci.appveyor.com/project/KmolYuan/pyslvs)
[![Build status](https://img.shields.io/travis/KmolYuan/pyslvs.svg?logo=travis)](https://travis-ci.org/KmolYuan/pyslvs)

# Pyslvs Libraries

A no-GUI module of mechanism synthesis system and a 2D geometric constraint solver.

## Build and Test

Enter directory and execute the Makefile. Then, run the unit test script after compiling.

```bash
make
python test_pyslvs.py
```

## Module parts

Pyslvs libraries is divided into two following sections:

+ Solver:

    Geometric solver and verification functions.

+ Synthesis:

    Mechanism synthesis system that including several random algorithm and enumeration algorithm, dependent with geometric solver.

Most of classes and functions can be work with a generic Python format (just like a list of coordinates or string-like expression), and you also can cut in at any step.

### Solver

In Pyslvs, expression is mainly to [PMKS](http://designengrlab.github.io/PMKS/), which is more 

+ **parser** module:

    Analysis expression from strings, turn into symbols object. There also including highlighting function with [Pygments](http://pygments.org/).

+ **pmks** library:

    Including PMKS expression object classes.

+ **bfgs** library:

    **preload**: *pmks*

    Python wrapper of [Sketchsolve](https://code.google.com/archive/p/sketchsolve/). A simple and fast constraint solver with BFGS algorithm.

+ **tinycadlib** library:

    **preload**: *pmks*, *bfgs*

    Particular solution takes more faster then constraint solving.

+ **triangulation** library:

    **preload**: *pmks*, *tinycadlib*

    Autometic configuration algorithm for particular solution function in "tinycadlib".

### Graph Synthesis

+ **number** library:

    Number synthesis function for searching solutions of the number of joints and links.

+ **graph** library:

    Graph expression came from NetworkX. Also contains graph verification functions. 

+ **planar_check** library:

    **preload**: *graph*

    Planar graph checking came from NetworkX.

+ **atlas** library:

    **preload**: *graph*, *planar_check*

    Graph combination algorithm.

#### Adesign (Dimensional Synthesis)

[Adesign](https://github.com/KmolYuan/Adesign) module: Cython algorithms libraries provide evolution designing.

+ **verify** library:

    Provide base fitness function class for algorithm.

+ **planarlinkage** library:

    **preload**: *pmks*, *bfgs*, *tinycadlib*, *verify*

    Dimensional synthesis verification function objects.

+ **rga** library:

    **preload**: *verify*

    Real-coded genetic algorithm for dimensional synthesis.

+ **firefly** library:

    **preload**: *verify*

    Firefly algorithm for dimensional synthesis.

+ **de** library:

    **preload**: *verify*

    Differential Evolution for dimensional synthesis.
