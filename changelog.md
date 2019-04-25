Change Log
===

Version 19.04
===

+ More methods of graph class.
+ Support multi-graph for "Planar check" function.
+ Add contracted graph enumeration algorithm:
    + Speed up performance for searching between link assortment.

Version 19.03
===

+ String `__version__` provided.
+ Triangulation data type changed for performance improvement.
+ Change "print" function as "logger.debug".
+ Embed the "Coordinate" object in return value to improve performance.
+ Use "is" instead of "==" in type comparison.
+ Enable C division as default for all Cython sources.

Version 19.02
===

+ Change the function name of triangle formulas.
+ Remove the external color format of the VPoint color attribute.
+ Remove the "four bar constraint" option of planar linkage class.
+ Grammar enhancement.
+ Module "pmks" renamed to "expression".
+ Module "planarlinkage" renamed to "planar_linkage".
+ Adjust limitation option of algorithm to "max_gen", "min_fit" and "max_time".
+ Add "cdef" static method for joint creation.
+ Change "links" argument of VPoint constructor into string iterable object.
+ Using normal exception instead of "RuntimeError".
+ Supported specified link length option.
+ Cancel the copy step of submodule pxd file.
+ Enhancement of mechanism grammar.
+ Add "-j" option to speed up compile time.
+ Many performance improvements.
+ Remove the definition of inner scope.
+ Fix the memory leak of BFGS solver.
+ Reformat BFGS cython wrapper for more readability.

Version 19.01
===

+ This version has been skipped.

Version 18.12
===

+ Outer loop layout for atlas.
+ Add "no global interpreter lock" decorator for pure C functions.
+ Added "edges_view" and "graph2vpoints" functions.
+ Update SHELL variable of Makefile on Windows platform.
+ Optimization of joint type of "VPoint" class with "VJoint" enumeration class.
+ Using C style type definitions of all class constructors.
+ Adjustment of input format for BFGS solver function.

Version 18.11
===

+ Split out contracted link function to number synthesis.
+ Add planar graph checking function from NetworkX.
+ Performance optimization for instance creation.
+ Reverse analysis functions of graph class.

Version 18.10
===

+ Improvements of "Adesign" kernel.
+ Three types of degenerate filter for 'topo' function.
+ High performance improvements of atlas algorithm.

Version 18.09.1
===

+ Add PYI type hints for libraries.

Version 18.09
===

+ Simulation corrections.
+ Add joint types of "VPoint" class.
+ Using "RuntimeError" instead of normal exception.
+ Change module name of "atlas".
+ Add support for slider variable input.
+ Add more information in atlas algorithm.
+ Change atlas algorithm "cancel" behavior to "skip".

Version 18.08
===

+ Add "same link" operator for "VPoint" class.
+ Fix bug of BFGS transform in "tinycadlib".
+ Remove test function of Sketch Solve.
+ Support for multi slider joint.

Version 18.07
===

+ Limitation format (upper and lower) of "Planar" class has been changed.
+ PMKS expression classes has been independent to new library.
+ Copy and comparison method of "VPoint" class.
+ New detect mechanism in VPoint distance method.
+ Remove "expr_path" function.
+ "Sketch Solve" kernel support PMKS expression.

Version 18.06
===

+ Submodule added.
+ Add a C++ tiny CAD solver "Sketch Solve" with BFGS algorithm:
    + Cython wrapper.
    + Use to make position analysis.
+ Base class for algorithm verification function object.
