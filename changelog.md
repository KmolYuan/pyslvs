Change Log
===

Version 18.11
===

+ Split out contracted link function to number synthesis.

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
