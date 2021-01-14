# Version 21.01

+ Add examples:
    + Stephenson III (slider)
    + Watt II (slider)
+ Remove "Coordinate" name alignment to "Coord".
+ Follows PEP517.

# Version 20.12

+ Add uniform four-bar linkage function.

# Version 20.11

+ Fix slider input function.

# Version 20.10

+ Unify "true" and "false" to lowercase in docstring.
+ Remove "data_collecting" function.
+ Remove "expr_parser" function.
+ Add "VPoint.is_slider" method.
+ Improve performance of meta-heuristic algorithms.
    + Release almost GIL, only left objective function.

# Version 20.09

+ This version has been skipped.

# Version 20.08

+ Change the definition of Path Signature: The mapping of stroke and curvature.

# Version 20.07

+ Rename "Coordinate" to "Coord"(keep old name).
+ Add "sx" & "sy" properties for VPoint slot coordinate.
+ Add "link_pos" and "points_pos" methods for VPoint and VLink respectively.
+ Add "to_coord" method for VPoint objects.
+ Reduce unused code of Sketch Solve and redesign its API.
+ Add PPP function to solve parallel linkage.
+ Rename the module "topo_config".

# Version 20.06

+ Support MSVC versions `1910`~`1916` and `1920`~`1926`.
+ Derivation will treat a path as a full stroke if difference is not the maximum.
+ Add function "PALP".

# Version 20.05

+ Remove "return none" annotation from stubs.
+ A shell version of "set_pycompiler".
+ Support incomplete path works with derivative function.

# Version 20.04

+ Fix the bug of topological algorithm.
+ Add `curvature`, `derivative` and `path_signature` functions.

# Version 20.03

+ Allow case insensitive color name.
+ Merge `altas` and `number` module as `structural`.
+ Rename `planar_check` module to `graph_planar`.
+ Add ordered simples option of planar linkage synthesis.
+ Change several variable names of backend.
+ Move graph related modules into `pyslvs.graph` package.
+ New graph searching method for polygon link.
    + Change triangular symbol "A" to "I" (input).
    + Add symbol "A" (angle of polygon link).

# Version 20.02

+ Change example list and collection list as a function.
+ Add an iterator for grabbing examples and collections names.
+ Function name typo corrections.
+ Add "get_include" function for setup.
+ "Shape only" synthesis mode.
+ DWT synthesis mode. (alpha)

# Version 20.01

+ Remove commas of link expression.
+ Turn on Cython "binding" option.

# Version 19.12

+ None.

# Version 19.11

+ All setting keys turn into lowercase.
+ Rename submodule "Adesign" to "metaheuristics".
+ Replace "time.perf_counter()" to "time.process_time()".
+ Correction of the random seed in Cython module.
+ Implement TLBO.
+ Add a simple test objective function.
+ Remove the entry point of unit test.

# Version 19.10

+ Manage Sketch Solve by CMake.
+ Use C++17 as default.
+ Make "Graph.degrees" public.
+ Add "Graph.adjacency_matrix" method.
+ Isomorphism:
    + More methods in "Graph" class.
    + Implement degree code for Graph class.
+ Remove the hints of Python objects in Cython sources.
+ Change the term "nodes" to "vertices".
+ Replace 'ground' string as "VLink.FRAME" class attribute.
+ Add placeholder "HOLDER" class variable for VPoint and VLink.
+ Add EFD algorithm.
+ Support MSVC compiler.
+ Support PEP 561.

# Version 19.09

+ Add test suit option of `setup.py`.
+ Replace "time.time" with "time.perf_counter".
+ Add Windows patch script `set_pycompiler`.

# Version 19.08

+ Add methods for Sketch Solve solver:
    + same_points
    + show_inputs
    + show_data
    + set_data
    + Add test cases for new methods.
+ Automatic fill up the link length options of data keys.
+ Add "Algorithm" class to share the same functions of three algorithms.
+ Add "slow down" option for algorithm stop limitations.
+ Member types optimization of algorithm classes.
+ Add "Horse leg" example.

# Version 19.07

+ Add "set attribute" methods for expression classes.
+ Add "duplicate" method for graph class.
+ Some corrections of stub files.
+ Adjust Cython "cimport" to related import instead of direct import.
+ Change Sketch Solve solver function to a class.

# Version 19.06

+ None.

# Version 19.05

+ Add "Six bar linkage mechanism" configuration.
+ Speed up performance of Cartesian product function.
+ Nested do loop method for contracted graph enumeration.
    + Speed up performance of searching contracted graphs.
+ Function renaming of "link_synthesis", "contracted_link_synthesis" and "conventional_graph".
+ Allow empty links in mechanism expression grammar.
+ Changed into standard PyPA module.

# Version 19.04

+ More methods of graph class.
+ Support multi-graph for "Planar check" function.
+ Add contracted graph enumeration algorithm:
    + Speed up performance for searching between link assortment.

# Version 19.03

+ String `__version__` provided.
+ Triangulation data type changed for performance improvement.
+ Change "print" function as "logger.debug".
+ Embed the "Coordinate" object in return value to improve performance.
+ Use "is" instead of "==" in type comparison.
+ Enable C division as default for all Cython sources.

# Version 19.02

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

# Version 19.01

+ This version has been skipped.

# Version 18.12

+ Outer loop layout for atlas.
+ Add "no global interpreter lock" decorator for pure C functions.
+ Added "edges_view" and "graph2vpoints" functions.
+ Update SHELL variable of Makefile on Windows platform.
+ Optimization of joint type of "VPoint" class with "VJoint" enumeration class.
+ Using C style type definitions of all class constructors.
+ Adjustment of input format for BFGS solver function.

# Version 18.11

+ Split out contracted link function to number synthesis.
+ Add planar graph checking function from NetworkX.
+ Performance optimization for instance creation.
+ Reverse analysis functions of graph class.

# Version 18.10

+ Improvements of "metaheuristics" kernel.
+ Three types of degenerate filter for 'topo' function.
+ High performance improvements of atlas algorithm.

# Version 18.09.1

+ Add PYI type hints for libraries.

# Version 18.09

+ Simulation corrections.
+ Add joint types of "VPoint" class.
+ Using "RuntimeError" instead of normal exception.
+ Change module name of "atlas".
+ Add support for slider variable input.
+ Add more information in atlas algorithm.
+ Change atlas algorithm "cancel" behavior to "skip".

# Version 18.08

+ Add "same link" operator for "VPoint" class.
+ Fix bug of BFGS transform in "tinycadlib".
+ Remove test function of Sketch Solve.
+ Support for multi slider joint.

# Version 18.07

+ Limitation format (upper and lower) of "Planar" class has been changed.
+ PMKS expression classes has been independent to new library.
+ Copy and comparison method of "VPoint" class.
+ New detect mechanism in VPoint distance method.
+ Remove "expr_path" function.
+ "Sketch Solve" kernel support PMKS expression.

# Version 18.06

+ Submodule added.
+ Add a C++ tiny CAD solver "Sketch Solve" with BFGS algorithm:
    + Cython wrapper.
    + Use to make position analysis.
+ Base class for algorithm verification function object.
