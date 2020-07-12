cmake_minimum_required(VERSION 3.0)
project(sketch_solve)

set(CMAKE_CXX_STANDARD 17)

add_executable(sketch_solve
    solve.h
    calc.h
    derivatives.cpp
    constraints.cpp
    solve.cpp
    calc.cpp
    main.cpp)
