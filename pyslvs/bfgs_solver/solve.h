#ifndef SOLVE_H
#define SOLVE_H

/*
 *  Created on: May 4, 2009
 *  Author: Jonathan
 *  Contributor: KmolYuan
 */

#ifdef __cplusplus
#include <cstddef>
extern "C" {
#else
#include <stdbool.h>
#include <stddef.h>
#endif

///////////////////////////////////////
/// Solver parameters
///////////////////////////////////////

#define ROUGH 0
// Solve exit codes
#define SUCCESS 1
#define NO_SOLUTION 0

///////////////////////////////////////
/// Position Expression data
///////////////////////////////////////

struct Point {
    double *x, *y;
};

struct Line {
    Point *p1, *p2;
};

struct Constraint {
    unsigned type;
    Point *point1, *point2;
    Line *line1, *line2;
    double *parameter;
};

///////////////////////////////////////
/// Constraint Functions (for safe access of members)
///////////////////////////////////////

// Geometric Constraints
Constraint point_on_point(Point *, Point *);
Constraint p2p_distance(Point *, Point *, double *);
#ifdef __cplusplus
[[maybe_unused]]
#endif
Constraint point_on_line(Point *, Line *);
Constraint internal_angle(Line *, Line *, double *);
#ifdef __cplusplus
[[maybe_unused]]
#endif
Constraint line_internal_angle(Line *, double *);

///////////////////////////////////////
/// Public Functions
///////////////////////////////////////

int solve(double **, size_t, Constraint *, size_t, bool);
void derivatives(double **, double *, size_t, Constraint *, size_t);

#ifdef __cplusplus
}
#endif

#endif  // SOLVE_H
