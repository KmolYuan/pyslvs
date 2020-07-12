#ifndef SOLVE_H
#define SOLVE_H

/*
 *  Created on: May 4, 2009
 *  Author: Jonathan
 *  Contributor: KmolYuan
 */

#include <cstddef>

///////////////////////////////////////
/// Solver parameters
///////////////////////////////////////

#define PertMag 1e-6
#define PertMin 1e-10
#define XConvergenceRough 1e-8
#define XConvergenceFine 1e-10
#define SmallF 1e-20
#define ValidSolutionFine 1e-12
#define ValidSoltuionRough 1e-4
#define Rough 0
// Note that the total number of iterations allowed is MaxIterations *xLength
#define MaxIterations 50

// Solve exit codes
#define Success 1
#define NoSolution 0

#ifdef __cplusplus
extern "C" {
#endif

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
    double *parameter;  // Radius, length, angle etc.
};

///////////////////////////////////////
/// Constraint Functions (for safe access of members)
///////////////////////////////////////

// Geometric Constraints
Constraint PointOnPointConstraint(Point *, Point *);
Constraint P2PDistanceConstraint(Point *, Point *, double *);
[[maybe_unused]] Constraint PointOnLineConstraint(Point *, Line *);
Constraint InternalAngleConstraint(Line *, Line *, double *);
[[maybe_unused]] Constraint LineInternalAngleConstraint(Line *, double *);

///////////////////////////////////////
/// Public Functions
///////////////////////////////////////

int solve(double **, size_t, Constraint *, size_t, bool);
void derivatives(double **, double *, size_t, Constraint *, size_t);

#ifdef __cplusplus
}
#endif

#endif  // SOLVE_H
