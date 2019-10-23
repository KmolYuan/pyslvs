#pragma once

/*
 *  Created on: May 4, 2009
 *  Author: Jonathan
 *  Contributor: KmolYuan
 */

#include <cstddef>

///////////////////////////////////////
/// BFGS Solver parameters
///////////////////////////////////////

#define PertMag 1e-6
#define PertMin 1e-10
#define XConvergenceRough 1e-8
#define XConvergenceFine 1e-10
#define SmallF 1e-20
#define ValidSolutionFine 1e-12
#define ValidSoltuionRough 1e-4
#define Rough 0
#define Fine 1
// Note that the total number of iterations allowed is MaxIterations *xLength
#define MaxIterations 50

///////////////////////////////////////
/// Solve exit codes
///////////////////////////////////////

#define Success 1
#define NoSolution 0

///////////////////////////////////////
/// Position Expression data
///////////////////////////////////////

struct Point {
    double *x, *y;
};

struct Line {
    Point *p1, *p2;
};

struct Arc {
    Point *center;
    double *startAngle, *endAngle, *rad;
};

struct Circle {
    Point *center;
    double *rad;
};

enum {
    // Geometric Constraint types
    PointOnPoint,
    PointToLine,
    PointOnLine,
    Horizontal,
    Vertical,
    InternalAngle,
    RadiusValue,
    TangentToArc,
    TangentToCircle,
    ArcRules,
    P2PDistance,
    P2PDistanceVert,
    P2PDistanceHorz,
    P2LDistance,
    P2LDistanceVert,
    P2LDistanceHorz,
    LineLength,
    EqualLegnth,
    ArcRadius,
    EqualRadiusArcs,
    EqualRadiusCircles,
    EqualRadiusCircArc,
    ConcentricArcs,
    ConcentricCircles,
    ConcentricCircArc,
    CircleRadius,
    ExternalAngle,
    Parallel,
    Perpendicular,
    Colinear,
    PointOnCircle,
    PointOnArc,
    PointOnLineMidpoint,
    PointOnArcMidpoint,
    PointOnCircleQuad,
    SymmetricPoints,
    SymmetricLines,
    SymmetricCircles,
    SymmetricArcs,
    LineInternalAngle,
    LineExternalAngle,
};

struct Constraint {
    unsigned long type;
    Point *point1, *point2;
    Line *line1, *line2, *SymLine;
    Circle *circle1, *circle2;
    Arc *arc1, *arc2;
    // RADIUS, length, angle etc...
    double *parameter;
};

///////////////////////////////////////
/// Constraint Functions (for safe access of mebers)
///////////////////////////////////////

// Geometric Constrants
Constraint PointOnPointConstraint(Point *, Point *);
Constraint P2PDistanceConstraint(Point *, Point *, double *);
Constraint P2PDistanceVertConstraint(Point *, Point *, double *);
Constraint P2PDistanceHorzConstraint(Point *, Point *, double *);
Constraint PointOnLineConstraint(Point *, Line *);
Constraint P2LDistanceConstraint(Point *, Line *, double *);
Constraint P2LDistanceVertConstraint(Point *, Line *, double *);
Constraint P2LDistanceHorzConstraint(Point *, Line *, double *);
Constraint VerticalConstraint(Line *);
Constraint HorizontalConstraint(Line *);
Constraint TangentToCircleConstraint(Line *, Circle *);
Constraint TangentToArcConstraint(Line *, Arc *);
Constraint ArcRulesConstraint(Arc *);
Constraint LineLengthConstraint(Line *, double *);
Constraint EqualLegnthConstraint(Line *, Line *);
Constraint ArcRadiusConstraint(Arc *, double *);
Constraint EqualRadiusArcsConstraint(Arc *, Arc *);
Constraint EqualRadiusCirclesConstraint(Circle *, Circle *);
Constraint EqualRadiusCircArcConstraint(Circle *, Arc *);
Constraint ConcentricArcsConstraint(Arc *, Arc *);
Constraint ConcentricCirclesConstraint(Circle *, Circle *);
Constraint ConcentricCircArcConstraint(Circle *, Arc *);
Constraint CircleRadiusConstraint(Circle *, double *);
Constraint InternalAngleConstraint(Line *, Line *, double *);
Constraint ExternalAngleConstraint(Line *, Line *, double *);
Constraint PerpendicularConstraint(Line *, Line *);
Constraint ParallelConstraint(Line *, Line *);
Constraint ColinearConstraint(Line *, Line *);
Constraint PointOnCircleConstraint(Point *, Circle *);
Constraint PointOnArcConstraint(Point *, Arc *);
Constraint PointOnLineMidpointConstraint(Point *, Line *);
Constraint PointOnArcMidpointConstraint(Point *, Arc *);
Constraint PointOnCircleQuadConstraint(Point *, Circle *, double *);
Constraint SymmetricPointsConstraint(Point *, Point *, Line *);
Constraint SymmetricLinesConstraint(Line *, Line *, Line *);
Constraint SymmetricCirclesConstraint(Circle *, Circle *, Line *);
Constraint SymmetricArcsConstraint(Arc *, Arc *, Line *);
Constraint LineInternalAngleConstraint(Line *, double *);
Constraint LineExternalAngleConstraint(Line *, double *);

///////////////////////////////////////
/// Public Functions
///////////////////////////////////////

int solve(double **, size_t, Constraint *, size_t, bool);
void derivatives(double **, double *, size_t, Constraint *, size_t);
