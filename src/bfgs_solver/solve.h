#pragma once

/*
 *  Created on: May 4, 2009
 *  Author: Jonathan
 *  Contributor: KmolYuan
 */

///////////////////////////////////////
/// BFGS Solver parameters
///////////////////////////////////////

#define PertMag            1e-6
#define PertMin            1e-10
#define XConvergenceRough  1e-8
#define XConvergenceFine   1e-10
#define SmallF             1e-20
#define ValidSolutionFine  1e-12
#define ValidSoltuionRough 1e-4
#define Rough              0
#define Fine               1
//Note that the total number of iterations allowed is MaxIterations *xLength
#define MaxIterations      50

///////////////////////////////////////
/// Solve exit codes
///////////////////////////////////////

#define Succsess   0
#define NoSolution 1

///////////////////////////////////////
/// Expression data
///////////////////////////////////////

struct Point { double *x, *y; };

struct Line { Point *p1, *p2; };

struct Arc {
    Point *center;
    double *startAngle, *endAngle, *rad;
};

struct Circle {
    Point *center;
    double *rad;
};

struct GeoConstraint {
    int type;
    Point *point1, *point2;
    Line *line1, *line2, *SymLine;
    Circle *circle1, *circle2;
    Arc *arc1, *arc2;
    double *parameter; //radius, length, angle etc...
};

//GeoConstraint Functions
GeoConstraint PointOnPointConstraint(Point *, Point *);
GeoConstraint P2PDistanceConstraint(Point *, Point *, double *);
GeoConstraint P2PDistanceVertConstraint(Point *, Point *, double *);
GeoConstraint P2PDistanceHorzConstraint(Point *, Point *, double *);
GeoConstraint PointOnLineConstraint(Point *, Line *);
GeoConstraint P2LDistanceConstraint(Point *, Line *, double *);
GeoConstraint P2LDistanceVertConstraint(Point *, Line *, double *);
GeoConstraint P2LDistanceHorzConstraint(Point *, Line *, double *);
GeoConstraint VerticalConstraint(Line *);
GeoConstraint HorizontalConstraint(Line *);
GeoConstraint TangentToCircleConstraint(Line *, Circle *);
GeoConstraint TangentToArcConstraint(Line *, Arc *);
GeoConstraint ArcRulesConstraint(Arc *);
GeoConstraint LineLengthConstraint(Line *, double *);
GeoConstraint EqualLegnthConstraint(Line *, Line *);
GeoConstraint ArcRadiusConstraint(Arc *, double *);
GeoConstraint EqualRadiusArcsConstraint(Arc *, Arc *);
GeoConstraint EqualRadiusCirclesConstraint(Circle *, Circle *);
GeoConstraint EqualRadiusCircArcConstraint(Circle *, Arc *);
GeoConstraint ConcentricArcsConstraint(Arc *, Arc *);
GeoConstraint ConcentricCirclesConstraint(Circle *, Circle *);
GeoConstraint ConcentricCircArcConstraint(Circle *, Arc *);
GeoConstraint CircleRadiusConstraint(Circle *, double *);
GeoConstraint InternalAngleConstraint(Line *, Line *, double *);
GeoConstraint ExternalAngleConstraint(Line *, Line *, double *);
GeoConstraint PerpendicularConstraint(Line *, Line *);
GeoConstraint ParallelConstraint(Line *, Line *);
GeoConstraint ColinearConstraint(Line *, Line *);
GeoConstraint PointOnCircleConstraint(Point *, Circle *);
GeoConstraint PointOnArcConstraint(Point *, Arc *);
GeoConstraint PointOnLineMidpointConstraint(Point *, Line *);
GeoConstraint PointOnArcMidpointConstraint(Point *, Arc *);
GeoConstraint PointOnCircleQuadConstraint(Point *, Circle *, double *);
GeoConstraint SymmetricPointsConstraint(Point *, Point *, Line *);
GeoConstraint SymmetricLinesConstraint(Line *, Line *, Line *);
GeoConstraint SymmetricCirclesConstraint(Circle *, Circle *, Line *);
GeoConstraint SymmetricArcsConstraint(Arc *, Arc *, Line *);
GeoConstraint LineInternalAngleConstraint(Line *, double *);
GeoConstraint LineExternalAngleConstraint(Line *, double *);

//Public Functions
int solve(double **, int, GeoConstraint *, int, int);
void derivatives(double **, double *, int, GeoConstraint *, int);
