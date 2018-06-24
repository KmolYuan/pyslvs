#pragma once

/*
 *  Created on: May 4, 2009
 *  Author: Jonathan
 *  Contributor: KmolYuan
 */

#ifndef WIN32
    #define _hypot hypot
#endif

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

///////////////////////////////////////
/// Constraint types
///////////////////////////////////////

#define PointOnPoint      0
#define PointToLine       1
#define PointOnLine       2
#define Horizontal        3
#define Vertical          4
#define InternalAngle     5
#define RadiusValue       6
#define TangentToArc      7
#define TangentToCircle   8
#define ArcRules          9
#define P2PDistance      10
#define P2PDistanceVert  11
#define P2PDistanceHorz  12
#define P2LDistance      13
#define P2LDistanceVert  14
#define P2LDistanceHorz  15
#define LineLength       16
#define EqualLegnth      17
#define ArcRadius        18
#define EqualRadiusArcs      19
#define EqualRadiusCircles   20
#define EqualRadiusCircArc   21
#define ConcentricArcs       22
#define ConcentricCircles    23
#define ConcentricCircArc    24
#define CircleRadius         25
#define ExternalAngle        26
#define Parallel             27
#define Perpendicular        28
#define Colinear             29
#define PointOnCircle        30
#define PointOnArc           31
#define PointOnLineMidpoint  32
#define PointOnArcMidpoint   33
#define PointOnCircleQuad    34
#define SymmetricPoints      35
#define SymmetricLines       36
#define SymmetricCircles     37
#define SymmetricArcs        38

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

///////////////////////////////////////////////////////////////////////
/// constraint defines (these make writing constraint equations easier
///////////////////////////////////////////////////////////////////////
#define P1_x         *cons[i].point1->x
#define P1_y         *cons[i].point1->y
#define P2_x         *cons[i].point2->x
#define P2_y         *cons[i].point2->y
#define L1_P1_x      *cons[i].line1->p1->x
#define L1_P1_y      *cons[i].line1->p1->y
#define L1_P2_x      *cons[i].line1->p2->x
#define L1_P2_y      *cons[i].line1->p2->y
#define L2_P1_x      *cons[i].line2->p1->x
#define L2_P1_y      *cons[i].line2->p1->y
#define L2_P2_x      *cons[i].line2->p2->x
#define L2_P2_y      *cons[i].line2->p2->y
#define C1_Center_x  *cons[i].circle1->center->x
#define C1_Center_y  *cons[i].circle1->center->y
#define C1_rad       *cons[i].circle1->rad
#define C2_Center_x  *cons[i].circle2->center->x
#define C2_Center_y  *cons[i].circle2->center->y
#define C2_rad       *cons[i].circle2->rad
#define A1_startA    *cons[i].arc1->startAngle
#define A1_endA      *cons[i].arc1->endAngle
#define A1_radius    *cons[i].arc1->rad
#define A1_Center_x  *cons[i].arc1->center->x
#define A1_Center_y  *cons[i].arc1->center->y
#define A2_startA    *cons[i].arc2->startAngle
#define A2_endA      *cons[i].arc2->endAngle
#define A2_radius    *cons[i].arc2->rad
#define A2_Center_x  *cons[i].arc2->center->x
#define A2_Center_y  *cons[i].arc2->center->y
#define A1_Start_x   (A1_Center_x + A1_radius * cos(A1_startA))
#define A1_Start_y   (A1_Center_y + A1_radius * sin(A1_startA))
#define A1_End_x     (A1_Center_x + A1_radius * cos(A1_endA))
#define A1_End_y     (A1_Center_y + A1_radius * sin(A1_endA))
#define A2_Start_x   (A1_Center_x + A2_radius * cos(A2_startA))
#define A2_Start_y   (A1_Center_y + A2_radius * sin(A2_startA))
#define A2_End_x     (A1_Center_x + A2_radius * cos(A2_endA))
#define A2_End_y     (A1_Center_y + A2_radius * sin(A2_endA))
#define length       *cons[i].parameter
#define distance     *cons[i].parameter
#define radius       *cons[i].parameter
#define angleP       *cons[i].parameter
#define quadIndex    *cons[i].parameter
#define Sym_P1_x     *cons[i].SymLine->p1->x
#define Sym_P1_y     *cons[i].SymLine->p1->y
#define Sym_P2_x     *cons[i].SymLine->p2->x
#define Sym_P2_y     *cons[i].SymLine->p2->y


struct Point { double *x, *y; };

struct Line { Point *p1, *p2; };

struct Arc {
    double *startAngle, *endAngle, *rad; //rad
    Point *center;
};

struct Circle {
    Point *center;
    double *rad;
};

struct Constraint {
    int type;
    Point *point1, *point2;
    Line *line1, *line2, *SymLine;
    Circle *circle1, *circle2;
    Arc *arc1, *arc2;
    double *parameter; //radius, length, angle etc...
};

//Constraint Functions
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
Constraint LineLengthConstraint(Line *);
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

//Public Functions
int solve(double **, int, Constraint *, int, int);
void derivatives(double **, double *, int, Constraint *, int);
