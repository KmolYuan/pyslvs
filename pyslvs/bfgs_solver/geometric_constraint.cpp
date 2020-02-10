/*
 * Create on: Jun 6, 2019
 * Author: KmolYuan
 * These functions make writing constraint easier
 */

#include "calc.h"

namespace {
Constraint P2PDistanceConstraint(unsigned type, Point *point1, Point *point2,
                                 double *value) {
    Constraint con = {};
    con.type = type;
    con.point1 = point1;
    con.point2 = point2;
    con.parameter = value;
    return con;
}

Constraint PointLineConstraint(unsigned type, Point *point1, Line *line1) {
    Constraint con = {};
    con.type = type;
    con.point1 = point1;
    con.line1 = line1;
    return con;
}

Constraint P2LDistanceConstraint(unsigned type, Point *point1, Line *line1,
                                 double *value) {
    Constraint con = {};
    con.type = type;
    con.point1 = point1;
    con.line1 = line1;
    con.parameter = value;
    return con;
}

Constraint VerticalHorizontalConstraint(unsigned type, Line *line1) {
    Constraint con = {};
    con.type = type;
    con.line1 = line1;
    return con;
}

Constraint LineConstraint(unsigned type, Line *line1, double *value) {
    Constraint con = {};
    con.type = type;
    con.line1 = line1;
    con.parameter = value;
    return con;
}

Constraint LinesConstraint(unsigned type, Line *line1, Line *line2) {
    Constraint con = {};
    con.type = type;
    con.line1 = line1;
    con.line2 = line2;
    return con;
}

Constraint ArcsConstraint(unsigned type, Arc *arc1, Arc *arc2) {
    Constraint con = {};
    con.type = type;
    con.arc1 = arc1;
    con.arc2 = arc2;
    return con;
}

Constraint CirclesConstraint(unsigned type, Circle *circle1, Circle *circle2) {
    Constraint con = {};
    con.type = type;
    con.circle1 = circle1;
    con.circle2 = circle2;
    return con;
}

Constraint CircleArcConstraint(unsigned type, Circle *circle1, Arc *arc1) {
    Constraint con = {};
    con.type = type;
    con.circle1 = circle1;
    con.arc1 = arc1;
    return con;
}

Constraint AngleConstraint(unsigned type, Line *line1, Line *line2,
                           double *value) {
    Constraint con = {};
    con.type = type;
    con.line1 = line1;
    con.line2 = line2;
    con.parameter = value;
    return con;
}

Constraint PointCircleConstraint(unsigned type, Point *point1,
                                 Circle *circle1) {
    Constraint con = {};
    con.type = type;
    con.point1 = point1;
    con.circle1 = circle1;
    return con;
}

Constraint PointArcConstraint(unsigned type, Point *point1, Arc *arc1) {
    Constraint con = {};
    con.type = type;
    con.point1 = point1;
    con.arc1 = arc1;
    return con;
}
}  // namespace

Constraint PointOnPointConstraint(Point *point1, Point *point2) {
    Constraint con = {};
    con.type = PointOnPoint;
    con.point1 = point1;
    con.point2 = point2;
    return con;
}

Constraint P2PDistanceConstraint(Point *point1, Point *point2, double *value) {
    return P2PDistanceConstraint(P2PDistance, point1, point2, value);
}

Constraint P2PDistanceVertConstraint(Point *point1, Point *point2,
                                     double *value) {
    return P2PDistanceConstraint(P2PDistanceVert, point1, point2, value);
}

Constraint P2PDistanceHConstraint(Point *point1, Point *point2,
                                     double *value) {
    return P2PDistanceConstraint(P2PDistanceHorz, point1, point2, value);
}

Constraint PointOnLineConstraint(Point *point1, Line *line1) {
    return PointLineConstraint(PointOnLine, point1, line1);
}

Constraint P2LDistanceConstraint(Point *point1, Line *line1, double *value) {
    return P2LDistanceConstraint(P2LDistance, point1, line1, value);
}

Constraint P2LDistanceVertConstraint(Point *point1, Line *line1,
                                     double *value) {
    return P2LDistanceConstraint(P2LDistanceVert, point1, line1, value);
}

Constraint P2LDistanceHConstraint(Point *point1, Line *line1,
                                     double *value) {
    return P2LDistanceConstraint(P2LDistanceHorz, point1, line1, value);
}

Constraint VerticalConstraint(Line *line1) {
    return VerticalHorizontalConstraint(Vertical, line1);
}

Constraint HorizontalConstraint(Line *line1) {
    return VerticalHorizontalConstraint(Horizontal, line1);
}

Constraint TangentToCircleConstraint(Line *line1, Circle *circle1) {
    Constraint con = {};
    con.type = TangentToCircle;
    con.line1 = line1;
    con.circle1 = circle1;
    return con;
}

Constraint TangentToArcConstraint(Line *line1, Arc *arc1) {
    Constraint con = {};
    con.type = TangentToCircle;
    con.line1 = line1;
    con.arc1 = arc1;
    return con;
}

Constraint ArcRulesConstraint(Arc *arc1) {
    Constraint con = {};
    con.type = ArcRules;
    con.arc1 = arc1;
    return con;
}

Constraint LineLengthConstraint(Line *line1, double *value) {
    return LineConstraint(LineLength, line1, value);
}

Constraint EqualLengthConstraint(Line *line1, Line *line2) {
    return LinesConstraint(EqualLegnth, line1, line2);
}

Constraint ArcRadiusConstraint(Arc *arc1, double *value) {
    Constraint con = {};
    con.type = ArcRadius;
    con.arc1 = arc1;
    con.parameter = value;
    return con;
}

Constraint EqualRadiusArcsConstraint(Arc *arc1, Arc *arc2) {
    return ArcsConstraint(EqualRadiusArcs, arc1, arc2);
}

Constraint EqualRadiusCirclesConstraint(Circle *circle1, Circle *circle2) {
    return CirclesConstraint(EqualRadiusCircles, circle1, circle2);
}

Constraint EqualRadiusCircleArcConstraint(Circle *circle1, Arc *arc1) {
    return CircleArcConstraint(EqualRadiusCircles, circle1, arc1);
}

Constraint ConcentricArcsConstraint(Arc *arc1, Arc *arc2) {
    return ArcsConstraint(ConcentricArcs, arc1, arc2);
}

Constraint ConcentricCirclesConstraint(Circle *circle1, Circle *circle2) {
    return CirclesConstraint(ConcentricCircles, circle1, circle2);
}

Constraint ConcentricCircleArcConstraint(Circle *circle1, Arc *arc1) {
    return CircleArcConstraint(ConcentricCircArc, circle1, arc1);
}

Constraint CircleRadiusConstraint(Circle *circle1, double *value) {
    Constraint con = {};
    con.type = CircleRadius;
    con.circle1 = circle1;
    con.parameter = value;
    return con;
}

Constraint InternalAngleConstraint(Line *line1, Line *line2, double *value) {
    return AngleConstraint(InternalAngle, line1, line2, value);
}

Constraint ExternalAngleConstraint(Line *line1, Line *line2, double *value) {
    return AngleConstraint(ExternalAngle, line1, line2, value);
}

Constraint PerpendicularConstraint(Line *line1, Line *line2) {
    return LinesConstraint(Perpendicular, line1, line2);
}

Constraint ParallelConstraint(Line *line1, Line *line2) {
    return LinesConstraint(Parallel, line1, line2);
}

Constraint CollinearConstraint(Line *line1, Line *line2) {
    return LinesConstraint(Collinear, line1, line2);
}

Constraint PointOnCircleConstraint(Point *point1, Circle *circle1) {
    return PointCircleConstraint(PointOnCircle, point1, circle1);
}

Constraint PointOnArcConstraint(Point *point1, Arc *arc1) {
    return PointArcConstraint(PointOnCircle, point1, arc1);
}

Constraint PointOnLineMidpointConstraint(Point *point1, Line *line1) {
    return PointLineConstraint(PointOnLineMidpoint, point1, line1);
}

Constraint PointOnArcMidpointConstraint(Point *point1, Arc *arc1) {
    return PointArcConstraint(PointOnArcMidpoint, point1, arc1);
}

Constraint PointOnCircleQuadConstraint(Point *point1, Circle *circle1,
                                       double *value) {
    Constraint con = {};
    con.type = PointOnCircleQuad;
    con.point1 = point1;
    con.circle1 = circle1;
    con.parameter = value;  // value only can be 0, 1, 2, 3, default 0
    return con;
}

Constraint SymmetricPointsConstraint(Point *point1, Point *point2, Line *sym) {
    Constraint con = {};
    con.type = SymmetricPoints;
    con.point1 = point1;
    con.point2 = point2;
    con.SymLine = sym;
    return con;
}

Constraint SymmetricLinesConstraint(Line *line1, Line *line2, Line *sym) {
    Constraint con = {};
    con.type = SymmetricLines;
    con.line1 = line1;
    con.line2 = line2;
    con.SymLine = sym;
    return con;
}

Constraint SymmetricCirclesConstraint(Circle *circle1, Circle *circle2,
                                      Line *sym) {
    Constraint con = {};
    con.type = SymmetricCircles;
    con.circle1 = circle1;
    con.circle2 = circle2;
    con.SymLine = sym;
    return con;
}

Constraint SymmetricArcsConstraint(Arc *arc1, Arc *arc2, Line *sym) {
    Constraint con = {};
    con.type = SymmetricArcs;
    con.arc1 = arc1;
    con.arc2 = arc2;
    con.SymLine = sym;
    return con;
}

Constraint LineInternalAngleConstraint(Line *line1, double *value) {
    return LineConstraint(LineInternalAngle, line1, value);
}

Constraint LineExternalAngleConstraint(Line *line1, double *value) {
    return LineConstraint(LineExternalAngle, line1, value);
}
