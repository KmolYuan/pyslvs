/*
 * Create on: Jun 6, 2018
 * Author: KmolYuan
 * These functions make writing constraint easier
 */

#include "solve.h"


Constraint PointOnPointConstraint(Point *point1, Point *point2) {
    Constraint con;
    con.type = Constraint::PointOnPoint;
    con.point1 = point1;
    con.point2 = point2;
    return con;
}


inline Constraint _P2PDistanceConstraint(
    Constraint::Types type,
    Point *point1,
    Point *point2,
    double *value
) {
    Constraint con;
    con.type = type;
    con.point1 = point1;
    con.point2 = point2;
    con.parameter = value;
    return con;
}


Constraint P2PDistanceConstraint(Point *point1, Point *point2, double *value) {
    return _P2PDistanceConstraint(Constraint::P2PDistance, point1, point2, value);
}


Constraint P2PDistanceVertConstraint(Point *point1, Point *point2, double *value) {
    return _P2PDistanceConstraint(Constraint::P2PDistanceVert, point1, point2, value);
}


Constraint P2PDistanceHorzConstraint(Point *point1, Point *point2, double *value) {
    return _P2PDistanceConstraint(Constraint::P2PDistanceHorz, point1, point2, value);
}


inline Constraint _PointLineConstraint(
    Constraint::Types type,
    Point *point1,
    Line *line1
) {
    Constraint con;
    con.type = type;
    con.point1 = point1;
    con.line1 = line1;
    return con;
}


Constraint PointOnLineConstraint(Point *point1, Line *line1) {
    return _PointLineConstraint(Constraint::PointOnLine, point1, line1);
}


inline Constraint _P2LDistanceConstraint(
    Constraint::Types type,
    Point *point1,
    Line *line1,
    double *value
) {
    Constraint con;
    con.type = type;
    con.point1 = point1;
    con.line1 = line1;
    con.parameter = value;
    return con;
}


Constraint P2LDistanceConstraint(Point *point1, Line *line1, double *value) {
    return _P2LDistanceConstraint(Constraint::P2LDistance, point1, line1, value);
}


Constraint P2LDistanceVertConstraint(Point *point1, Line *line1, double *value) {
    return _P2LDistanceConstraint(Constraint::P2LDistanceVert, point1, line1, value);
}


Constraint P2LDistanceHorzConstraint(Point *point1, Line *line1, double *value) {
    return _P2LDistanceConstraint(Constraint::P2LDistanceHorz, point1, line1, value);
}


inline Constraint _VerticalHorizontalConstraint(
    Constraint::Types type,
    Line *line1
) {
    Constraint con;
    con.type = type;
    con.line1 = line1;
    return con;
}


Constraint VerticalConstraint(Line *line1) {
    return _VerticalHorizontalConstraint(Constraint::Vertical, line1);
}


Constraint HorizontalConstraint(Line *line1) {
    return _VerticalHorizontalConstraint(Constraint::Horizontal, line1);
}


Constraint TangentToCircleConstraint(Line *line1, Circle *circle1) {
    Constraint con;
    con.type = Constraint::TangentToCircle;
    con.line1 = line1;
    con.circle1 = circle1;
    return con;
}


Constraint TangentToArcConstraint(Line *line1, Arc *arc1) {
    Constraint con;
    con.type = Constraint::TangentToCircle;
    con.line1 = line1;
    con.arc1 = arc1;
    return con;
}


Constraint ArcRulesConstraint(Arc *arc1) {
    Constraint con;
    con.type = Constraint::ArcRules;
    con.arc1 = arc1;
    return con;
}


inline Constraint _LineConstraint(
    Constraint::Types type,
    Line *line1,
    double *value
) {
    Constraint con;
    con.type = type;
    con.line1 = line1;
    con.parameter = value;
    return con;
}


Constraint LineLengthConstraint(Line *line1, double *value) {
    return _LineConstraint(Constraint::LineLength, line1, value);
}


inline Constraint _LinesConstraint(
    Constraint::Types type,
    Line *line1,
    Line *line2
) {
    Constraint con;
    con.type = type;
    con.line1 = line1;
    con.line2 = line2;
    return con;
}


Constraint EqualLegnthConstraint(Line *line1, Line *line2) {
    return _LinesConstraint(Constraint::EqualLegnth, line1, line2);
}


Constraint ArcRadiusConstraint(Arc *arc1, double *value) {
    Constraint con;
    con.type = Constraint::ArcRadius;
    con.arc1 = arc1;
    con.parameter = value;
    return con;
}


inline Constraint _ArcsConstraint(Constraint::Types type, Arc *arc1, Arc *arc2) {
    Constraint con;
    con.type = type;
    con.arc1 = arc1;
    con.arc2 = arc2;
    return con;
}


Constraint EqualRadiusArcsConstraint(Arc *arc1, Arc *arc2) {
    return _ArcsConstraint(Constraint::EqualRadiusArcs, arc1, arc2);
}


inline Constraint _CirclesConstraint(
    Constraint::Types type,
    Circle *circle1,
    Circle *circle2
) {
    Constraint con;
    con.type = type;
    con.circle1 = circle1;
    con.circle2 = circle2;
    return con;
}


inline Constraint _CircArcConstraint(
    Constraint::Types type,
    Circle *circle1,
    Arc *arc1
) {
    Constraint con;
    con.type = type;
    con.circle1 = circle1;
    con.arc1 = arc1;
    return con;
}


Constraint EqualRadiusCirclesConstraint(Circle *circle1, Circle *circle2) {
    return _CirclesConstraint(Constraint::EqualRadiusCircles, circle1, circle2);
}


Constraint EqualRadiusCircArcConstraint(Circle *circle1, Arc *arc1) {
    return _CircArcConstraint(Constraint::EqualRadiusCircles, circle1, arc1);
}


Constraint ConcentricArcsConstraint(Arc *arc1, Arc *arc2) {
    return _ArcsConstraint(Constraint::ConcentricArcs, arc1, arc2);
}


Constraint ConcentricCirclesConstraint(Circle *circle1, Circle *circle2) {
    return _CirclesConstraint(Constraint::ConcentricCircles, circle1, circle2);
}


Constraint ConcentricCircArcConstraint(Circle *circle1, Arc *arc1) {
    return _CircArcConstraint(Constraint::ConcentricCircArc, circle1, arc1);
}


Constraint CircleRadiusConstraint(Circle *circle1, double *value) {
    Constraint con;
    con.type = Constraint::CircleRadius;
    con.circle1 = circle1;
    con.parameter = value;
    return con;
}


inline Constraint _AngleConstraint(
    Constraint::Types type,
    Line *line1,
    Line *line2,
    double *value
) {
    Constraint con;
    con.type = type;
    con.line1 = line1;
    con.line2 = line2;
    con.parameter = value;
    return con;
}


Constraint InternalAngleConstraint(Line *line1, Line *line2, double *value) {
    return _AngleConstraint(Constraint::InternalAngle, line1, line2, value);
}


Constraint ExternalAngleConstraint(Line *line1, Line *line2, double *value) {
    return _AngleConstraint(Constraint::ExternalAngle, line1, line2, value);
}


Constraint PerpendicularConstraint(Line *line1, Line *line2) {
    return _LinesConstraint(Constraint::Perpendicular, line1, line2);
}


Constraint ParallelConstraint(Line *line1, Line *line2) {
    return _LinesConstraint(Constraint::Parallel, line1, line2);
}


Constraint ColinearConstraint(Line *line1, Line *line2) {
    return _LinesConstraint(Constraint::Colinear, line1, line2);
}


inline Constraint _PointCircleConstraint(
    Constraint::Types type,
    Point *point1,
    Circle *circle1
) {
    Constraint con;
    con.type = type;
    con.point1 = point1;
    con.circle1 = circle1;
    return con;
}


Constraint PointOnCircleConstraint(Point *point1, Circle *circle1) {
    return _PointCircleConstraint(Constraint::PointOnCircle, point1, circle1);
}


inline Constraint _PointArcConstraint(Constraint::Types type, Point *point1, Arc *arc1) {
    Constraint con;
    con.type = type;
    con.point1 = point1;
    con.arc1 = arc1;
    return con;
}


Constraint PointOnArcConstraint(Point *point1, Arc *arc1) {
    return _PointArcConstraint(Constraint::PointOnCircle, point1, arc1);
}


Constraint PointOnLineMidpointConstraint(Point *point1, Line *line1) {
    return _PointLineConstraint(Constraint::PointOnLineMidpoint, point1, line1);
}


Constraint PointOnArcMidpointConstraint(Point *point1, Arc *arc1) {
    return _PointArcConstraint(Constraint::PointOnArcMidpoint, point1, arc1);
}


Constraint PointOnCircleQuadConstraint(Point *point1, Circle *circle1, double *value) {
    Constraint con;
    con.type = Constraint::PointOnCircleQuad;
    con.point1 = point1;
    con.circle1 = circle1;
    con.parameter = value; //value only can be 0, 1, 2, 3, default 0.
    return con;
}


Constraint SymmetricPointsConstraint(Point *point1, Point *point2, Line *sym) {
    Constraint con;
    con.type = Constraint::SymmetricPoints;
    con.point1 = point1;
    con.point2 = point2;
    con.SymLine = sym;
    return con;
}


Constraint SymmetricLinesConstraint(Line *line1, Line *line2, Line *sym) {
    Constraint con;
    con.type = Constraint::SymmetricLines;
    con.line1 = line1;
    con.line2 = line2;
    con.SymLine = sym;
    return con;
}


Constraint SymmetricCirclesConstraint(Circle *circle1, Circle *circle2, Line *sym) {
    Constraint con;
    con.type = Constraint::SymmetricCircles;
    con.circle1 = circle1;
    con.circle2 = circle2;
    con.SymLine = sym;
    return con;
}


Constraint SymmetricArcsConstraint(Arc *arc1, Arc *arc2, Line *sym) {
    Constraint con;
    con.type = Constraint::SymmetricArcs;
    con.arc1 = arc1;
    con.arc2 = arc2;
    con.SymLine = sym;
    return con;
}


Constraint LineInternalAngleConstraint(Line *line1, double *value) {
    return _LineConstraint(Constraint::LineInternalAngle, line1, value);
}


Constraint LineExternalAngleConstraint(Line *line1, double *value) {
    return _LineConstraint(Constraint::LineExternalAngle, line1, value);
}
