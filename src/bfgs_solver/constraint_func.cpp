/*
 * Create on: Jun 6, 2018
 * Author: KmolYuan
 * These functions make writing constraint easier
 */

#include "solve.h"


GeoConstraint PointOnPointConstraint(Point *point1, Point *point2) {
    GeoConstraint con;
    con.type = GeoConstraint::PointOnPoint;
    con.point1 = point1;
    con.point2 = point2;
    return con;
}


inline GeoConstraint _P2PDistanceConstraint(
    GeoConstraint::Types type,
    Point *point1,
    Point *point2,
    double *value
) {
    GeoConstraint con;
    con.type = type;
    con.point1 = point1;
    con.point2 = point2;
    con.parameter = value;
    return con;
}


GeoConstraint P2PDistanceConstraint(Point *point1, Point *point2, double *value) {
    return _P2PDistanceConstraint(GeoConstraint::P2PDistance, point1, point2, value);
}


GeoConstraint P2PDistanceVertConstraint(Point *point1, Point *point2, double *value) {
    return _P2PDistanceConstraint(GeoConstraint::P2PDistanceVert, point1, point2, value);
}


GeoConstraint P2PDistanceHorzConstraint(Point *point1, Point *point2, double *value) {
    return _P2PDistanceConstraint(GeoConstraint::P2PDistanceHorz, point1, point2, value);
}


inline GeoConstraint _PointLineConstraint(
    GeoConstraint::Types type,
    Point *point1,
    Line *line1
) {
    GeoConstraint con;
    con.type = type;
    con.point1 = point1;
    con.line1 = line1;
    return con;
}


GeoConstraint PointOnLineConstraint(Point *point1, Line *line1) {
    return _PointLineConstraint(GeoConstraint::PointOnLine, point1, line1);
}


inline GeoConstraint _P2LDistanceConstraint(
    GeoConstraint::Types type,
    Point *point1,
    Line *line1,
    double *value
) {
    GeoConstraint con;
    con.type = type;
    con.point1 = point1;
    con.line1 = line1;
    con.parameter = value;
    return con;
}


GeoConstraint P2LDistanceConstraint(Point *point1, Line *line1, double *value) {
    return _P2LDistanceConstraint(GeoConstraint::P2LDistance, point1, line1, value);
}


GeoConstraint P2LDistanceVertConstraint(Point *point1, Line *line1, double *value) {
    return _P2LDistanceConstraint(GeoConstraint::P2LDistanceVert, point1, line1, value);
}


GeoConstraint P2LDistanceHorzConstraint(Point *point1, Line *line1, double *value) {
    return _P2LDistanceConstraint(GeoConstraint::P2LDistanceHorz, point1, line1, value);
}


inline GeoConstraint _VerticalHorizontalConstraint(
    GeoConstraint::Types type,
    Line *line1
) {
    GeoConstraint con;
    con.type = type;
    con.line1 = line1;
    return con;
}


GeoConstraint VerticalConstraint(Line *line1) {
    return _VerticalHorizontalConstraint(GeoConstraint::Vertical, line1);
}


GeoConstraint HorizontalConstraint(Line *line1) {
    return _VerticalHorizontalConstraint(GeoConstraint::Horizontal, line1);
}


GeoConstraint TangentToCircleConstraint(Line *line1, Circle *circle1) {
    GeoConstraint con;
    con.type = GeoConstraint::TangentToCircle;
    con.line1 = line1;
    con.circle1 = circle1;
    return con;
}


GeoConstraint TangentToArcConstraint(Line *line1, Arc *arc1) {
    GeoConstraint con;
    con.type = GeoConstraint::TangentToCircle;
    con.line1 = line1;
    con.arc1 = arc1;
    return con;
}


GeoConstraint ArcRulesConstraint(Arc *arc1) {
    GeoConstraint con;
    con.type = GeoConstraint::ArcRules;
    con.arc1 = arc1;
    return con;
}


inline GeoConstraint _LineConstraint(
    GeoConstraint::Types type,
    Line *line1,
    double *value
) {
    GeoConstraint con;
    con.type = type;
    con.line1 = line1;
    con.parameter = value;
    return con;
}


GeoConstraint LineLengthConstraint(Line *line1, double *value) {
    return _LineConstraint(GeoConstraint::LineLength, line1, value);
}


inline GeoConstraint _LinesConstraint(
    GeoConstraint::Types type,
    Line *line1,
    Line *line2
) {
    GeoConstraint con;
    con.type = type;
    con.line1 = line1;
    con.line2 = line2;
    return con;
}


GeoConstraint EqualLegnthConstraint(Line *line1, Line *line2) {
    return _LinesConstraint(GeoConstraint::EqualLegnth, line1, line2);
}


GeoConstraint ArcRadiusConstraint(Arc *arc1, double *value) {
    GeoConstraint con;
    con.type = GeoConstraint::ArcRadius;
    con.arc1 = arc1;
    con.parameter = value;
    return con;
}


inline GeoConstraint _ArcsConstraint(GeoConstraint::Types type, Arc *arc1, Arc *arc2) {
    GeoConstraint con;
    con.type = type;
    con.arc1 = arc1;
    con.arc2 = arc2;
    return con;
}


GeoConstraint EqualRadiusArcsConstraint(Arc *arc1, Arc *arc2) {
    return _ArcsConstraint(GeoConstraint::EqualRadiusArcs, arc1, arc2);
}


inline GeoConstraint _CirclesConstraint(
    GeoConstraint::Types type,
    Circle *circle1,
    Circle *circle2
) {
    GeoConstraint con;
    con.type = type;
    con.circle1 = circle1;
    con.circle2 = circle2;
    return con;
}


inline GeoConstraint _CircArcConstraint(
    GeoConstraint::Types type,
    Circle *circle1,
    Arc *arc1
) {
    GeoConstraint con;
    con.type = type;
    con.circle1 = circle1;
    con.arc1 = arc1;
    return con;
}


GeoConstraint EqualRadiusCirclesConstraint(Circle *circle1, Circle *circle2) {
    return _CirclesConstraint(GeoConstraint::EqualRadiusCircles, circle1, circle2);
}


GeoConstraint EqualRadiusCircArcConstraint(Circle *circle1, Arc *arc1) {
    return _CircArcConstraint(GeoConstraint::EqualRadiusCircles, circle1, arc1);
}


GeoConstraint ConcentricArcsConstraint(Arc *arc1, Arc *arc2) {
    return _ArcsConstraint(GeoConstraint::ConcentricArcs, arc1, arc2);
}


GeoConstraint ConcentricCirclesConstraint(Circle *circle1, Circle *circle2) {
    return _CirclesConstraint(GeoConstraint::ConcentricCircles, circle1, circle2);
}


GeoConstraint ConcentricCircArcConstraint(Circle *circle1, Arc *arc1) {
    return _CircArcConstraint(GeoConstraint::ConcentricCircArc, circle1, arc1);
}


GeoConstraint CircleRadiusConstraint(Circle *circle1, double *value) {
    GeoConstraint con;
    con.type = GeoConstraint::CircleRadius;
    con.circle1 = circle1;
    con.parameter = value;
    return con;
}


inline GeoConstraint _AngleConstraint(
    GeoConstraint::Types type,
    Line *line1,
    Line *line2,
    double *value
) {
    GeoConstraint con;
    con.type = type;
    con.line1 = line1;
    con.line2 = line2;
    con.parameter = value;
    return con;
}


GeoConstraint InternalAngleConstraint(Line *line1, Line *line2, double *value) {
    return _AngleConstraint(GeoConstraint::InternalAngle, line1, line2, value);
}


GeoConstraint ExternalAngleConstraint(Line *line1, Line *line2, double *value) {
    return _AngleConstraint(GeoConstraint::ExternalAngle, line1, line2, value);
}


GeoConstraint PerpendicularConstraint(Line *line1, Line *line2) {
    return _LinesConstraint(GeoConstraint::Perpendicular, line1, line2);
}


GeoConstraint ParallelConstraint(Line *line1, Line *line2) {
    return _LinesConstraint(GeoConstraint::Parallel, line1, line2);
}


GeoConstraint ColinearConstraint(Line *line1, Line *line2) {
    return _LinesConstraint(GeoConstraint::Colinear, line1, line2);
}


inline GeoConstraint _PointCircleConstraint(
    GeoConstraint::Types type,
    Point *point1,
    Circle *circle1
) {
    GeoConstraint con;
    con.type = type;
    con.point1 = point1;
    con.circle1 = circle1;
    return con;
}


GeoConstraint PointOnCircleConstraint(Point *point1, Circle *circle1) {
    return _PointCircleConstraint(GeoConstraint::PointOnCircle, point1, circle1);
}


inline GeoConstraint _PointArcConstraint(GeoConstraint::Types type, Point *point1, Arc *arc1) {
    GeoConstraint con;
    con.type = type;
    con.point1 = point1;
    con.arc1 = arc1;
    return con;
}


GeoConstraint PointOnArcConstraint(Point *point1, Arc *arc1) {
    return _PointArcConstraint(GeoConstraint::PointOnCircle, point1, arc1);
}


GeoConstraint PointOnLineMidpointConstraint(Point *point1, Line *line1) {
    return _PointLineConstraint(GeoConstraint::PointOnLineMidpoint, point1, line1);
}


GeoConstraint PointOnArcMidpointConstraint(Point *point1, Arc *arc1) {
    return _PointArcConstraint(GeoConstraint::PointOnArcMidpoint, point1, arc1);
}


GeoConstraint PointOnCircleQuadConstraint(Point *point1, Circle *circle1, double *value) {
    GeoConstraint con;
    con.type = GeoConstraint::PointOnCircleQuad;
    con.point1 = point1;
    con.circle1 = circle1;
    con.parameter = value; //value only can be 0, 1, 2, 3, default 0.
    return con;
}


GeoConstraint SymmetricPointsConstraint(Point *point1, Point *point2, Line *sym) {
    GeoConstraint con;
    con.type = GeoConstraint::SymmetricPoints;
    con.point1 = point1;
    con.point2 = point2;
    con.SymLine = sym;
    return con;
}


GeoConstraint SymmetricLinesConstraint(Line *line1, Line *line2, Line *sym) {
    GeoConstraint con;
    con.type = GeoConstraint::SymmetricLines;
    con.line1 = line1;
    con.line2 = line2;
    con.SymLine = sym;
    return con;
}


GeoConstraint SymmetricCirclesConstraint(Circle *circle1, Circle *circle2, Line *sym) {
    GeoConstraint con;
    con.type = GeoConstraint::SymmetricCircles;
    con.circle1 = circle1;
    con.circle2 = circle2;
    con.SymLine = sym;
    return con;
}


GeoConstraint SymmetricArcsConstraint(Arc *arc1, Arc *arc2, Line *sym) {
    GeoConstraint con;
    con.type = GeoConstraint::SymmetricArcs;
    con.arc1 = arc1;
    con.arc2 = arc2;
    con.SymLine = sym;
    return con;
}


GeoConstraint LineInternalAngleConstraint(Line *line1, double *value) {
    return _LineConstraint(GeoConstraint::LineInternalAngle, line1, value);
}


GeoConstraint LineExternalAngleConstraint(Line *line1, double *value) {
    return _LineConstraint(GeoConstraint::LineExternalAngle, line1, value);
}
