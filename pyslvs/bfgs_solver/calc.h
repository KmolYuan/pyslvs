#ifndef CLAC_H
#define CLAC_H

/*
 *  Created on: July 24, 2019
 *  Author: KmolYuan
 */

#include "solve.h"

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
    Collinear,
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

double calc(Constraint *, size_t);

#endif  // CLAC_H
