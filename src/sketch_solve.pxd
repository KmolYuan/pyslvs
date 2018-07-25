# -*- coding: utf-8 -*-
# cython: language_level=3

"""Wrapper header of BFGS algorithm."""

# __author__ = "Yuan Chang"
# __copyright__ = "Copyright (C) 2016-2018"
# __license__ = "AGPL"
# __email__ = "pyslvs@gmail.com"

cdef extern from "solve.h":
    
    cdef int Rough
    cdef int Fine
    cdef int MaxIterations
    
    cdef int Succsess
    cdef int NoSolution
    
    struct Point:
        double *x
        double *y
    
    struct Line:
        Point *p1
        Point *p2
    
    struct Arc:
        double *startAngle
        double *endAngle
        double *rad
        Point *center
    
    struct Circle:
        Point *center
        double *rad
    
    struct GeoConstraint:
        pass
    
    GeoConstraint PointOnPointConstraint(Point *, Point *);
    GeoConstraint P2PDistanceConstraint(Point *, Point *, double *)
    GeoConstraint P2PDistanceVertConstraint(Point *, Point *, double *)
    GeoConstraint P2PDistanceHorzConstraint(Point *, Point *, double *)
    GeoConstraint PointOnLineConstraint(Point *, Line *)
    GeoConstraint P2LDistanceConstraint(Point *, Line *, double *)
    GeoConstraint P2LDistanceVertConstraint(Point *, Line *, double *)
    GeoConstraint P2LDistanceHorzConstraint(Point *, Line *, double *)
    GeoConstraint VerticalConstraint(Line *)
    GeoConstraint HorizontalConstraint(Line *)
    GeoConstraint TangentToCircleConstraint(Line *, Circle *)
    GeoConstraint TangentToArcConstraint(Line *, Arc *)
    GeoConstraint ArcRulesConstraint(Arc *)
    GeoConstraint LineLengthConstraint(Line *)
    GeoConstraint EqualLegnthConstraint(Line *, Line *)
    GeoConstraint ArcRadiusConstraint(Arc *, double *)
    GeoConstraint EqualRadiusArcsConstraint(Arc *, Arc *)
    GeoConstraint EqualRadiusCirclesConstraint(Circle *, Circle *)
    GeoConstraint EqualRadiusCircArcConstraint(Circle *, Arc *)
    GeoConstraint ConcentricArcsConstraint(Arc *, Arc *)
    GeoConstraint ConcentricCirclesConstraint(Circle *, Circle *)
    GeoConstraint ConcentricCircArcConstraint(Circle *, Arc *)
    GeoConstraint CircleRadiusConstraint(Circle *, double *)
    GeoConstraint InternalAngleConstraint(Line *, Line *, double *)
    GeoConstraint ExternalAngleConstraint(Line *, Line *, double *)
    GeoConstraint PerpendicularConstraint(Line *, Line *)
    GeoConstraint ParallelConstraint(Line *, Line *)
    GeoConstraint ColinearConstraint(Line *, Line *)
    GeoConstraint PointOnCircleConstraint(Point *, Circle *)
    GeoConstraint PointOnArcConstraint(Point *, Arc *)
    GeoConstraint PointOnLineMidpointConstraint(Point *, Line *)
    GeoConstraint PointOnArcMidpointConstraint(Point *, Arc *)
    GeoConstraint PointOnCircleQuadConstraint(Point *, Circle *, double *)
    GeoConstraint SymmetricPointsConstraint(Point *, Point *, Line *)
    GeoConstraint SymmetricLinesConstraint(Line *, Line *, Line *)
    GeoConstraint SymmetricCirclesConstraint(Circle *, Circle *, Line *)
    GeoConstraint SymmetricArcsConstraint(Arc *, Arc *, Line *)
    GeoConstraint LineInternalAngleConstraint(Line *, double *)
    GeoConstraint LineExternalAngleConstraint(Line *, double *)
    
    int solve(double **, int, GeoConstraint *, int, int)
    void derivatives(double **, double *, int, GeoConstraint *, int)
