/*
 *  Created on: July 24, 2019
 *  Author: KmolYuan
 */

#include "calc.h"
#include <cmath>
#ifdef DEBUG
#include <iostream>
using namespace std;
#endif

///////////////////////////////////////////////////////////////////////
/// constraint defines (these make writing constraint equations easier)
///////////////////////////////////////////////////////////////////////

#define P1 con->point1
#define P1_x (*P1->x)
#define P1_y (*P1->y)
#define P2 con->point2
#define P2_x (*P2->x)
#define P2_y (*P2->y)
#define L1 con->line1
#define L1_P1 L1->p1
#define L1_P1_x (*L1_P1->x)
#define L1_P1_y (*L1_P1->y)
#define L1_P2 L1->p2
#define L1_P2_x (*L1_P2->x)
#define L1_P2_y (*L1_P2->y)
#define L2 con->line2
#define L2_P1 L2->p1
#define L2_P1_x (*L2_P1->x)
#define L2_P1_y (*L2_P1->y)
#define L2_P2 L2->p2
#define L2_P2_x (*L2_P2->x)
#define L2_P2_y (*L2_P2->y)
#define C1 con->circle1
#define C1_Center C1->center
#define C1_Center_x (*C1_Center->x)
#define C1_Center_y (*C1_Center->y)
#define C1_rad (*C1->rad)
#define C2 con->circle2
#define C2_Center C2->center
#define C2_Center_x (*C2_Center->x)
#define C2_Center_y (*C2_Center->y)
#define C2_rad (*C2->rad)
#define A1 con->arc1
#define A1_startA (*A1->startAngle)
#define A1_endA (*A1->endAngle)
#define A1_radius (*A1->rad)
#define A1_Center A1->center
#define A1_Center_x (*A1_Center->x)
#define A1_Center_y (*A1_Center->y)
#define A2 con->arc2
#define A2_startA (*A2->startAngle)
#define A2_endA (*A2->endAngle)
#define A2_radius (*A2->rad)
#define A2_Center A2->center
#define A2_Center_x (*A2_Center->x)
#define A2_Center_y (*A2_Center->y)
#define A1_Start_x (A1_Center_x + A1_radius * cos(A1_startA))
#define A1_Start_y (A1_Center_y + A1_radius * sin(A1_startA))
#define A1_End_x (A1_Center_x + A1_radius * cos(A1_endA))
#define A1_End_y (A1_Center_y + A1_radius * sin(A1_endA))
#define A2_Start_x (A1_Center_x + A2_radius * cos(A2_startA))
#define A2_Start_y (A1_Center_y + A2_radius * sin(A2_startA))
#define A2_End_x (A1_Center_x + A2_radius * cos(A2_endA))
#define A2_End_y (A1_Center_y + A2_radius * sin(A2_endA))
#define RADIUS (*con->parameter)
#define DISTANCE fabs(RADIUS)
#define QUAD (int)RADIUS
#define SYM con->SymLine
#define SYM_P1 SYM->p1
#define SYM_P1_x (*SYM_P1->x)
#define SYM_P1_y (*SYM_P1->y)
#define SYM_P2 SYM->p2
#define SYM_P2_x (*SYM_P2->x)
#define SYM_P2_y (*SYM_P2->y)

namespace {
double point_on_point(Constraint *con) {
    // Hopefully avoid this constraint,
    // make coincident points use the same parameters
    double dx = P2_x - P1_x;
    double dy = P2_y - P1_y;
    return dx * dx + dy * dy;
}

double p2p_distance(Constraint *con) {
    double temp = hypot(P2_x - P1_x, P2_y - P1_y) - DISTANCE;
    return temp * temp * 100;
}

double p2p_distance_vert(Constraint *con) {
    double temp = hypot(0, P2_y - P1_y) - DISTANCE;
    return temp * temp * 100;
}

double p2p_distance_horz(Constraint *con) {
    double temp = hypot(P2_x - P1_x, 0) - DISTANCE;
    return temp * temp * 100;
}

double point_on_line(Constraint *con) {
    double dx = L1_P2_x - L1_P1_x;
    double dy = L1_P2_y - L1_P1_y;

    double t =
        ((P1_x - L1_P1_x) * dx + (P1_y - L1_P1_y) * dy) / (dx * dx + dy * dy);
    double temp = hypot((L1_P1_x + dx * t) - P1_x, (L1_P1_y + dy * t) - P1_y);
    return temp * temp / 1000;
}

double p2l_distance(Constraint *con) {
    double dx = L1_P2_x - L1_P1_x;
    double dy = L1_P2_y - L1_P1_y;

    double t = -(L1_P1_x * dx - P1_x * dx + L1_P1_y * dy - P1_y * dy) /
               (dx * dx + dy * dy);
    double temp =
        hypot(P1_x - (L1_P1_x + dx * t), P1_y - (L1_P1_y + dy * t)) - DISTANCE;
    return temp * temp / 100;
}

double p2l_distance_vert(Constraint *con) {
    double dx = L1_P2_x - L1_P1_x;
    double dy = L1_P2_y - L1_P1_y;

    double t = (P1_x - L1_P1_x) / dx;
    double temp = fabs(P1_y - (L1_P1_y + dy * t)) - DISTANCE;
    return temp * temp / 100;
}

double p2l_distance_horz(Constraint *con) {
    double dx = L1_P2_x - L1_P1_x;
    double dy = L1_P2_y - L1_P1_y;

    double t = (P1_y - L1_P1_y) / dy;
    double temp = fabs(P1_x - (L1_P1_x + dx * t)) - DISTANCE;
    return temp * temp / 100;
}

double line_length(Constraint *con) {
    double temp = hypot(L1_P2_x - L1_P1_x, L1_P2_y - L1_P1_y) - DISTANCE;
    return temp * temp * 100;
}

double equal_length(Constraint *con) {
    double temp = hypot(L1_P2_x - L1_P1_x, L1_P2_y - L1_P1_y) -
                  hypot(L2_P2_x - L2_P1_x, L2_P2_y - L2_P1_y);
    return temp * temp * 100;
}

double vertical(Constraint *con) {
    double dx = L1_P2_x - L1_P1_x;
    return dx * dx * 100;
}

double horizontal(Constraint *con) {
    double dy = L1_P2_y - L1_P1_y;
    return dy * dy * 100;
}

double tangent_to_circle(Constraint *con) {
    double dx = L1_P2_x - L1_P1_x;
    double dy = L1_P2_y - L1_P1_y;
    double hyp1 = hypot(dx, dy);
    // Calculate the expected tangent intersection points
    double temp = (-dy * (C1_Center_x - dy / hyp1 * C1_rad) +
                   dx * (C1_Center_y + dx / hyp1 * C1_rad) +
                   (L1_P1_x * L1_P2_y - L1_P2_x * L1_P1_y)) /
                  hyp1;
    double temp2 = (-dy * (C1_Center_x + dy / hyp1 * C1_rad) +
                    dx * (C1_Center_y - dx / hyp1 * C1_rad) +
                    (L1_P1_x * L1_P2_y - L1_P2_x * L1_P1_y)) /
                   hyp1;
    return (temp < temp2) ? temp * temp : temp2 * temp2;
}

double tangent_to_arc(Constraint *con) {
    double dx = L1_P2_x - L1_P1_x;
    double dy = L1_P2_y - L1_P1_y;
    double t =
        -(L1_P1_x * dx - A1_Center_x * dx + L1_P1_y * dy - A1_Center_y * dy) /
        (dx * dx + dy * dy);
    double Ex = A1_Center_x - (L1_P1_x + dx * t);
    double Ey = A1_Center_y - (L1_P1_y + dy * t);
    double temp = Ex * Ex + Ey * Ey -
                  (A1_Center_x - A1_Start_x) * (A1_Center_x - A1_Start_x) -
                  (A1_Center_y - A1_Start_y) * (A1_Center_y - A1_Start_y);
    return temp * temp;
}

double arc_rules(Constraint *con) {
    double dx = A1_End_x * A1_End_x;
    double dy = A1_End_y * A1_End_y;
    double rad1 = A1_Start_x * A1_Start_x;
    double rad2 = A1_Start_y * A1_Start_y;
    double temp = -2 * A1_Center_x * A1_End_x + dx -
                  2 * A1_Center_y * A1_End_y + dy +
                  2 * A1_Center_x * A1_Start_x - rad1 +
                  2 * A1_Center_y * A1_Start_y - rad2;
    return temp * temp /
           (4. * dx + dy - 2 * A1_End_x * A1_Start_x + rad1 -
            2 * A1_End_y * A1_Start_y + rad2);
}

double arc_radius(Constraint *con) {
    double rad1 = hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
    double rad2 = hypot(A1_Center_x - A1_End_x, A1_Center_y - A1_End_y);
    double temp = rad1 - rad2 - RADIUS;
    return temp * temp;
}

double equal_radius_arcs(Constraint *con) {
    double temp = hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y) -
                  hypot(A2_Center_x - A2_Start_x, A2_Center_y - A2_Start_y);
    return temp * temp;
}

double equal_radius_circles(Constraint *con) {
    double temp = C1_rad - C2_rad;
    return temp * temp;
}

double equal_radius_circ_arc(Constraint *con) {
    double rad1 = hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
    double temp = rad1 - C1_rad;
    return temp * temp;
}

double concentric_arcs(Constraint *con) {
    double temp = hypot(A1_Center_x - A2_Center_x, A1_Center_y - A2_Center_y);
    return temp * temp;
}

double concentric_circles(Constraint *con) {
    double temp = hypot(C1_Center_x - C2_Center_x, C1_Center_y - C2_Center_y);
    return temp * temp;
}

double concentric_circle_arc(Constraint *con) {
    double temp = hypot(A1_Center_x - C1_Center_x, A1_Center_y - C1_Center_y);
    return temp * temp;
}

double circle_radius(Constraint *con) {
    return (C1_rad - RADIUS) * (C1_rad - RADIUS);
}

double internal_angle(Constraint *con) {
    double temp = atan2(L2_P2_y - L2_P1_y, L2_P2_x - L2_P1_x) -
                  atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) - RADIUS;
    return temp * temp;
}

double external_angle(Constraint *con) {
    double temp = atan2(L2_P2_y - L2_P1_y, L2_P2_x - L2_P1_x) -
                  atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) - (M_PI - RADIUS);
    return temp * temp;
}

double line_internal_angle(Constraint *con) {
    double temp = sin(atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) - RADIUS);
    return temp * temp;
}

double line_external_angle(Constraint *con) {
    double temp =
        sin(atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) - (M_PI - RADIUS));
    return temp * temp;
}

double perpendicular(Constraint *con) {
    double dx = L1_P2_x - L1_P1_x;
    double dy = L1_P2_y - L1_P1_y;
    double dx2 = L2_P2_x - L2_P1_x;
    double dy2 = L2_P2_y - L2_P1_y;

    double hyp1 = hypot(dx, dy);
    double hyp2 = hypot(dx2, dy2);

    dx /= hyp1;
    dy /= hyp1;
    dx2 /= hyp2;
    dy2 /= hyp2;

    double temp = dx * dx2 + dy * dy2;
    return temp * temp;
}

double parallel(Constraint *con) {
    double dx = L1_P2_x - L1_P1_x;
    double dy = L1_P2_y - L1_P1_y;
    double dx2 = L2_P2_x - L2_P1_x;
    double dy2 = L2_P2_y - L2_P1_y;

    double hyp1 = hypot(dx, dy);
    double hyp2 = hypot(dx2, dy2);

    dx /= hyp1;
    dy /= hyp1;
    dx2 /= hyp2;
    dy2 /= hyp2;

    double temp = dy * dx2 - dx * dy2;
    return temp * temp;
}

double colinear(Constraint *con) {
    double dx = L1_P2_x - L1_P1_x;
    double dy = L1_P2_y - L1_P1_y;

    double m = dy / dx;
    double n = dx / dy;
    // Calculate the error between the expected intersection point
    // and the true point of the second lines two end points on the
    // first line
    if (m <= 1 && m > -1) {
        // Calculate the expected y point given the x coordinate of the point
        double e1 = L1_P1_y + m * (L2_P1_x - L1_P1_x);
        double e2 = L1_P1_y + m * (L2_P2_x - L1_P1_x);
        return (e1 - L2_P1_y) * (e1 - L2_P1_y) +
               (e2 - L2_P2_y) * (e2 - L2_P2_y);
    } else {
        // Calculate the expected x point given the y coordinate of the point
        double e1 = L1_P1_x + n * (L2_P1_y - L1_P1_y);
        double e2 = L1_P1_x + n * (L2_P2_y - L1_P1_y);
        return (e1 - L2_P1_x) * (e1 - L2_P1_x) +
               (e2 - L2_P2_x) * (e2 - L2_P2_x);
    }
}

double point_on_circle(Constraint *con) {
    // see what the current RADIUS to the point is
    double rad1 = hypot(C1_Center_x - P1_x, C1_Center_y - P1_y);
    // Compare this RADIUS to the RADIUS of the circle, return the error squared
    double temp = rad1 - C1_rad;
    return temp * temp;
}

double point_on_arc(Constraint *con) {
    // see what the current RADIUS to the point is
    double rad1 = hypot(A1_Center_x - P1_x, A1_Center_y - P1_y);
    double rad2 = hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
    // Compare this RADIUS to the RADIUS of the circle, return the error squared
    double temp = rad1 - rad2;
    return temp * temp;
}

double point_on_line_midpoint(Constraint *con) {
    double Ex = (L1_P1_x + L1_P2_x) / 2 - P1_x;
    double Ey = (L1_P1_y + L1_P2_y) / 2 - P1_y;
    return Ex * Ex + Ey * Ey;
}

double point_on_arc_midpoint(Constraint *con) {
    double rad1 = hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
    double temp = atan2(A1_Start_y - A1_Center_y, A1_Start_x - A1_Center_x);
    double temp2 = atan2(A1_End_y - A1_Center_y, A1_End_x - A1_Center_x);
    double Ex = A1_Center_x + rad1 * cos((temp2 + temp) / 2);
    double Ey = A1_Center_y + rad1 * sin((temp2 + temp) / 2);
    temp = Ex - P1_x;
    temp2 = Ey - P1_y;
    return temp * temp + temp2 * temp2;
}

double point_on_circle_quad(Constraint *con) {
    double Ex = C1_Center_x;
    double Ey = C1_Center_y;
    switch (QUAD) {
        case 1:
            Ey += C1_rad;
            break;
        case 2:
            Ex -= C1_rad;
            break;
        case 3:
            Ey -= C1_rad;
            break;
        case 4:
        default:
            Ex += C1_rad;
            break;
    }
    double temp = Ex - P1_x;
    double temp2 = Ey - P1_y;
    return temp * temp + temp2 * temp2;
}

double symmetric_points(Constraint *con) {
    double dx = SYM_P2_x - SYM_P1_x;
    double dy = SYM_P2_y - SYM_P1_y;
    double t = -(dy * P1_x - dx * P1_y - dy * SYM_P1_x + dx * SYM_P1_y) /
               (dx * dx + dy * dy);
    double temp = P1_x + dy * t * 2 - P2_x;
    double temp2 = P1_y - dx * t * 2 - P2_y;
    return temp * temp + temp2 * temp2;
}

double symmetric_lines(Constraint *con) {
    double dx = SYM_P2_x - SYM_P1_x;
    double dy = SYM_P2_y - SYM_P1_y;
    double hyp1 = dx * dx + dy * dy;
    double m = -dy * SYM_P1_x + dx * SYM_P1_y;

    double t = -(dy * L1_P1_x - dx * L1_P1_y + m) / hyp1;
    double Ex = L1_P1_x + dy * t * 2;
    double Ey = L1_P1_y - dx * t * 2;
    double temp = Ex - L2_P1_x;
    double temp2 = Ey - L2_P1_y;
    double error = temp * temp + temp2 * temp2;
    t = -(dy * L1_P2_x - dx * L1_P2_y + m) / hyp1;
    Ex = L1_P2_x + dy * t * 2;
    Ey = L1_P2_y - dx * t * 2;
    temp = Ex - L2_P2_x;
    temp2 = Ey - L2_P2_y;
    return error + temp * temp + temp2 * temp2;
}

double symmetric_circles(Constraint *con) {
    double dx = SYM_P2_x - SYM_P1_x;
    double dy = SYM_P2_y - SYM_P1_y;
    double t =
        -(dy * C1_Center_x - dx * C1_Center_y - dy * SYM_P1_x + dx * SYM_P1_y) /
        (dx * dx + dy * dy);
    double Ex = C1_Center_x + dy * t * 2;
    double Ey = C1_Center_y - dx * t * 2;
    double temp = Ex - C2_Center_x;
    double temp2 = Ey - C2_Center_y;
    double error = temp * temp + temp2 * temp2;
    temp = C1_rad - C2_rad;
    return error + temp * temp;
}

double symmetric_arcs(Constraint *con) {
    double dx = SYM_P2_x - SYM_P1_x;
    double dy = SYM_P2_y - SYM_P1_y;
    double hyp1 = dx * dx + dy * dy;
    double m = -dy * SYM_P1_x + dx * SYM_P1_y;

    double t = -(dy * A1_Start_x - dx * A1_Start_y + m) / hyp1;
    double Ex = A1_Start_x + dy * t * 2;
    double Ey = A1_Start_y - dx * t * 2;
    double temp = Ex - A2_Start_x;
    double temp2 = Ey - A2_Start_y;
    double error = temp * temp + temp2 * temp2;
    t = -(dy * A1_End_x - dx * A1_End_y + m) / hyp1;
    Ex = A1_End_x + dy * t * 2;
    Ey = A1_End_y - dx * t * 2;
    temp = Ex - A2_End_x;
    temp2 = Ey - A2_End_y;
    error += temp * temp + temp2 * temp2;
    t = -(dy * A1_Center_x - dx * A1_Center_y + m) / hyp1;
    Ex = A1_Center_x + dy * t * 2;
    Ey = A1_Center_y - dx * t * 2;
    temp = Ex - A2_Center_x;
    temp2 = Ey - A2_Center_y;
    return error + temp * temp + temp2 * temp2;
}
}  // namespace

double calc(Constraint *cons, size_t cons_len) {
    double error = 0;
    for (size_t i = 0; i < cons_len; i++) {
        Constraint *con = cons + i;
        switch (con->type) {
            case PointOnPoint:
                error += point_on_point(con);
                break;
            case P2PDistance:
                error += p2p_distance(con);
                break;
            case P2PDistanceVert:
                error += p2p_distance_vert(con);
                break;
            case P2PDistanceHorz:
                error += p2p_distance_horz(con);
                break;
            case PointOnLine:
                error += point_on_line(con);
                break;
            case P2LDistance:
                error += p2l_distance(con);
                break;
            case P2LDistanceVert:
                error += p2l_distance_vert(con);
                break;
            case P2LDistanceHorz:
                error += p2l_distance_horz(con);
                break;
            case LineLength:
                error += line_length(con);
                break;
            case EqualLegnth:
                error += equal_length(con);
                break;
            case Vertical:
                error += vertical(con);
                break;
            case Horizontal:
                error += horizontal(con);
                break;
            case TangentToCircle:
                error += tangent_to_circle(con);
                break;
            case TangentToArc:
                error += tangent_to_arc(con);
                break;
            case ArcRules:
                error += arc_rules(con);
                break;
            case ArcRadius:
                error += arc_radius(con);
                break;
            case EqualRadiusArcs:
                error += equal_radius_arcs(con);
                break;
            case EqualRadiusCircles:
                error += equal_radius_circles(con);
                break;
            case EqualRadiusCircArc:
                error += equal_radius_circ_arc(con);
                break;
            case ConcentricArcs:
                error += concentric_arcs(con);
                break;
            case ConcentricCircles:
                error += concentric_circles(con);
                break;
            case ConcentricCircArc:
                error += concentric_circle_arc(con);
                break;
            case CircleRadius:
                error += circle_radius(con);
                break;
            case InternalAngle:
                error += internal_angle(con);
                break;
            case ExternalAngle:
                error += external_angle(con);
                break;
            case LineInternalAngle:
                error += line_internal_angle(con);
                break;
            case LineExternalAngle:
                error += line_external_angle(con);
                break;
            case Perpendicular:
                error += perpendicular(con);
                break;
            case Parallel:
                error += parallel(con);
                break;
            case Collinear:
                error += colinear(con);
                break;
            case PointOnCircle:
                error += point_on_circle(con);
                break;
            case PointOnArc:
                error += point_on_arc(con);
                break;
            case PointOnLineMidpoint:
                error += point_on_line_midpoint(con);
                break;
            case PointOnArcMidpoint:
                error += point_on_arc_midpoint(con);
                break;
            case PointOnCircleQuad:
                error += point_on_circle_quad(con);
                break;
            case SymmetricPoints:
                error += symmetric_points(con);
                break;
            case SymmetricLines:
                error += symmetric_lines(con);
                break;
            case SymmetricCircles:
                error += symmetric_circles(con);
                break;
            case SymmetricArcs:
                error += symmetric_arcs(con);
                break;
        }
    }
#ifdef DEBUG
    cout << "Error: " << error << endl;
#endif
    return error;
}
