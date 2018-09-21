/*
 *  Created on: July 24, 2018
 *  Author: KmolYuan
 */

#include <cmath>
#include "calc.h"
#ifdef DEBUG
#include <iostream>
using namespace std;
#endif

///////////////////////////////////////////////////////////////////////
/// constraint defines (these make writing constraint equations easier)
///////////////////////////////////////////////////////////////////////

#define CON_i        cons[i]
#define P1           CON_i.point1
#define P1_x         *P1->x
#define P1_y         *P1->y
#define P2           CON_i.point2
#define P2_x         *P2->x
#define P2_y         *P2->y
#define L1           CON_i.line1
#define L1_P1        L1->p1
#define L1_P1_x      *L1_P1->x
#define L1_P1_y      *L1_P1->y
#define L1_P2        L1->p2
#define L1_P2_x      *L1_P2->x
#define L1_P2_y      *L1_P2->y
#define L2           CON_i.line2
#define L2_P1        L2->p1
#define L2_P1_x      *L2_P1->x
#define L2_P1_y      *L2_P1->y
#define L2_P2        L2->p2
#define L2_P2_x      *L2_P2->x
#define L2_P2_y      *L2_P2->y
#define C1           CON_i.circle1
#define C1_Center    C1->center
#define C1_Center_x  *C1_Center->x
#define C1_Center_y  *C1_Center->y
#define C1_rad       *C1->rad
#define C2           CON_i.circle2
#define C2_Center    C2->center
#define C2_Center_x  *C2_Center->x
#define C2_Center_y  *C2_Center->y
#define C2_rad       *C2->rad
#define A1           CON_i.arc1
#define A1_startA    *A1->startAngle
#define A1_endA      *A1->endAngle
#define A1_radius    *A1->rad
#define A1_Center    A1->center
#define A1_Center_x  *A1_Center->x
#define A1_Center_y  *A1_Center->y
#define A2           CON_i.arc2
#define A2_startA    *A2->startAngle
#define A2_endA      *A2->endAngle
#define A2_radius    *A2->rad
#define A2_Center    A2->center
#define A2_Center_x  *A2_Center->x
#define A2_Center_y  *A2_Center->y
#define A1_Start_x   (A1_Center_x + A1_radius * cos(A1_startA))
#define A1_Start_y   (A1_Center_y + A1_radius * sin(A1_startA))
#define A1_End_x     (A1_Center_x + A1_radius * cos(A1_endA))
#define A1_End_y     (A1_Center_y + A1_radius * sin(A1_endA))
#define A2_Start_x   (A1_Center_x + A2_radius * cos(A2_startA))
#define A2_Start_y   (A1_Center_y + A2_radius * sin(A2_startA))
#define A2_End_x     (A1_Center_x + A2_radius * cos(A2_endA))
#define A2_End_y     (A1_Center_y + A2_radius * sin(A2_endA))
#define PARAMETER    *CON_i.parameter
#define distance     fabs(PARAMETER)
#define radius       PARAMETER
#define angleP       PARAMETER
#define quadIndex    (int)PARAMETER
#define Sym          CON_i.SymLine
#define Sym_P1       Sym->p1
#define Sym_P1_x     *Sym_P1->x
#define Sym_P1_y     *Sym_P1->y
#define Sym_P2       Sym->p2
#define Sym_P2_x     *Sym_P2->x
#define Sym_P2_y     *Sym_P2->y


double calc(Constraint *cons, const int consLength) {
    double error = 0;
    double temp, temp2, dx, dy, m, n, Ex, Ey, rad1, rad2, t, dx2, dy2, hyp1, hyp2;
    for(int i = 0; i < consLength; i++) {
        switch((int)cons[i].type) {

        case Constraint::PointOnPoint:
            // Hopefully avoid this constraint,
            // make coincident points use the same parameters
            dx = P1_x - P2_x;
            dy = P1_y - P2_y;
            error += dx * dx + dy * dy;
            break;

        case Constraint::P2PDistance:
            temp = _hypot(P2_x - P1_x, P2_y - P1_y) - distance;
            error += temp * temp * 100;
            break;

        case Constraint::P2PDistanceVert:
            temp = _hypot(0, P2_y - P1_y) - distance;
            error += temp * temp * 100;
            break;

        case Constraint::P2PDistanceHorz:
            temp = _hypot(P2_x - P1_x, 0) - distance;
            error += temp * temp * 100;
            break;

        case Constraint::PointOnLine:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;

            t = ((P1_x - L1_P1_x) * dx + (P1_y - L1_P1_y) * dy) / (dx * dx + dy * dy);
            temp = _hypot((L1_P1_x + dx * t) - P1_x, (L1_P1_y + dy * t) - P1_y);
            error += temp * temp / 1000;
            break;

        case Constraint::P2LDistance:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;

            t = -(L1_P1_x * dx - P1_x * dx + L1_P1_y * dy - P1_y * dy) / (dx * dx + dy * dy);
            temp = _hypot(P1_x - (L1_P1_x + dx * t), P1_y - (L1_P1_y + dy * t)) - distance;
            error += temp * temp / 100;
            break;

        case Constraint::P2LDistanceVert:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;

            t = (P1_x - L1_P1_x) / dx;
            temp = fabs(P1_y - (L1_P1_y + dy * t)) - distance;
            error += temp * temp / 100;
            break;

        case Constraint::P2LDistanceHorz:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;

            t = (P1_y - L1_P1_y) / dy;
            temp = fabs(P1_x - (L1_P1_x + dx * t)) - distance;
            error += temp * temp / 100;
            break;

        case Constraint::LineLength:
            temp = _hypot(L1_P2_x - L1_P1_x, L1_P2_y - L1_P1_y) - distance;
            error += temp * temp * 100;
            break;

        case Constraint::EqualLegnth:
            temp = _hypot(L1_P2_x - L1_P1_x, L1_P2_y - L1_P1_y) -
                   _hypot(L2_P2_x - L2_P1_x, L2_P2_y - L2_P1_y);
            error += temp * temp * 100;
            break;

        case Constraint::Vertical:
            dx = L1_P2_x - L1_P1_x;
            error += dx * dx * 100;
            break;

        case Constraint::Horizontal:
            dy = L1_P2_y - L1_P1_y;
            error += dy * dy * 100;
            break;

        case Constraint::TangentToCircle:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;
            hyp1 = _hypot(dx, dy);
            // Calculate the expected tangent intersection points
            temp = (-dy * (C1_Center_x - dy / hyp1 * C1_rad) +
                    dx * (C1_Center_y + dx / hyp1 * C1_rad) +
                    (L1_P1_x * L1_P2_y - L1_P2_x * L1_P1_y)) / hyp1;
            temp2 = (-dy * (C1_Center_x + dy / hyp1 * C1_rad) +
                     dx * (C1_Center_y - dx / hyp1 * C1_rad) +
                     (L1_P1_x * L1_P2_y - L1_P2_x * L1_P1_y)) / hyp1;
            error += (temp < temp2) ? temp * temp : temp2 * temp2;
            break;

        case Constraint::TangentToArc:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;
            t = -(L1_P1_x * dx -
                  A1_Center_x * dx +
                  L1_P1_y * dy -
                  A1_Center_y * dy) / (dx * dx + dy * dy);
            Ex = A1_Center_x - (L1_P1_x + dx * t);
            Ey = A1_Center_y - (L1_P1_y + dy * t);
            temp = Ex * Ex + Ey * Ey -
                   (A1_Center_x - A1_Start_x) * (A1_Center_x - A1_Start_x) -
                   (A1_Center_y - A1_Start_y) * (A1_Center_y - A1_Start_y);
            error += temp * temp;
            break;

        case Constraint::ArcRules:
            dx = A1_End_x * A1_End_x;
            dy = A1_End_y * A1_End_y;
            rad1 = A1_Start_x * A1_Start_x;
            rad2 = A1_Start_y * A1_Start_y;
            temp = -2 * A1_Center_x * A1_End_x + dx -
                    2 * A1_Center_y * A1_End_y + dy +
                    2 * A1_Center_x * A1_Start_x - rad1 +
                    2 * A1_Center_y * A1_Start_y - rad2;
            error += temp * temp / (4. * dx + dy -
                                    2 * A1_End_x * A1_Start_x + rad1 -
                                    2 * A1_End_y * A1_Start_y + rad2);
            break;

        case Constraint::ArcRadius:
            rad1 = _hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
            rad2 = _hypot(A1_Center_x - A1_End_x, A1_Center_y - A1_End_y);
            temp= rad1 - radius;
            error += temp * temp;
            break;

        case Constraint::EqualRadiusArcs:
            temp = _hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y) -
                   _hypot(A2_Center_x - A2_Start_x, A2_Center_y - A2_Start_y);
            error += temp * temp;
            break;

        case Constraint::EqualRadiusCircles:
            temp = C1_rad - C2_rad;
            error += temp * temp;
            break;

        case Constraint::EqualRadiusCircArc:
            rad1 = _hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
            temp = rad1 - C1_rad;
            error += temp * temp;
            break;

        case Constraint::ConcentricArcs:
            temp = _hypot(A1_Center_x - A2_Center_x, A1_Center_y - A2_Center_y);
            error += temp * temp;
            break;

        case Constraint::ConcentricCircles:
            temp = _hypot(C1_Center_x - C2_Center_x, C1_Center_y - C2_Center_y);
            error += temp * temp;
            break;

        case Constraint::ConcentricCircArc:
            temp = _hypot(A1_Center_x - C1_Center_x, A1_Center_y - C1_Center_y);
            error += temp * temp;
            break;

        case Constraint::CircleRadius:
            error += (C1_rad - radius) * (C1_rad - radius);
            break;

        case Constraint::InternalAngle:
            temp = atan2(L2_P2_y - L2_P1_y, L2_P2_x - L2_P1_x) -
                   atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) - angleP;
            error += temp * temp;
            break;

        case Constraint::ExternalAngle:
            temp = atan2(L2_P2_y - L2_P1_y, L2_P2_x - L2_P1_x) -
                   atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) - (M_PI - angleP);
            error += temp * temp;
            break;

        case Constraint::LineInternalAngle:
            temp = sin(atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) - angleP);
            error += temp * temp;
            break;

        case Constraint::LineExternalAngle:
            temp = sin(atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) - (M_PI - angleP));
            error += temp * temp;
            break;

        case Constraint::Perpendicular:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;
            dx2 = L2_P2_x - L2_P1_x;
            dy2 = L2_P2_y - L2_P1_y;

            hyp1 = _hypot(dx, dy);
            hyp2 = _hypot(dx2, dy2);

            dx = dx / hyp1;
            dy = dy / hyp1;
            dx2 = dx2 / hyp2;
            dy2 = dy2 / hyp2;

            temp = dx * dx2 + dy * dy2;
            error += temp * temp;
            break;

        case Constraint::Parallel:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;
            dx2 = L2_P2_x - L2_P1_x;
            dy2 = L2_P2_y - L2_P1_y;

            hyp1 = _hypot(dx, dy);
            hyp2 = _hypot(dx2, dy2);

            dx = dx / hyp1;
            dy = dy / hyp1;
            dx2 = dx2 / hyp2;
            dy2 = dy2 / hyp2;

            temp = dy * dx2 - dx * dy2;
            error += temp * temp;
            break;

        case Constraint::Colinear:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;

            m = dy / dx;
            n = dx / dy;
            // Calculate the error between the expected intersection point
            // and the true point of the second lines two end points on the
            // first line
            if(m <= 1 && m > -1) {
                // Calculate the expected y point given the x coordinate of the point
                Ey = L1_P1_y + m * (L2_P1_x - L1_P1_x);
                error += (Ey - L2_P1_y) * (Ey - L2_P1_y);

                Ey = L1_P1_y + m * (L2_P2_x - L1_P1_x);
                error += (Ey - L2_P2_y) * (Ey - L2_P2_y);
            } else {
                // Calculate the expected x point given the y coordinate of the point
                Ex = L1_P1_x + n * (L2_P1_y - L1_P1_y);
                error += (Ex - L2_P1_x) * (Ex - L2_P1_x);

                Ex = L1_P1_x + n * (L2_P2_y - L1_P1_y);
                error += (Ex - L2_P2_x) * (Ex - L2_P2_x);
            }
            break;

        case Constraint::PointOnCircle:
            // see what the current radius to the point is
            rad1 = _hypot(C1_Center_x - P1_x, C1_Center_y - P1_y);
            // Compare this radius to the radius of the circle, return the error squared
            temp = rad1 - C1_rad;
            error += temp * temp;
            break;

        case Constraint::PointOnArc:
            // see what the current radius to the point is
            rad1 = _hypot(A1_Center_x - P1_x, A1_Center_y - P1_y);
            rad2 = _hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
            // Compare this radius to the radius of the circle, return the error squared
            temp = rad1 - rad2;
            error += temp * temp;
            break;

        case Constraint::PointOnLineMidpoint:
            Ex = (L1_P1_x + L1_P2_x) / 2 - P1_x;
            Ey = (L1_P1_y + L1_P2_y) / 2 - P1_y;
            error += Ex * Ex + Ey * Ey;
            break;

        case Constraint::PointOnArcMidpoint:
            rad1 = _hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
            temp = atan2(A1_Start_y - A1_Center_y, A1_Start_x - A1_Center_x);
            temp2 = atan2(A1_End_y - A1_Center_y, A1_End_x - A1_Center_x);
            Ex = A1_Center_x + rad1 * cos((temp2 + temp) / 2);
            Ey = A1_Center_y + rad1 * sin((temp2 + temp) / 2);
            temp = Ex - P1_x;
            temp2 = Ey - P1_y;
            error += temp * temp + temp2 * temp2;
            break;

        case Constraint::PointOnCircleQuad:
            Ex = C1_Center_x;
            Ey = C1_Center_y;
            switch(quadIndex) {
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
            temp = Ex - P1_x;
            temp2 = Ey - P1_y;
            error += temp * temp + temp2 * temp2;
            break;

        case Constraint::SymmetricPoints:
            dx = Sym_P2_x - Sym_P1_x;
            dy = Sym_P2_y - Sym_P1_y;
            t = -(dy * P1_x -
                  dx * P1_y -
                  dy * Sym_P1_x +
                  dx * Sym_P1_y) / (dx * dx + dy * dy);
            Ex = P1_x + dy * t * 2;
            Ey = P1_y - dx * t * 2;
            temp = Ex - P2_x;
            temp2 = Ey - P2_y;
            error += temp * temp + temp2 * temp2;
            break;

        case Constraint::SymmetricLines:
            dx = Sym_P2_x - Sym_P1_x;
            dy = Sym_P2_y - Sym_P1_y;
            hyp1 = dx * dx + dy * dy;
            m = -dy * Sym_P1_x + dx * Sym_P1_y;

            t = -(dy * L1_P1_x - dx * L1_P1_y + m) / hyp1;
            Ex = L1_P1_x + dy * t * 2;
            Ey = L1_P1_y - dx * t * 2;
            temp = Ex - L2_P1_x;
            temp2 = Ey - L2_P1_y;
            error += temp * temp + temp2 * temp2;
            t = -(dy * L1_P2_x - dx * L1_P2_y + m) / hyp1;
            Ex = L1_P2_x + dy * t * 2;
            Ey = L1_P2_y - dx * t * 2;
            temp = Ex - L2_P2_x;
            temp2 = Ey - L2_P2_y;
            error += temp * temp + temp2 * temp2;
            break;

        case Constraint::SymmetricCircles:
            dx = Sym_P2_x - Sym_P1_x;
            dy = Sym_P2_y - Sym_P1_y;
            t = -(dy * C1_Center_x -
                  dx * C1_Center_y -
                  dy * Sym_P1_x +
                  dx * Sym_P1_y) / (dx * dx + dy * dy);
            Ex = C1_Center_x + dy * t * 2;
            Ey = C1_Center_y - dx * t * 2;
            temp = Ex - C2_Center_x;
            temp2 = Ey - C2_Center_y;
            error += temp * temp + temp2 * temp2;
            temp = C1_rad - C2_rad;
            error += temp * temp;
            break;

        case Constraint::SymmetricArcs:
            dx = Sym_P2_x - Sym_P1_x;
            dy = Sym_P2_y - Sym_P1_y;
            hyp1 = dx * dx + dy * dy;
            m = - dy * Sym_P1_x + dx * Sym_P1_y;

            t = -(dy * A1_Start_x - dx * A1_Start_y + m) / hyp1;
            Ex = A1_Start_x + dy * t * 2;
            Ey = A1_Start_y - dx * t * 2;
            temp = Ex - A2_Start_x;
            temp2 = Ey - A2_Start_y;
            error += temp * temp + temp2 * temp2;
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
            error += temp * temp + temp2 * temp2;
            break;
        }
#ifdef DEBUG
    cout << "Type: " << cons[i].type << " Error: " << error << endl;
#endif
    }
    return error;
}
