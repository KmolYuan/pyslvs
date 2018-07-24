/*
 *  Created on: July 24, 2018
 *  Author: KmolYuan
 */

#include "position.h"
#include "constrain_func.h"

double calc(Constraint *cons, const int consLength) {
    double error = 0;
    double temp, temp2, dx, dy, m, n, Ex, Ey, rad1, rad2, t, dx2, dy2, hyp1, hyp2;
    for(int i = 0; i < consLength; i++) {
        switch(cons[i].type) {
        case PointOnPoint:
            //Hopefully avoid this constraint, make coincident points use the same parameters
            dx = P1_x - P2_x;
            dy = P1_y - P2_y;
            error += dx * dx + dy * dy;
            break;

        case P2PDistance:
            temp = _hypot(P2_x - P1_x, P2_y - P1_y) - distance;
            error += temp * temp * 100;
            break;

        case P2PDistanceVert:
            temp = _hypot(0, P2_y - P1_y) - distance;
            error += temp * temp * 100;
            break;

        case P2PDistanceHorz:
            temp = _hypot(P2_x - P1_x, 0) - distance;
            error += temp * temp * 100;
            break;

        case PointOnLine:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;

            t = ((P1_x - L1_P1_x) * dx + (P1_y - L1_P1_y) * dy) / (dx * dx + dy * dy);
            temp = _hypot((L1_P1_x + dx * t) - P1_x, (L1_P1_y + dy * t) - P1_y);
            error += temp * temp / 100;
            break;

        case P2LDistance:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;

            t = -(L1_P1_x * dx - P1_x * dx + L1_P1_y * dy - P1_y * dy) / (dx * dx + dy * dy);
            temp = _hypot(P1_x - (L1_P1_x + dx * t), P1_y - (L1_P1_y + dy * t)) - distance;
            error += temp * temp / 100;
            break;

        case P2LDistanceVert:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;

            t = (P1_x - L1_P1_x) / dx;
            temp = fabs(P1_y - (L1_P1_y + dy * t)) - distance;
            error += temp * temp / 100;
            break;

        case P2LDistanceHorz:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;

            t = (P1_y - L1_P1_y) / dy;
            temp = fabs(P1_x - (L1_P1_x + dx * t)) - distance;
            error += temp * temp / 100;
            break;

        case Vertical:
            dx = L1_P2_x - L1_P1_x;
            error += dx * dx * 100;
            break;

        case Horizontal:
            dy = L1_P2_y - L1_P1_y;
            error += dy * dy * 100;
            break;

        case TangentToCircle:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;
            hyp1 = _hypot(dx, dy);
            //Calculate the expected tangent intersection points
            temp = (-dy * (C1_Center_x - dy / hyp1 * C1_rad) +
                    dx * (C1_Center_y + dx / hyp1 * C1_rad) +
                    (L1_P1_x * L1_P2_y - L1_P2_x * L1_P1_y)) / hyp1;
            temp2 = (-dy * (C1_Center_x + dy / hyp1 * C1_rad) +
                     dx * (C1_Center_y - dx / hyp1 * C1_rad) +
                     (L1_P1_x * L1_P2_y - L1_P2_x * L1_P1_y)) / hyp1;
            error += (temp < temp2) ? temp * temp : temp2 * temp2;
            break;

        case TangentToArc:
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

        case ArcRules:
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

        case LineLength:
            temp = _hypot(L1_P2_x - L1_P1_x, L1_P2_y - L1_P1_y) - length;
            error += temp * temp * 100;
            break;

        case EqualLegnth:
            temp = _hypot(L1_P2_x - L1_P1_x, L1_P2_y - L1_P1_y) -
                   _hypot(L2_P2_x - L2_P1_x, L2_P2_y - L2_P1_y);
            error += temp * temp * 100;
            break;

        case ArcRadius:
            rad1 = _hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
            rad2 = _hypot(A1_Center_x - A1_End_x, A1_Center_y - A1_End_y);
            temp= rad1 - radius;
            error += temp * temp;
            break;

        case EqualRadiusArcs:
            temp = _hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y) -
                   _hypot(A2_Center_x - A2_Start_x, A2_Center_y - A2_Start_y);
            error += temp * temp;
            break;

        case EqualRadiusCircles:
            temp = C1_rad - C2_rad;
            error += temp * temp;
            break;

        case EqualRadiusCircArc:
            rad1 = _hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
            temp = rad1 - C1_rad;
            error += temp * temp;
            break;

        case ConcentricArcs:
            temp = _hypot(A1_Center_x - A2_Center_x, A1_Center_y - A2_Center_y);
            error += temp * temp;
            break;

        case ConcentricCircles:
            temp = _hypot(C1_Center_x - C2_Center_x, C1_Center_y - C2_Center_y);
            error += temp * temp;
            break;

        case ConcentricCircArc:
            temp = _hypot(A1_Center_x - C1_Center_x, A1_Center_y - C1_Center_y);
            error += temp * temp;
            break;

        case CircleRadius:
            error += (C1_rad - radius) * (C1_rad - radius);
            break;

        case InternalAngle:
            temp = atan2(L2_P2_y - L2_P1_y, L2_P2_x - L2_P1_x) -
                   atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) - angleP;
            error += temp * temp;
            break;

        case ExternalAngle:
            temp = atan2(L2_P2_y - L2_P1_y, L2_P2_x - L2_P1_x) -
                   atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) - (M_PI - angleP);
            error += temp * temp;
            break;

        case LineInternalAngle:
            temp = sin(atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) - angleP);
            error += temp * temp;
            break;

        case LineExternalAngle:
            temp = sin(atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) - (M_PI - angleP));
            error += temp * temp;
            break;

        case Perpendicular:
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

        case Parallel:
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

        case Colinear:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;

            m = dy / dx;
            n = dx / dy;
            // Calculate the error between the expected intersection point
            // and the true point of the second lines two end points on the
            // first line
            if(m <= 1 && m > -1) {
                //Calculate the expected y point given the x coordinate of the point
                Ey = L1_P1_y + m * (L2_P1_x - L1_P1_x);
                error += (Ey - L2_P1_y) * (Ey - L2_P1_y);

                Ey = L1_P1_y + m * (L2_P2_x - L1_P1_x);
                error += (Ey - L2_P2_y) * (Ey - L2_P2_y);
            } else {
                //Calculate the expected x point given the y coordinate of the point
                Ex = L1_P1_x + n * (L2_P1_y - L1_P1_y);
                error += (Ex - L2_P1_x) * (Ex - L2_P1_x);

                Ex = L1_P1_x + n * (L2_P2_y - L1_P1_y);
                error += (Ex - L2_P2_x) * (Ex - L2_P2_x);
            }
            break;

        case PointOnCircle:
            //see what the current radius to the point is
            rad1 = _hypot(C1_Center_x - P1_x, C1_Center_y - P1_y);
            //Compare this radius to the radius of the circle, return the error squared
            temp = rad1 - C1_rad;
            error += temp * temp;
            break;

        case PointOnArc:
            //see what the current radius to the point is
            rad1 = _hypot(A1_Center_x - P1_x, A1_Center_y - P1_y);
            rad2 = _hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
            //Compare this radius to the radius of the circle, return the error squared
            temp = rad1 - rad2;
            error += temp * temp;
            break;

        case PointOnLineMidpoint:
            Ex = (L1_P1_x + L1_P2_x) / 2 - P1_x;
            Ey = (L1_P1_y + L1_P2_y) / 2 - P1_y;
            error += Ex * Ex + Ey * Ey;
            break;

        case PointOnArcMidpoint:
            rad1 = _hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
            temp = atan2(A1_Start_y - A1_Center_y, A1_Start_x - A1_Center_x);
            temp2 = atan2(A1_End_y - A1_Center_y, A1_End_x - A1_Center_x);
            Ex = A1_Center_x + rad1 * cos((temp2 + temp) / 2);
            Ey = A1_Center_y + rad1 * sin((temp2 + temp) / 2);
            temp = Ex - P1_x;
            temp2 = Ey - P1_y;
            error += temp * temp + temp2 * temp2;
            break;

        case PointOnCircleQuad:
            Ex = C1_Center_x;
            Ey = C1_Center_y;
            switch((int)quadIndex) {
            case 1:
                Ey += C1_rad;
                break;
            case 2:
                Ex -= C1_rad;
                break;
            case 3:
                Ey -= C1_rad;
                break;
            default:
                Ex += C1_rad;
                break;
            }
            temp = Ex - P1_x;
            temp2 = Ey - P1_y;
            error += temp * temp + temp2 * temp2;
            break;

        case SymmetricPoints:
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

        case SymmetricLines:
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

        case SymmetricCircles:
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

        case SymmetricArcs:
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
    }
    return error;
}
