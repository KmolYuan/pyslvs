/*
 *  Created on: Jun 8, 2009
 *  Author: jgeorge
 *  Contributor: KmolYuan
 */

#include "calc.h"

using namespace std;

void derivatives(double **x, double *grad, size_t x_len, Constraint *cons,
                 size_t cons_len) {
    for (size_t i = 0; i < cons_len; i++) {
        size_t position;
        Constraint *con = cons + i;
        switch (con->type) {
            //////////////////////////////////////
            /// Point on Point GeoConstraint derivative
            //////////////////////////////////////
            case PointOnPoint:
                // Derivative with respect to p1x
                position = con->point1->x - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += 2 * (*con->point1->x - *con->point2->x);

                // Derivative with respect to p1y
                position = con->point1->y - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += 2 * (*con->point1->y - *con->point2->y);

                // Derivative with respect to p2x
                position = con->point2->x - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += -2 * (*con->point1->x - *con->point2->x);

                // Derivative with respect to p2y
                position = con->point2->y - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += -2 * (*con->point1->y - *con->point2->y);
                break;
            //////////////////////////////////////
            /// Point to Point Distance GeoConstraint derivative
            //////////////////////////////////////
            case P2PDistance:
                // Derivative with respect to p1x
                position = con->point1->x - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += 2 * (*con->point1->x - *con->point2->x);

                // Derivative with respect to p1y
                position = con->point1->y - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += 2 * (*con->point1->y - *con->point2->y);

                // Derivative with respect to p2x
                position = con->point2->x - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += -2 * (*con->point1->x - *con->point2->x);

                // Derivative with respect to p2y
                position = con->point2->y - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += -2 * (*con->point1->y - *con->point2->y);

                // Derivative with respect to DISTANCE
                position = con->parameter - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += -2 * *con->parameter;
                break;
            //////////////////////////////////////
            /// Point to Point Distance Vert GeoConstraint derivative
            //////////////////////////////////////
            case P2PDistanceVert:
                // Derivative with respect to p1y
                position = con->point1->y - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += 2 * (*con->point1->y - *con->point2->y);

                // Derivative with respect to p2y
                position = con->point2->y - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += -2 * (*con->point1->y - *con->point2->y);

                // Derivative with respect to DISTANCE
                position = con->parameter - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += -2 * *con->parameter;
                break;
            //////////////////////////////////////
            /// Point to Point Horz Distance GeoConstraint derivative
            //////////////////////////////////////
            case P2PDistanceHorz:
                // Derivative with respect to p1x
                position = con->point1->x - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += 2 * (*con->point1->x - *con->point2->x);

                // Derivative with respect to p2x
                position = con->point2->x - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += -2 * (*con->point1->x - *con->point2->x);

                // Derivative with respect to DISTANCE
                position = con->parameter - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += -2 * *con->parameter;
                break;
            //////////////////////////////////////
            /// Point on line GeoConstraint derivatives
            //////////////////////////////////////
            case PointOnLine:
                // Derivative with respect to p1x
                position = con->point1->x - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += 2 * (*con->point1->x - *con->point2->x);

                // Derivative with respect to p2x
                position = con->point2->x - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += -2 * (*con->point1->x - *con->point2->x);

                // Derivative with respect to DISTANCE
                position = con->parameter - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += -2 * *con->parameter;
                break;
        }
    }
}
