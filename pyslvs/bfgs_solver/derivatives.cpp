/*
 *  Created on: Jun 8, 2009
 *  Author: jgeorge
 *  Contributor: KmolYuan
 */

#include "derivatives.h"

void derivatives(double **x, double *grad, size_t x_len, Constraint *cons,
                 size_t cons_len) {
    for (auto i = 0u; i < cons_len; i++) {
        auto con = cons + i;
        switch (con->type) {
            //////////////////////////////////////
            /// Point on Point GeoConstraint derivative
            //////////////////////////////////////
            case POINT_ON_POINT: {
                // Derivative with respect to p1x
                auto position = size_t(con->point1->x - x[0]);
                if (position >= 0 && position < x_len)
                    grad[position] += 2 * (*con->point1->x - *con->point2->x);

                // Derivative with respect to p1y
                position = con->point1->y - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += 2 * (*con->point1->y - *con->point2->y);

                // Derivative with respect to p2x
                position = con->point2->x - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] -= 2 * (*con->point1->x - *con->point2->x);

                // Derivative with respect to p2y
                position = con->point2->y - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] -= 2 * (*con->point1->y - *con->point2->y);
                break;
            }
            //////////////////////////////////////
            /// Point to Point Distance GeoConstraint derivative
            //////////////////////////////////////
            case P2P_DISTANCE: {
                // Derivative with respect to p1x
                auto position = size_t(con->point1->x - x[0]);
                if (position >= 0 && position < x_len)
                    grad[position] += 2 * (*con->point1->x - *con->point2->x);

                // Derivative with respect to p1y
                position = con->point1->y - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] += 2 * (*con->point1->y - *con->point2->y);

                // Derivative with respect to p2x
                position = con->point2->x - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] -= 2 * (*con->point1->x - *con->point2->x);

                // Derivative with respect to p2y
                position = con->point2->y - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] -= 2 * (*con->point1->y - *con->point2->y);

                // Derivative with respect to DISTANCE
                position = con->parameter - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] -= 2 * *con->parameter;
                break;
            }
            //////////////////////////////////////
            /// Point on line GeoConstraint derivatives
            //////////////////////////////////////
            case POINT_ON_LINE: {
                // Derivative with respect to p1x
                auto position = size_t(con->point1->x - x[0]);
                if (position >= 0 && position < x_len)
                    grad[position] += 2 * (*con->point1->x - *con->point2->x);

                // Derivative with respect to p2x
                position = con->point2->x - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] -= 2 * (*con->point1->x - *con->point2->x);

                // Derivative with respect to DISTANCE
                position = con->parameter - x[0];
                if (position >= 0 && position < x_len)
                    grad[position] -= 2 * *con->parameter;
                break;
            }
        }
    }
}
