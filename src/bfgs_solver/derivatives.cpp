/*
 *  Created on: Jun 8, 2009
 *  Author: jgeorge
 *  Contributor: KmolYuan
 */

#include <iostream>
#include "solve.h"

using namespace std;


void derivatives(
    double **x,
    double *gradF,
    const int xLength,
    Constraint *cons,
    const int consLength
) {
    int position;
    for (int i = 0; i < consLength; i++)
        switch((int)cons[i].type) {
        //////////////////////////////////////
        // Point on Point GeoConstraint derivative
        //////////////////////////////////////
        case PointOnPoint:
            // Derivative with respect to p1x
            position = cons[i].point1->x - x[0];
            if (position >= 0 && position < xLength)
                gradF[position] += 2 * (*cons[i].point1->x - *cons[i].point2->x);

            // Derivative with respect to p1y
            position = cons[i].point1->y - x[0];
            if (position >= 0 && position < xLength)
                gradF[position] += 2 * (*cons[i].point1->y - *cons[i].point2->y);

            // Derivative with respect to p2x
            position = cons[i].point2->x - x[0];
            if (position >= 0 && position < xLength)
                gradF[position] += -2 * (*cons[i].point1->x - *cons[i].point2->x);

            // Derivative with respect to p2y
            position = cons[i].point2->y - x[0];
            if (position >=0 && position<xLength)
                gradF[position] += -2 * (*cons[i].point1->y - *cons[i].point2->y);
            break;

        //////////////////////////////////////
        // Point to Point Distance GeoConstraint derivative
        //////////////////////////////////////
        case P2PDistance:
            // Derivative with respect to p1x
            position = cons[i].point1->x - x[0];
            if (position >= 0 && position < xLength)
                gradF[position] += 2 * (*cons[i].point1->x - *cons[i].point2->x);

            // Derivative with respect to p1y
            position = cons[i].point1->y - x[0];
            if (position >= 0 && position < xLength)
                gradF[position] += 2 * (*cons[i].point1->y - *cons[i].point2->y);

            // Derivative with respect to p2x
            position = cons[i].point2->x - x[0];
            if (position >= 0 && position < xLength)
                gradF[position] += -2 * (*cons[i].point1->x - *cons[i].point2->x);

            // Derivative with respect to p2y
            position = cons[i].point2->y - x[0];
            if (position >=0 && position < xLength)
                gradF[position] += -2 * (*cons[i].point1->y - *cons[i].point2->y);

            // Derivative with respect to distance
            position = cons[i].parameter - x[0];
            if (position >= 0 && position < xLength)
                gradF[position] += -2 * *cons[i].parameter;
            break;

        //////////////////////////////////////
        //Point to Point Distance Vert GeoConstraint derivative
        //////////////////////////////////////
        case P2PDistanceVert:
            // Derivative with respect to p1y
            position = cons[i].point1->y - x[0];
            if (position >= 0 && position < xLength)
                gradF[position] += 2 * (*cons[i].point1->y - *cons[i].point2->y);

            // Derivative with respect to p2y
            position = cons[i].point2->y - x[0];
            if (position >= 0 && position < xLength)
                gradF[position] += -2 * (*cons[i].point1->y - *cons[i].point2->y);

            // Derivative with respect to distance
            position = cons[i].parameter - x[0];
            if (position >= 0 && position < xLength)
                gradF[position] += -2 * *cons[i].parameter;
            break;

        //////////////////////////////////////
        // Point to Point Horz Distance GeoConstraint derivative
        //////////////////////////////////////
        case P2PDistanceHorz:
            // Derivative with respect to p1x
            position = cons[i].point1->x - x[0];
            if (position >= 0 && position < xLength)
                gradF[position] += 2 * (*cons[i].point1->x - *cons[i].point2->x);

            // Derivative with respect to p2x
            position = cons[i].point2->x - x[0];
            if (position >= 0 && position < xLength)
                gradF[position] += -2 * (*cons[i].point1->x - *cons[i].point2->x);

            // Derivative with respect to distance
            position = cons[i].parameter - x[0];
            if (position >= 0 && position < xLength)
                gradF[position] += -2 * *cons[i].parameter;
            break;

        //////////////////////////////////////
        // Point on line GeoConstraint derivatives
        //////////////////////////////////////
        case PointOnLine:
            // Derivative with respect to p1x
            position = cons[i].point1->x - x[0];
            if (position >= 0 && position < xLength)
                gradF[position] += 2 * (*cons[i].point1->x - *cons[i].point2->x);

            // Derivative with respect to p2x
            position = cons[i].point2->x - x[0];
            if (position >= 0 && position < xLength)
                gradF[position] += -2 * (*cons[i].point1->x - *cons[i].point2->x);

            // Derivative with respect to distance
            position = cons[i].parameter - x[0];
            if (position >= 0 && position < xLength)
                gradF[position] += -2 * *cons[i].parameter;
            break;
        }
}
