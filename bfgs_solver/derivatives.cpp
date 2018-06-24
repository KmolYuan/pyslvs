/*
 *  Created on: Jun 8, 2009
 *  Author: jgeorge
 *  Contributor: KmolYuan
 */

#include <iostream>
#include "solve.h"

using namespace std;


void derivatives(double **x,
    double *gradF,
    const int xLength,
    Constraint *cons,
    const int consLength
) {
    int position;
    for(int i = 0; i < consLength; i++)
        switch(cons[i].type) {
        //////////////////////////////////////
        //Point on Point Constraint derivative
        //////////////////////////////////////
        case PointOnPoint:
            // Derivative with respect to p1x
            position = &P1_x - x[0];
            if(position >= 0 && position < xLength)
                gradF[position] += 2 * (P1_x - P2_x);

            // Derivative with respect to p1y
            position = &P1_y - x[0];
            if(position >= 0 && position < xLength)
                gradF[position] += 2 * (P1_y - P2_y);

            // Derivative with respect to p2x
            position = &P2_x - x[0];
            if(position >= 0 && position < xLength)
                gradF[position] += -2 * (P1_x - P2_x);

            // Derivative with respect to p2y
            position = &P2_y - x[0];
            if(position >=0 && position<xLength)
                gradF[position] += -2 * (P1_y - P2_y);
            break;

        //////////////////////////////////////
        //Point to Point Distance Constraint derivative
        //////////////////////////////////////
        case P2PDistance:
            // Derivative with respect to p1x
            position = &P1_x - x[0];
            if(position >= 0 && position < xLength)
                gradF[position] += 2 * (P1_x - P2_x);

            // Derivative with respect to p1y
            position = &P1_y - x[0];
            if(position >= 0 && position < xLength)
                gradF[position] += 2 * (P1_y - P2_y);

            // Derivative with respect to p2x
            position = &P2_x - x[0];
            if(position >= 0 && position < xLength)
                gradF[position] += -2 * (P1_x - P2_x);

            // Derivative with respect to p2y
            position = &P2_y - x[0];
            if(position >=0 && position < xLength)
                gradF[position] += -2 * (P1_y - P2_y);

            // Derivative with respect to distance
            position = &distance - x[0];
            if(position >= 0 && position < xLength)
                gradF[position] += -2 * distance;
            break;

        //////////////////////////////////////
        //Point to Point Distance Vert Constraint derivative
        //////////////////////////////////////
        case P2PDistanceVert:
            // Derivative with respect to p1y
            position = &P1_y - x[0];
            if(position >=0 && position < xLength)
                gradF[position] += 2 * (P1_y - P2_y);

            // Derivative with respect to p2y
            position = &P2_y - x[0];
            if(position >= 0 && position < xLength)
                gradF[position] += -2 * (P1_y - P2_y);

            // Derivative with respect to distance
            position = &distance - x[0];
            if(position >= 0 && position < xLength)
                gradF[position] += -2 * distance;
            break;

        //////////////////////////////////////
        //Point to Point Horz Distance Constraint derivative
        //////////////////////////////////////
        case P2PDistanceHorz:
            // Derivative with respect to p1x
            position = &P1_x - x[0];
            if(position >= 0 && position < xLength)
                gradF[position] += 2 * (P1_x - P2_x);

            // Derivative with respect to p2x
            position = &P2_x - x[0];
            if(position >= 0 && position < xLength)
                gradF[position] += -2 * (P1_x - P2_x);

            // Derivative with respect to distance
            position = &distance - x[0];
            if(position >= 0 && position < xLength)
                gradF[position] += -2 * distance;
            break;

        //////////////////////////////////////
        //Point on line Constraint derivatives
        //////////////////////////////////////
        case PointOnLine:
            // Derivative with respect to p1x
            position = &P1_x - x[0];
            if(position >= 0 && position < xLength)
                gradF[position] += 2 * (P1_x - P2_x);

            // Derivative with respect to p2x
            position = &P2_x - x[0];
            if(position >= 0 && position < xLength)
                gradF[position] += -2 * (P1_x - P2_x);

            // Derivative with respect to distance
            position = &distance - x[0];
            if(position >= 0 && position < xLength)
                gradF[position] += -2 * distance;
            break;
        }
}
