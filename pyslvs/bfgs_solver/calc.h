#ifndef CLAC_H
#define CLAC_H

/*
 *  Created on: July 24, 2019
 *  Author: KmolYuan
 */

#include "solve.h"

enum {
    // Geometric Constraint types
    POINT_ON_POINT,
    POINT_ON_LINE,
    INTERNAL_ANGLE,
    P2P_DISTANCE,
    LINE_INTERNAL_ANGLE,
};

double calc(Constraint *, size_t);

#endif  // CLAC_H
