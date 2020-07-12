/*
 *  Created on: July 24, 2019
 *  Author: KmolYuan
 */

#include "calc.h"
#include <cmath>

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
#define ANGLE (*con->parameter)
#define DISTANCE fabs(ANGLE)

namespace {
auto point_on_point(Constraint *con) -> double {
    // Hopefully avoid this constraint,
    // make coincident points use the same parameters
    auto dx = P2_x - P1_x;
    auto dy = P2_y - P1_y;
    return dx * dx + dy * dy;
}

auto p2p_distance(Constraint *con) -> double {
    auto tmp = hypot(P2_x - P1_x, P2_y - P1_y) - DISTANCE;
    return tmp * tmp * 100;
}

auto point_on_line(Constraint *con) -> double {
    auto dx = L1_P2_x - L1_P1_x;
    auto dy = L1_P2_y - L1_P1_y;
    auto t =
        ((P1_x - L1_P1_x) * dx + (P1_y - L1_P1_y) * dy) / (dx * dx + dy * dy);
    dx = (L1_P1_x + dx * t) - P1_x;
    dy = (L1_P1_y + dy * t) - P1_y;
    return (dx * dx + dy * dy) * 1e-3;
}

auto internal_angle(Constraint *con) -> double {
    auto tmp = atan2(L2_P2_y - L2_P1_y, L2_P2_x - L2_P1_x)
               - atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) - ANGLE;
    return tmp * tmp;
}

auto line_internal_angle(Constraint *con) -> double {
    auto tmp = sin(atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) - ANGLE);
    return tmp * tmp;
}
}  // namespace

auto calc(Constraint *cons, size_t cons_len) -> double {
    auto error = 0.;
    for (auto i = 0u; i < cons_len; i++) {
        auto con = cons + i;
        switch (con->type) {
            case POINT_ON_POINT:
                error += point_on_point(con);
                break;
            case P2P_DISTANCE:
                error += p2p_distance(con);
                break;
            case POINT_ON_LINE:
                error += point_on_line(con);
                break;
            case INTERNAL_ANGLE:
                error += internal_angle(con);
                break;
            case LINE_INTERNAL_ANGLE:
                error += line_internal_angle(con);
                break;
        }
    }
    return error;
}
