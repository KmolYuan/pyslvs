/*
 * Create on: Jun 6, 2019
 * Author: KmolYuan
 * These functions make writing constraint easier
 */

#include "calc.h"

namespace {
auto p2p_distance_constraint(unsigned type, Point *point1, Point *point2,
                             double *value) -> Constraint {
    auto con = Constraint{};
    con.type = type;
    con.point1 = point1;
    con.point2 = point2;
    con.parameter = value;
    return con;
}

auto point_line_constraint(unsigned type, Point *point1, Line *line1)
    -> Constraint {
    auto con = Constraint{};
    con.type = type;
    con.point1 = point1;
    con.line1 = line1;
    return con;
}

auto line_constraint(unsigned type, Line *line1, double *value) -> Constraint {
    auto con = Constraint{};
    con.type = type;
    con.line1 = line1;
    con.parameter = value;
    return con;
}

auto angle_constraint(unsigned type, Line *line1, Line *line2, double *value)
    -> Constraint {
    auto con = Constraint{};
    con.type = type;
    con.line1 = line1;
    con.line2 = line2;
    con.parameter = value;
    return con;
}
}  // namespace

auto point_on_point(Point *point1, Point *point2) -> Constraint {
    auto con = Constraint{};
    con.type = POINT_ON_POINT;
    con.point1 = point1;
    con.point2 = point2;
    return con;
}

auto p2p_distance(Point *point1, Point *point2, double *value) -> Constraint {
    return p2p_distance_constraint(P2P_DISTANCE, point1, point2, value);
}

[[maybe_unused]] auto point_on_line(Point *point1, Line *line1) -> Constraint {
    return point_line_constraint(POINT_ON_LINE, point1, line1);
}

auto internal_angle(Line *line1, Line *line2, double *value) -> Constraint {
    return angle_constraint(INTERNAL_ANGLE, line1, line2, value);
}

[[maybe_unused]] auto line_internal_angle(Line *line1, double *value)
    -> Constraint {
    return line_constraint(LINE_INTERNAL_ANGLE, line1, value);
}
