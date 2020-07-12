/*
 * Create on: Jun 6, 2019
 * Author: KmolYuan
 * These functions make writing constraint easier
 */

#include "calc.h"

namespace {
auto P2PDistanceConstraint(unsigned type, Point *point1, Point *point2,
                           double *value) -> Constraint {
    auto con = Constraint{};
    con.type = type;
    con.point1 = point1;
    con.point2 = point2;
    con.parameter = value;
    return con;
}

auto PointLineConstraint(unsigned type, Point *point1, Line *line1)
    -> Constraint {
    auto con = Constraint{};
    con.type = type;
    con.point1 = point1;
    con.line1 = line1;
    return con;
}

auto LineConstraint(unsigned type, Line *line1, double *value) -> Constraint {
    auto con = Constraint{};
    con.type = type;
    con.line1 = line1;
    con.parameter = value;
    return con;
}

auto AngleConstraint(unsigned type, Line *line1, Line *line2, double *value)
    -> Constraint {
    auto con = Constraint{};
    con.type = type;
    con.line1 = line1;
    con.line2 = line2;
    con.parameter = value;
    return con;
}
}  // namespace

auto PointOnPointConstraint(Point *point1, Point *point2) -> Constraint {
    auto con = Constraint{};
    con.type = POINT_ON_POINT;
    con.point1 = point1;
    con.point2 = point2;
    return con;
}

auto P2PDistanceConstraint(Point *point1, Point *point2, double *value)
    -> Constraint {
    return P2PDistanceConstraint(P2P_DISTANCE, point1, point2, value);
}

[[maybe_unused]] auto PointOnLineConstraint(Point *point1, Line *line1)
    -> Constraint {
    return PointLineConstraint(POINT_ON_LINE, point1, line1);
}

auto InternalAngleConstraint(Line *line1, Line *line2, double *value)
    -> Constraint {
    return AngleConstraint(INTERNAL_ANGLE, line1, line2, value);
}

[[maybe_unused]] auto LineInternalAngleConstraint(Line *line1, double *value)
    -> Constraint {
    return LineConstraint(LINE_INTERNAL_ANGLE, line1, value);
}
