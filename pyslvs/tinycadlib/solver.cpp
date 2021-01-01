//
// Created by KmolYuan on 2021/10/12.
//

#include "solver.h"
#include <cmath>
#include <utility>

namespace {
    auto distance(double x1, double y1, double x2, double y2) -> double {
        return hypot(x1 - x2, y1 - y2);
    }

    auto slope_angle(double x1, double y1, double x2, double y2) -> double {
        return atan2(y1 - y2, x1 - x2);
    }
}

auto cpxy(CCoord c1, double x, double y) -> CCoord {
    return {c1.x + x, c1.y + y};
}

auto cppp(CCoord c1, CCoord c2, CCoord c3) -> CCoord {
    auto length = distance(c1.x, c1.y, c2.x, c2.y);
    auto alpha = slope_angle(c2.x, c2.y, c1.x, c1.y);
    return {c3.x + length * cos(alpha), c3.y + length * sin(alpha)};
}

auto cplap(CCoord c1, double d0, double a0, CCoord c2, bool inverse) -> CCoord {
    auto a1 = atan2(c2.y - c1.y, c2.x - c1.x);
    a1 += inverse ? -a0 : a0;
    return {c1.x + d0 * cos(a1), c1.y + d0 * sin(a1)};
}

auto cpllp(CCoord c1, double d0, double d1, CCoord c2, bool inverse) -> CCoord {
    auto dx = c2.x - c1.x;
    auto dy = c2.y - c1.y;
    auto d = distance(c1.x, c1.y, c2.x, c2.y);
    if (d > d0 + d1)
        return {NAN, NAN};
    if (d < fabs(d0 - d1))
        return {NAN, NAN};
    if (d == 0 && d0 == d1)
        return {NAN, NAN};
    auto a = (d0 * d0 - d1 * d1 + d * d) / (2 * d);
    auto h = sqrt(d0 * d0 - a * a);
    auto xm = c1.x + a * dx / d;
    auto ym = c1.y + a * dy / d;
    if (inverse)
        return {xm + h * dy / d, ym - h * dx / d};
    else
        return {xm - h * dy / d, ym + h * dx / d};
}

auto cplpp(CCoord c1, double d0, CCoord c2, CCoord c3, bool inverse) -> CCoord {
    auto line_mag = distance(c2.x, c2.y, c3.x, c3.y);
    auto dx = c3.x - c2.x;
    auto dy = c3.y - c2.y;
    auto u = ((c1.x - c2.x) * dx + (c1.y - c2.y) * dy) / (line_mag * line_mag);
    auto inter = CCoord{c2.x + u * dx, c2.y + u * dy};
    auto d = distance(c1.x, c1.y, inter.x, inter.y);
    if (d > d0)
        return {NAN, NAN};
    else if (d == d0)
        return inter;
    d = sqrt(d0 * d0 - d * d) / line_mag;
    dx *= d;
    dy *= d;
    if (inverse)
        return {inter.x - dx, inter.y - dy};
    else
        return {inter.x + dx, inter.y + dy};
}

auto cpalp(CCoord c1, double a0, double d0, CCoord c2, bool inverse) -> CCoord {
    a0 += slope_angle(c2.x, c2.y, c1.x, c1.y);
    auto tan_a = tan(a0);
    auto tan2_a = tan_a * tan_a;
    auto tan2_a1 = tan2_a + 1;
    auto c1l = c1.x - c1.y / tan_a;
    auto c1c2x = c1.x - c2.x;
    auto c1c2y = c1.y - c2.y;
    auto sq = sqrt(d0 * d0 * tan2_a1 - c1c2x * c1c2x * tan2_a - c1c2y * c1c2y + 2 * tan_a * c1c2y * c1c2x);
    auto cx = c1l - (c1l - c2.y * tan_a - c2.x + (inverse ? -sq : sq)) / tan2_a1;
    return {cx, tan_a * (cx - c1.x) + c1.y};
}


/// Constructure of solver core.
/// Check the link length override options.
ExprSolver::ExprSolver(Stack stack, JointPos joint_pos, Param param) :
    stack(std::move(stack)),
    param(std::move(param)),
    joint_pos(std::move(joint_pos)) {}

auto ExprSolver::solve() -> bool {
    for (auto expr : stack) {
        auto &c = joint_pos[expr.target];
        switch (expr.func) {
            case PXY:
                c = cpxy(joint_pos[expr.c1],
                         param[expr.v1],
                         param[expr.v2]);
                break;
            case PPP:
                c = cppp(joint_pos[expr.c1],
                         joint_pos[expr.c2],
                         joint_pos[expr.c3]);
                break;
            case PLA:
                c = cplap(joint_pos[expr.c1],
                          param[expr.v1],
                          param[expr.v2],
                          joint_pos[expr.c1],
                          expr.op);
                break;
            case PLAP:
                c = cplap(joint_pos[expr.c1],
                          param[expr.v1],
                          param[expr.v2],
                          joint_pos[expr.c2],
                          expr.op);
                break;
            case PLLP:
                c = cpllp(joint_pos[expr.c1],
                          param[expr.v1],
                          param[expr.v2],
                          joint_pos[expr.c2],
                          expr.op);
                break;
            case PLPP:
                c = cplpp(joint_pos[expr.c1],
                          param[expr.v1],
                          joint_pos[expr.c2],
                          joint_pos[expr.c3],
                          expr.op);
                break;
            case PALP:
                c = cpalp(joint_pos[expr.c1],
                          param[expr.v1],
                          param[expr.v2],
                          joint_pos[expr.c2],
                          expr.op);
        }
        if (std::isnan(c.x) || std::isnan(c.y))
            return false;
    }
    return true;
}
