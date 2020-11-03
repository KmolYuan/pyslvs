#ifndef SOLVER_H
#define SOLVER_H

#include <utility>
#include <vector>
#include <map>
#include "swappable_pair.hpp"

enum Label {
    P_LABEL,
    L_LABEL,
    I_LABEL,
    A_LABEL,
    S_LABEL,
};

enum Func {
    PXY,
    PPP,
    PLA,
    PLAP,
    PLLP,
    PLPP,
    PALP,
};

using Sym = std::pair<int, int>;

struct CCoord {
    double x, y;
};

auto cpxy(CCoord c1, double x, double y) -> CCoord;
auto cppp(CCoord c1, CCoord c2, CCoord c3) -> CCoord;
auto cplap(CCoord c1, double d0, double a0, CCoord c2, bool inverse) -> CCoord;
auto cpllp(CCoord c1, double d0, double d1, CCoord c2, bool inverse) -> CCoord;
auto cplpp(CCoord c1, double d0, CCoord c2, CCoord c3, bool inverse) -> CCoord;
auto cpalp(CCoord c1, double a0, double d0, CCoord c2, bool inverse) -> CCoord;

struct Expr {
    bool op;
    Func func;
    Sym v1, v2, c1, c2, c3, target;
};

using Stack = std::vector<Expr>;
using JointPos = std::map<Sym, CCoord>;
using LinkLen = std::map<SwappablePair, double>;
using Param = std::map<Sym, double>;

class ExprSolver {
    Stack stack;
    // Link length for P_LABEL.
    Param param;

public:
    JointPos joint_pos;

    ExprSolver() = default;
    ExprSolver(Stack stack, JointPos joint_pos, Param param);

    auto solve() -> bool;
};

#endif //SOLVER_H
