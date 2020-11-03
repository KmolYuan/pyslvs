#include <iostream>
#include <cmath>
#include "solver.h"

int main() {
    auto stack = Stack{
        {false, PLA, Sym{L_LABEL, 0}, Sym{I_LABEL, 0},
            Sym{P_LABEL, 0}, Sym{}, Sym{}, Sym{P_LABEL, 1}},
    };
    auto join_pos = JointPos{
        {{P_LABEL, 0}, {0, 0}},
        {{P_LABEL, 1}, {0, 0}},
    };
    auto param = Param{
        {{I_LABEL, 0}, 90. / 180 * M_PI},
        {{L_LABEL, 0}, 10.},
    };
    auto solver = ExprSolver(stack, join_pos, param);
    solver.solve();
    for (auto[_, pos] : solver.joint_pos)
        std::cout << pos.x << ", " << pos.y << std::endl;
}
