/* ============================================================================
 * Name        : SolverPointers.cpp
 * Author      : Jonathan George
 * Contributor : KmolYuan
 *============================================================================*/

#include <cmath>
#include <iostream>
#include "derivatives.h"

using namespace std;

// Show coordinates
auto printpoints(Point *points, size_t N) {
    for (auto i = 0u; i < N; i++)
        cout << "Point " << i << ": (" << *points[i].x << ", " << *points[i].y
             << ")" << endl;
}

auto main() -> int {
    // Input a parameter list
    double parameters[] = {0, 0, 5, 0, 6, 5, 6, 5};
    const size_t param_count = sizeof(parameters) / sizeof(*parameters);
    // Make a list of pointers of parameters
    auto pparameters = new double *[param_count];
    for (int i = 0; i < param_count; i++)
        pparameters[i] = &parameters[i];
    // Input a constant parameter list
    double constants[] = {30, 10, 24, M_PI * 0.5};

    // Create geometric objects and constraints with pointers
    Point points[] = {
        {pparameters[0], pparameters[1]}, {pparameters[2], pparameters[3]},
        {pparameters[4], pparameters[5]}, {pparameters[6], pparameters[7]},
        {&constants[0], &constants[1]},
    };
    Line lines[] = {
        {&points[0], &points[1]},
        {&points[1], &points[2]},
    };
    Constraint cons[] = {
        internal_angle(&lines[0], &lines[1], &constants[3]),
        point_on_point(&points[2], &points[4]),
        point_on_point(&points[3], &points[4]),
        p2p_distance(&points[1], &points[2], &constants[2]),
    };
    const auto point_count = sizeof(points) / sizeof(*points);
    const auto cons_count = sizeof(cons) / sizeof(*cons);
    printpoints(points, point_count);

    // Solve
    if (solve(pparameters, param_count, cons, cons_count, false))
        cout << "A good Solution was found." << endl;
    else
        cout << "No valid Solutions were found from this start point." << endl;
    printpoints(points, point_count);

    double gradF[point_count] = {0};
    derivatives(pparameters, gradF, point_count, cons, cons_count);
    for (int i = 0; i < point_count; i++)
        cout << "GradF[" << i << "]: " << gradF[i] << endl;
    delete[] pparameters;
    return 0;
}
