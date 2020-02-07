/*
 * ============================================================================
 * Name        : SolverPointers.cpp
 * Author      : Jonathan George
 * Contributor : KmolYuan
 *============================================================================*/

#include <iostream>
#include "solve.h"

using namespace std;

// Show coordinates
void printpoints(Point *points, size_t N) {
    for (size_t i = 0; i < N; i++)
        cout << "Point " << i << ": (" << *points[i].x << ", " << *points[i].y
             << ")" << endl;
}

int main() {
    // Input a parameter list
    double parameters[] = {0, 0, 5, 0, 6, 5, 6, 5};
    const size_t param_count = sizeof(parameters) / sizeof(*parameters);
    // Make a list of pointers of parameters
    auto pparameters = new double *[param_count];
    for (int i = 0; i < param_count; i++)
        pparameters[i] = &parameters[i];
    // Input a constant parameter list
    double constants[] = {30, 10, 24};

    // Create geometric objects and constraints with pointers
    Point points[] = {
        {pparameters[0], pparameters[1]}, {pparameters[2], pparameters[3]},
        {pparameters[4], pparameters[5]}, {pparameters[6], pparameters[7]},
        {&constants[0], &constants[1]},
    };
    const size_t point_count = sizeof(points) / sizeof(*points);

    Line lines[] = {
        {&points[0], &points[1]},
        {&points[1], &points[2]},
    };

    Constraint cons[] = {
        HorizontalConstraint(&lines[0]),
        PointOnPointConstraint(&points[2], &points[4]),
        PointOnPointConstraint(&points[3], &points[4]),
        P2PDistanceConstraint(&points[1], &points[2], &constants[2]),
    };
    const int cons_count = sizeof(cons) / sizeof(*cons);
    printpoints(points, point_count);

    // Solve
    if (solve(pparameters, param_count, cons, cons_count, Rough) == Success)
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
