/*
 * ============================================================================
 * Name        : SolverPointers.cpp
 * Author      : Jonathan George
 * Contributor : KmolYuan
 *============================================================================*/

#include <iostream>
#include "solve.h"

using namespace std;

//Show coordinates.
template<size_t N> inline void printpoints(Point (&points)[N]) {
    for(int i = 0; i < (int)N; i++)
        cout << "Point " << i << ": (" << *points[i].x << ", " << *points[i].y << ")" << endl;
}


int main(void) {
    //Input a parameter list.
    double parameters[] = {
        0, 0,
        5, 0,
        6, 5,
        6, 5,
    };
    const int param_count = sizeof(parameters) / sizeof(*parameters);
    //Make a list of pointers of parameters.
    double *pparameters[param_count];
    for(int i = 0; i < param_count; i++)
        pparameters[i] = &parameters[i];
    //Input a constant parameter list.
    double constants[] = {30, 10, 24,};

    //Create geometric objects and constraints with pointers.
    Point points[] = {
        {pparameters[0], pparameters[1]},
        {pparameters[2], pparameters[3]},
        {pparameters[4], pparameters[5]},
        {pparameters[6], pparameters[7]},
        {&constants[0], &constants[1]},
    };
    const int point_count = sizeof(points) / sizeof(*points);

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

    printpoints<point_count>(points);

    //Solve
    if (solve(pparameters, param_count, cons, cons_count, Rough) == Succsess)
        cout << "A good Solution was found." << endl;
    else
        cout << "No valid Solutions were found from this start point." << endl;

    printpoints<point_count>(points);

    double gradF[point_count] = {0};
    derivatives(pparameters, gradF, point_count, cons, cons_count);
    for(int i = 0; i < point_count; i++)
        cout << "GradF[" << i << "]: " << gradF[i] << endl;

    return 0;
}
