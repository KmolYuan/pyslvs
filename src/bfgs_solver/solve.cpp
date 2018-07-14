/*
 * solve.cpp
 *
 *  Created on: May 4, 2009
 *  Author: Jonathan George
 *  Contributor: KmolYuan
 */

#include <iostream>
#include "solve.h"
#include <cmath>

#define _hypot hypot

///////////////////////////////////////////////////////////////////////
/// constraint defines (these make writing constraint equations easier
///////////////////////////////////////////////////////////////////////

#define P1_x         *cons[i].point1->x
#define P1_y         *cons[i].point1->y
#define P2_x         *cons[i].point2->x
#define P2_y         *cons[i].point2->y
#define L1_P1_x      *cons[i].line1->p1->x
#define L1_P1_y      *cons[i].line1->p1->y
#define L1_P2_x      *cons[i].line1->p2->x
#define L1_P2_y      *cons[i].line1->p2->y
#define L2_P1_x      *cons[i].line2->p1->x
#define L2_P1_y      *cons[i].line2->p1->y
#define L2_P2_x      *cons[i].line2->p2->x
#define L2_P2_y      *cons[i].line2->p2->y
#define C1_Center_x  *cons[i].circle1->center->x
#define C1_Center_y  *cons[i].circle1->center->y
#define C1_rad       *cons[i].circle1->rad
#define C2_Center_x  *cons[i].circle2->center->x
#define C2_Center_y  *cons[i].circle2->center->y
#define C2_rad       *cons[i].circle2->rad
#define A1_startA    *cons[i].arc1->startAngle
#define A1_endA      *cons[i].arc1->endAngle
#define A1_radius    *cons[i].arc1->rad
#define A1_Center_x  *cons[i].arc1->center->x
#define A1_Center_y  *cons[i].arc1->center->y
#define A2_startA    *cons[i].arc2->startAngle
#define A2_endA      *cons[i].arc2->endAngle
#define A2_radius    *cons[i].arc2->rad
#define A2_Center_x  *cons[i].arc2->center->x
#define A2_Center_y  *cons[i].arc2->center->y
#define A1_Start_x   (A1_Center_x + A1_radius * cos(A1_startA))
#define A1_Start_y   (A1_Center_y + A1_radius * sin(A1_startA))
#define A1_End_x     (A1_Center_x + A1_radius * cos(A1_endA))
#define A1_End_y     (A1_Center_y + A1_radius * sin(A1_endA))
#define A2_Start_x   (A1_Center_x + A2_radius * cos(A2_startA))
#define A2_Start_y   (A1_Center_y + A2_radius * sin(A2_startA))
#define A2_End_x     (A1_Center_x + A2_radius * cos(A2_endA))
#define A2_End_y     (A1_Center_y + A2_radius * sin(A2_endA))
#define length       *cons[i].parameter
#define distance     *cons[i].parameter
#define radius       *cons[i].parameter
#define angleP       *cons[i].parameter
#define quadIndex    *cons[i].parameter
#define Sym_P1_x     *cons[i].SymLine->p1->x
#define Sym_P1_y     *cons[i].SymLine->p1->y
#define Sym_P2_x     *cons[i].SymLine->p2->x
#define Sym_P2_y     *cons[i].SymLine->p2->y

using namespace std;


double calc(Constraint *, const int);

int solve(
    double **param_ptr,
    const int xLength,
    Constraint *cons,
    const int consLength,
    int const isFine
) {
    //Save the original parameters for later.
    double *origSolution = new double[xLength];
    for(int i = 0; i < xLength; i++)
        origSolution[i] = *param_ptr[i];

    double convergence = (isFine > 0) ? XConvergenceFine : XConvergenceRough;

    //integer to keep track of how many times calc is called
    int ftimes = 0;

    //Calculate Function at the starting point:
    double f0 = calc(cons, consLength);
    if(f0 < SmallF)
        return Succsess;
    ftimes++;

    //Calculate the gradient at the starting point:
    double *grad = new double[xLength]; //The gradient vector (1xn)
    double f1, f2, f3, alpha1, alpha2, alpha3, alphaStar;
    double norm = 0;
    double pert = f0 * PertMag;

    double first, second, temper; //The norm of the gradient vector
    for(int j = 0; j < xLength; j++) {
        temper = *param_ptr[j];
        *param_ptr[j] = temper - pert;
        first = calc(cons, consLength);
        *param_ptr[j] = temper + pert;
        second = calc(cons, consLength);
        grad[j] = 0.5 * (second - first) / pert;
        ftimes++;
        *param_ptr[j] = temper;
        norm = norm + grad[j] * grad[j];
    }
    norm = sqrt(norm);

    //Estimate the norm of N
    //Initialize N and calculate s
    double *s = new double[xLength]; //The current search direction
    double **N = new double *[xLength];
    for(int i = 0; i < xLength; i++)
        N[i] = new double[xLength]; //The estimate of the Hessian inverse
    for(int i = 0; i < xLength; i++) {
        for(int j = 0; j < xLength; j++) {
            if(i == j) {
                //N[i][j]=norm; //Calculate a scaled identity matrix as a Hessian inverse estimate
                //N[i][j]=grad[i]/(norm+.001);
                N[i][j] = 1;
                s[i] = -grad[i]; //Calculate the initial search vector
            } else
                N[i][j] = 0;
        }
    }
    double fnew = f0 + 1; //make fnew greater than fold
    double alpha = 1; //Initial search vector multiplier

    double *xold = new double[xLength]; //Storage for the previous design variables
    for(int i = 0; i < xLength; i++)
        xold[i] = *param_ptr[i];//Copy last values to xold

    ///////////////////////////////////////////////////////
    /// Start of line search
    ///////////////////////////////////////////////////////

    //Make the initial position alpha1
    alpha1 = 0;
    f1 = f0;

    //Take a step of alpha=1 as alpha2
    alpha2 = 1;
    for(int i = 0; i < xLength; i++)
        *param_ptr[i] = xold[i] + alpha2 * s[i];//calculate the new x

    f2 = calc(cons, consLength);
    ftimes++;

    //Take a step of alpha 3 that is 2*alpha2
    alpha3 = alpha * 2;
    for(int i = 0; i < xLength; i++)
        *param_ptr[i] = xold[i] + alpha3 * s[i];//calculate the new x

    f3 = calc(cons, consLength);
    ftimes++;

    //Now reduce or lengthen alpha2 and alpha3 until the minimum is
    //Bracketed by the triplet f1>f2<f3
    while(f2 > f1 || f2 > f3) {
        if(f2 > f1) {
            //If f2 is greater than f1 then we shorten alpha2 and alpha3 closer to f1
            //Effectively both are shortened by a factor of two.
            alpha3 = alpha2;
            f3 = f2;
            alpha2 = alpha2 / 2;
            for(int i = 0; i < xLength; i++)
                *param_ptr[i] = xold[i] + alpha2 * s[i];//calculate the new x
            f2 = calc(cons, consLength);
            ftimes++;
        } else if(f2 > f3) {
            //If f2 is greater than f3 then we length alpah2 and alpha3 closer to f1
            //Effectively both are lengthened by a factor of two.
            alpha2 = alpha3;
            f2 = f3;
            alpha3 = alpha3 * 2;
            for(int i = 0; i < xLength; i++)
                *param_ptr[i] = xold[i] + alpha3 * s[i];//calculate the new x
            f3 = calc(cons, consLength);
            ftimes++;
        }
    }
    // get the alpha for the minimum f of the quadratic approximation
    alphaStar = alpha2 + (alpha2 - alpha1) * (f1 - f3) / (3*(f1 - 2 * f2 + f3));

    //Guarantee that the new alphaStar is within the bracket
    if(alphaStar > alpha3 || alphaStar < alpha1)
        alphaStar = alpha2;
    if(alphaStar != alphaStar)
        alphaStar = 0.001;//Fix nan problem

    /// Set the values to alphaStar
    for(int i = 0; i < xLength; i++)
        *param_ptr[i] = xold[i] + alphaStar * s[i];//calculate the new x
    fnew = calc(cons, consLength);
    ftimes++;

    /////////////////////////////////////
    ///end of line search
    /////////////////////////////////////

    double *deltaX = new double[xLength];
    double *gradnew = new double[xLength];
    double *gamma = new double[xLength];
    double *gammatDotN = new double[xLength];
    double **FirstSecond = new double *[xLength];
    double **deltaXDotGammatDotN = new double *[xLength];
    double **gammatDotDeltaXt = new double *[xLength];
    double **NDotGammaDotDeltaXt = new double *[xLength];
    for(int i = 0; i < xLength; i++) {
        FirstSecond[i] = new double[xLength];
        deltaXDotGammatDotN[i] = new double[xLength];
        gammatDotDeltaXt[i] = new double[xLength];
        NDotGammaDotDeltaXt[i] = new double[xLength];
    }

    ///Calculate deltaX
    for(int i = 0; i < xLength; i++)
        deltaX[i] = *param_ptr[i] - xold[i]; //Calculate the difference in x for the Hessian update

    int iterations = 1;
    double maxIterNumber = MaxIterations * xLength;
    double deltaXnorm = 1;
    int steps;
    double deltaXtDotGamma, gammatDotNDotGamma, firstTerm, bottom;
    while(deltaXnorm > convergence && fnew > SmallF && iterations < maxIterNumber) {
        //////////////////////////////////////////////////////////////////////
        ///Start of main loop!!!!
        //////////////////////////////////////////////////////////////////////
        bottom = 0;
        deltaXtDotGamma = 0;
        pert = fnew * PertMag;
        if(pert < PertMin)
            pert = PertMin;
        for(int i = 0; i < xLength; i++) {
            //Calculate the new gradient vector
            temper = *param_ptr[i];
            *param_ptr[i] = temper - pert;
            first = calc(cons, consLength);
            *param_ptr[i] = temper + pert;
            second = calc(cons, consLength);
            gradnew[i] = 0.5 * (second - first) / pert;
            ftimes++;
            *param_ptr[i] = temper;
            //Calculate the change in the gradient
            gamma[i] = gradnew[i] - grad[i];
            bottom += deltaX[i] * gamma[i];

            deltaXtDotGamma += deltaX[i] * gamma[i];
        }

        //make sure that bottom is never 0
        if (bottom == 0)
            bottom = 0.0000000001;

        //calculate all (1xn).(nxn)
        for(int i = 0; i < xLength; i++) {
            gammatDotN[i] = 0;
            for(int j = 0; j < xLength; j++)
                gammatDotN[i] += gamma[j] * N[i][j];//This is gammatDotN transpose
        }

        //calculate all (1xn).(nx1)
        gammatDotNDotGamma = 0;
        for(int i = 0; i < xLength; i++)
            gammatDotNDotGamma += gammatDotN[i] * gamma[i];

        //Calculate the first term

        firstTerm = 0;
        firstTerm = 1 + gammatDotNDotGamma / bottom;

        //Calculate all (nx1).(1xn) matrices
        for(int i = 0; i < xLength; i++) {
            for(int j = 0; j < xLength; j++) {
                FirstSecond[i][j] = deltaX[j] * deltaX[i] / bottom * firstTerm;
                deltaXDotGammatDotN[i][j] = deltaX[i] * gammatDotN[j];
                gammatDotDeltaXt[i][j] = gamma[i] * deltaX[j];
            }
        }

        //Calculate all (nxn).(nxn) matrices
        for(int i = 0; i < xLength; i++)
            for(int j = 0; j < xLength; j++) {
                NDotGammaDotDeltaXt[i][j] = 0;
                for(int k = 0; k < xLength; k++)
                    NDotGammaDotDeltaXt[i][j] += N[i][k] * gammatDotDeltaXt[k][j];
            }

        //Now calculate the BFGS update on N
        //cout<<"N:"<<endl;
        for(int i = 0; i < xLength; i++)
            for(int j = 0; j < xLength; j++)
                N[i][j] = N[i][j] + FirstSecond[i][j] -
                        (deltaXDotGammatDotN[i][j] + NDotGammaDotDeltaXt[i][j]) / bottom;

        //Calculate s
        for(int i = 0; i < xLength; i++) {
            s[i] = 0;
            for(int j = 0; j < xLength; j++)
                s[i] += -N[i][j] * gradnew[j];
        }

        alpha = 1; //Initial search vector multiplier

        //copy newest values to the xold
        for(int i = 0; i < xLength; i++)
            xold[i] = *param_ptr[i]; //Copy last values to xold
        steps = 0;

        ///////////////////////////////////////////////////////
        /// Start of line search
        ///////////////////////////////////////////////////////

        //Make the initial position alpha1
        alpha1 = 0;
        f1 = fnew;

        //Take a step of alpha=1 as alpha2
        alpha2 = 1;
        for(int i = 0; i < xLength; i++)
            *param_ptr[i] = xold[i] + alpha2 * s[i]; //calculate the new x
        f2 = calc(cons, consLength);
        ftimes++;

        //Take a step of alpha 3 that is 2*alpha2
        alpha3 = alpha2 * 2;
        for(int i = 0; i < xLength; i++)
            *param_ptr[i] = xold[i] + alpha3 * s[i]; //calculate the new x
        f3 = calc(cons, consLength);
        ftimes++;

        //Now reduce or lengthen alpha2 and alpha3 until the minimum is
        //Bracketed by the triplet f1>f2<f3
        steps = 0;
        while(f2 > f1 || f2 > f3) {
            if(f2 > f1) {
                //If f2 is greater than f1 then we shorten alpha2 and alpha3 closer to f1
                //Effectively both are shortened by a factor of two.
                alpha3 = alpha2;
                f3 = f2;
                alpha2 /= 2;
                for(int i = 0; i < xLength; i++)
                    *param_ptr[i] = xold[i] + alpha2 * s[i]; //calculate the new x
                f2 = calc(cons, consLength);
                ftimes++;
            } else if(f2 > f3) {
                //If f2 is greater than f3 then we length alpah2 and alpha3 closer to f1
                //Effectively both are lengthened by a factor of two.
                alpha2 = alpha3;
                f2 = f3;
                alpha3 = alpha3 * 2;
                for(int i = 0; i < xLength; i++)
                    *param_ptr[i] = xold[i] + alpha3 * s[i]; //calculate the new x
                f3 = calc(cons, consLength);
                ftimes++;
            }
            steps++;
        }

        // get the alpha for the minimum f of the quadratic approximation
        alphaStar = alpha2 + (alpha2 - alpha1) * (f1 - f3) / (3 * (f1 - 2 * f2 + f3));

        //Guarantee that the new alphaStar is within the bracket
        if(alphaStar >= alpha3 || alphaStar <= alpha1)
            alphaStar = alpha2;
        if(alphaStar != alphaStar)
            alphaStar = 0;

        /// Set the values to alphaStar
        for(int i = 0; i < xLength; i++)
            *param_ptr[i] = xold[i] + alphaStar * s[i]; //calculate the new x
        fnew = calc(cons, consLength);
        ftimes++;

        /////////////////////////////////////
        ///end of line search
        ////////////////////////////////////

        deltaXnorm = 0;
        for(int i = 0; i < xLength; i++) {
            //Calculate the difference in x for the hessian update
            deltaX[i] = *param_ptr[i] - xold[i];
            deltaXnorm += deltaX[i] * deltaX[i];
            grad[i] = gradnew[i];
        }
        deltaXnorm = sqrt(deltaXnorm);
        iterations++;
        /////////////////////////////////////////////////////////////
        ///End of Main loop
        /////////////////////////////////////////////////////////////
    }

#ifdef DEBUG
    cout << "Fnew: " << fnew << endl;
    cout << "Number of Iterations: " << iterations << endl;
    cout << "Number of function calls: " << ftimes << endl;
#endif

    delete s;
    delete [] N;
    delete [] FirstSecond;
    delete [] deltaXDotGammatDotN;
    delete [] gammatDotDeltaXt;
    delete [] NDotGammaDotDeltaXt;
    delete origSolution;

    delete grad;
    delete xold;
    delete gammatDotN;

    ///End of function
    if(fnew < ((isFine == 1) ? ValidSolutionFine : ValidSoltuionRough))
        return Succsess;
    else {
        //Replace the bad numbers with the last result
        for(int i = 0; i < xLength; i++)
            *param_ptr[i] = origSolution[i];
        return NoSolution;
    }
}


double calc(Constraint *cons, const int consLength) {
    double error = 0;
    double temp, temp2, dx, dy, m, n, Ex, Ey, rad1, rad2, t, dx2, dy2, hyp1, hyp2;
    for(int i = 0; i < consLength; i++) {
        switch(cons[i].type) {
        case PointOnPoint:
            //Hopefully avoid this constraint, make coincident points use the same parameters
            dx = P1_x - P2_x;
            dy = P1_y - P2_y;
            error += dx * dx + dy * dy;
            break;

        case P2PDistance:
            temp = _hypot(P2_x - P1_x, P2_y - P1_y) - distance;
            error += temp * temp * 10000;
            break;

        case P2PDistanceVert:
            temp = _hypot(0, P2_y - P1_y) - distance;
            error += temp * temp * 10000;
            break;

        case P2PDistanceHorz:
            temp = _hypot(P2_x - P1_x, 0) - distance;
            error += temp * temp * 10000;
            break;

        case PointOnLine:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;

            m = dy / dx;
            n = dx / dy;
            if(m <= 1 && m >= -1) {
                //Calculate the expected y point given the x coordinate of the point
                Ey = L1_P1_y + m * (P1_x - L1_P1_x);
                error += (Ey - P1_y) * (Ey - P1_y);
            } else {
                //Calculate the expected x point given the y coordinate of the point
                Ex = L1_P1_x + n * (P1_y - L1_P1_y);
                error += (Ex - P1_x) * (Ex - P1_x);
            }
            break;

        case P2LDistance:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;

            t = -(L1_P1_x * dx - P1_x * dx + L1_P1_y * dy - P1_y * dy) / (dx * dx + dy * dy);
            temp = _hypot(P1_x - (L1_P1_x + dx * t), P1_y - (L1_P1_y + dy * t)) - distance;
            error += temp * temp / 10;
            break;

        case P2LDistanceVert:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;

            t = (P1_x - L1_P1_x) / dx;
            temp = fabs(P1_y - (L1_P1_y + dy * t)) - distance;
            error += temp * temp;
            break;

        case P2LDistanceHorz:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;

            t = (P1_y - L1_P1_y) / dy;
            temp = fabs(P1_x - (L1_P1_x + dx * t)) - distance;
            error += temp * temp / 10;
            break;

        case Vertical:
            dx = L1_P2_x - L1_P1_x;
            error += dx * dx * 10000;
            break;

        case Horizontal:
            dy = L1_P2_y - L1_P1_y;
            error += dy * dy * 10000;
            break;

        case TangentToCircle:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;
            hyp1 = _hypot(dx, dy);
            {
                //Calculate the expected tangent intersection points
                double Rpx = C1_Center_x - dy / hyp1 * C1_rad;
                double Rpy = C1_Center_y + dx / hyp1 * C1_rad;
                double RpxN = C1_Center_x + dy / hyp1 * C1_rad;
                double RpyN = C1_Center_y - dx / hyp1 * C1_rad;

                double error1 = (-dy * Rpx +
                                 dx * Rpy +
                                 (L1_P1_x * L1_P2_y - L1_P2_x * L1_P1_y)) / hyp1;
                double error2 = (-dy * RpxN +
                                 dx * RpyN +
                                 (L1_P1_x * L1_P2_y - L1_P2_x * L1_P1_y)) / hyp1;
                error += (error1 < error2) ? error1 * error1 : error2 * error2;
            }
            break;

        case TangentToArc:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;
            t = -(L1_P1_x * dx -
                  A1_Center_x * dx +
                  L1_P1_y * dy -
                  A1_Center_y * dy) / (dx * dx + dy * dy);
            Ex = L1_P1_x + dx * t;
            Ey = L1_P1_y + dy * t;
            temp = (A1_Center_x - Ex) *
                    (A1_Center_x - Ex) +
                    (A1_Center_y - Ey) *
                    (A1_Center_y - Ey) -
                    (A1_Center_x - A1_Start_x) * (A1_Center_x - A1_Start_x) -
                    (A1_Center_y - A1_Start_y) * (A1_Center_y - A1_Start_y);
            error += temp * temp;
            break;

        case ArcRules:
            dx = A1_End_x * A1_End_x;
            dy = A1_End_y * A1_End_y;
            rad1 = A1_Start_x * A1_Start_x;
            rad2 = A1_Start_y * A1_Start_y;
            temp = -2 * A1_Center_x * A1_End_x + dx -
                    2 * A1_Center_y * A1_End_y + dy +
                    2 * A1_Center_x * A1_Start_x - rad1 +
                    2 * A1_Center_y * A1_Start_y - rad2;
            error += temp * temp / (4. * dx + dy -
                                    2 * A1_End_x * A1_Start_x + rad1 -
                                    2 * A1_End_y * A1_Start_y + rad2);
            break;

        case LineLength:
            temp = _hypot(L1_P2_x - L1_P1_x, L1_P2_y - L1_P1_y) - length;
            error += temp * temp * 10000;
            break;

        case EqualLegnth:
            temp = _hypot(L1_P2_x - L1_P1_x, L1_P2_y - L1_P1_y) -
                   _hypot(L2_P2_x - L2_P1_x, L2_P2_y - L2_P1_y);
            error += temp * temp * 10000;
            break;

        case ArcRadius:
            rad1 = _hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
            rad2 = _hypot(A1_Center_x - A1_End_x, A1_Center_y - A1_End_y);
            temp= rad1 - radius;
            error += temp * temp;
            break;

        case EqualRadiusArcs:
            rad1 = _hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
            rad2 = _hypot(A2_Center_x - A2_Start_x, A2_Center_y - A2_Start_y);
            temp = rad1 - rad2;
            error += temp * temp;
            break;

        case EqualRadiusCircles:
            temp = C1_rad - C2_rad;
            error += temp * temp;
            break;

        case EqualRadiusCircArc:
            rad1 = _hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
            temp = rad1 - C1_rad;
            error += temp * temp;
            break;

        case ConcentricArcs:
            temp = _hypot(A1_Center_x - A2_Center_x, A1_Center_y - A2_Center_y);
            error += temp * temp;
            break;

        case ConcentricCircles:
            temp = _hypot(C1_Center_x - C2_Center_x, C1_Center_y - C2_Center_y);
            error += temp * temp;
            break;

        case ConcentricCircArc:
            temp = _hypot(A1_Center_x - C1_Center_x, A1_Center_y - C1_Center_y);
            error += temp * temp;
            break;

        case CircleRadius:
            error += (C1_rad - radius) * (C1_rad - radius);
            break;

        case InternalAngle:
            temp = atan2(L2_P2_y - L2_P1_y, L2_P2_x - L2_P1_x) -
                   atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) -
                   angleP;
            error += temp * temp;
            break;

        case ExternalAngle:
            temp = atan2(L2_P2_y - L2_P1_y, L2_P2_x - L2_P1_x) -
                   atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) -
                   (M_PI - angleP);
            error += temp * temp;
            break;

        case LineInternalAngle:
            temp = sin(atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) - angleP);
            error += temp * temp;
            break;

        case LineExternalAngle:
            temp = sin(atan2(L1_P2_y - L1_P1_y, L1_P2_x - L1_P1_x) - (M_PI - angleP));
            error += temp * temp;
            break;

        case Perpendicular:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;
            dx2 = L2_P2_x - L2_P1_x;
            dy2 = L2_P2_y - L2_P1_y;

            hyp1 = _hypot(dx, dy);
            hyp2 = _hypot(dx2, dy2);

            dx = dx / hyp1;
            dy = dy / hyp1;
            dx2 = dx2 / hyp2;
            dy2 = dy2 / hyp2;

            temp = dx * dx2 + dy * dy2;
            error += temp * temp;
            break;

        case Parallel:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;
            dx2 = L2_P2_x - L2_P1_x;
            dy2 = L2_P2_y - L2_P1_y;

            hyp1 = _hypot(dx, dy);
            hyp2 = _hypot(dx2, dy2);

            dx = dx / hyp1;
            dy = dy / hyp1;
            dx2 = dx2 / hyp2;
            dy2 = dy2 / hyp2;

            temp = dy * dx2 - dx * dy2;
            error += temp * temp;
            break;

        case Colinear:
            dx = L1_P2_x - L1_P1_x;
            dy = L1_P2_y - L1_P1_y;

            m = dy / dx;
            n = dx / dy;
            // Calculate the error between the expected intersection point
            // and the true point of the second lines two end points on the
            // first line
            if(m <= 1 && m > -1) {
                //Calculate the expected y point given the x coordinate of the point
                Ey = L1_P1_y + m * (L2_P1_x - L1_P1_x);
                error += (Ey - L2_P1_y) * (Ey - L2_P1_y);

                Ey = L1_P1_y + m * (L2_P2_x - L1_P1_x);
                error += (Ey - L2_P2_y) * (Ey - L2_P2_y);
            } else {
                //Calculate the expected x point given the y coordinate of the point
                Ex = L1_P1_x + n * (L2_P1_y - L1_P1_y);
                error += (Ex - L2_P1_x) * (Ex - L2_P1_x);

                Ex = L1_P1_x + n * (L2_P2_y - L1_P1_y);
                error += (Ex - L2_P2_x) * (Ex - L2_P2_x);
            }
            break;

        case PointOnCircle:
            //see what the current radius to the point is
            rad1 = _hypot(C1_Center_x - P1_x, C1_Center_y - P1_y);
            //Compare this radius to the radius of the circle, return the error squared
            temp = rad1 - C1_rad;
            error += temp * temp;
            break;

        case PointOnArc:
            //see what the current radius to the point is
            rad1 = _hypot(A1_Center_x - P1_x, A1_Center_y - P1_y);
            rad2 = _hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
            //Compare this radius to the radius of the circle, return the error squared
            temp = rad1 - rad2;
            error += temp * temp;
            break;

        case PointOnLineMidpoint:
            Ex = (L1_P1_x + L1_P2_x) / 2 - P1_x;
            Ey = (L1_P1_y + L1_P2_y) / 2 - P1_y;
            error += Ex * Ex + Ey * Ey;
            break;

        case PointOnArcMidpoint:
            rad1 = _hypot(A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
            temp = atan2(A1_Start_y - A1_Center_y, A1_Start_x - A1_Center_x);
            temp2 = atan2(A1_End_y - A1_Center_y, A1_End_x - A1_Center_x);
            Ex = A1_Center_x + rad1 * cos((temp2 + temp) / 2);
            Ey = A1_Center_y + rad1 * sin((temp2 + temp) / 2);
            temp = Ex - P1_x;
            temp2 = Ey - P1_y;
            error += temp * temp + temp2 * temp2;
            break;

        case PointOnCircleQuad:
            Ex = C1_Center_x;
            Ey = C1_Center_y;
            switch((int)quadIndex) {
            case 1:
                Ey += C1_rad;
                break;
            case 2:
                Ex -= C1_rad;
                break;
            case 3:
                Ey -= C1_rad;
                break;
            default:
                Ex += C1_rad;
                break;
            }
            temp = Ex - P1_x;
            temp2 = Ey - P1_y;
            error += temp * temp + temp2 * temp2;
            break;

        case SymmetricPoints:
            dx = Sym_P2_x - Sym_P1_x;
            dy = Sym_P2_y - Sym_P1_y;
            t = -(dy * P1_x -
                  dx * P1_y -
                  dy * Sym_P1_x +
                  dx * Sym_P1_y) / (dx * dx + dy * dy);
            Ex = P1_x + dy * t * 2;
            Ey = P1_y - dx * t * 2;
            temp = Ex - P2_x;
            temp2 = Ey - P2_y;
            error += temp * temp + temp2 * temp2;
            break;

        case SymmetricLines:
            dx = Sym_P2_x - Sym_P1_x;
            dy = Sym_P2_y - Sym_P1_y;
            hyp1 = dx * dx + dy * dy;
            m = -dy * Sym_P1_x + dx * Sym_P1_y;

            t = -(dy * L1_P1_x - dx * L1_P1_y + m) / hyp1;
            Ex = L1_P1_x + dy * t * 2;
            Ey = L1_P1_y - dx * t * 2;
            temp = Ex - L2_P1_x;
            temp2 = Ey - L2_P1_y;
            error += temp * temp + temp2 * temp2;
            t = -(dy * L1_P2_x - dx * L1_P2_y + m) / hyp1;
            Ex = L1_P2_x + dy * t * 2;
            Ey = L1_P2_y - dx * t * 2;
            temp = Ex - L2_P2_x;
            temp2 = Ey - L2_P2_y;
            error += temp * temp + temp2 * temp2;
            break;

        case SymmetricCircles:
            dx = Sym_P2_x - Sym_P1_x;
            dy = Sym_P2_y - Sym_P1_y;
            t = -(dy * C1_Center_x -
                  dx * C1_Center_y -
                  dy * Sym_P1_x +
                  dx * Sym_P1_y) / (dx * dx + dy * dy);
            Ex = C1_Center_x + dy * t * 2;
            Ey = C1_Center_y - dx * t * 2;
            temp = Ex - C2_Center_x;
            temp2 = Ey - C2_Center_y;
            error += temp * temp + temp2 * temp2;
            temp = C1_rad - C2_rad;
            error += temp * temp;
            break;

        case SymmetricArcs:
            dx = Sym_P2_x - Sym_P1_x;
            dy = Sym_P2_y - Sym_P1_y;
            hyp1 = dx * dx + dy * dy;
            m = - dy * Sym_P1_x + dx * Sym_P1_y;

            t = -(dy * A1_Start_x - dx * A1_Start_y + m) / hyp1;
            Ex = A1_Start_x + dy * t * 2;
            Ey = A1_Start_y - dx * t * 2;
            temp = Ex - A2_Start_x;
            temp2 = Ey - A2_Start_y;
            error += temp * temp + temp2 * temp2;
            t = -(dy * A1_End_x - dx * A1_End_y + m) / hyp1;
            Ex = A1_End_x + dy * t * 2;
            Ey = A1_End_y - dx * t * 2;
            temp = Ex - A2_End_x;
            temp2 = Ey - A2_End_y;
            error += temp * temp + temp2 * temp2;
            t = -(dy * A1_Center_x - dx * A1_Center_y + m) / hyp1;
            Ex = A1_Center_x + dy * t * 2;
            Ey = A1_Center_y - dx * t * 2;
            temp = Ex - A2_Center_x;
            temp2 = Ey - A2_Center_y;
            error += temp * temp + temp2 * temp2;
            break;
        }
    }
    return error;
}
