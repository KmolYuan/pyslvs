/*
 * solve.cpp
 *
 *  Created on: May 4, 2009
 *  Author: Jonathan George
 *  Contributor: KmolYuan
 */

#include <cmath>
#ifdef DEBUG
#include <iostream>
#endif
#include "calc.h"
#include "solve.h"

using namespace std;

int solve(double **param_ptr, const size_t xLength, Constraint *cons,
          const size_t consLength, const bool isFine) {
    // Save the original parameters for later.
    auto origSolution = new double[xLength];
    for (size_t i = 0; i < xLength; i++)
        origSolution[i] = *param_ptr[i];

    double convergence = isFine ? XConvergenceFine : XConvergenceRough;

    // integer to keep track of how many times calc is called
    int ftimes = 0;

    // Calculate Function at the starting point:
    double f0 = calc(cons, consLength);
    if (f0 < SmallF)
        return Success;
    ftimes++;

    // Calculate the gradient at the starting point:
    // The gradient vector (1xn)
    auto grad = new double[xLength];
    double f1, f2, f3, alpha1, alpha2, alpha3, alphaStar;
    double pert = f0 * PertMag;

    // The norm of the gradient vector
    for (size_t j = 0; j < xLength; j++) {
        double tmp = *param_ptr[j];
        *param_ptr[j] = tmp - pert;
        double first = calc(cons, consLength);
        *param_ptr[j] = tmp + pert;
        double second = calc(cons, consLength);
        grad[j] = 0.5 * (second - first) / pert;
        *param_ptr[j] = tmp;
        ftimes++;
    }

    // Estimate the norm of N
    // Initialize N and calculate s
    // The current search direction
    auto s = new double[xLength];
    auto N = new double *[xLength];
    for (size_t i = 0; i < xLength; i++)
        // The estimate of the Hessian inverse
        N[i] = new double[xLength];
    for (size_t i = 0; i < xLength; i++)
        for (size_t j = 0; j < xLength; j++) {
            if (i == j) {
                N[i][j] = 1;
                // Calculate the initial search vector
                s[i] = -grad[i];
                continue;
            }
            N[i][j] = 0;
        }

    // Initial search vector multiplier
    double alpha = 1;

    // Storage for the previous design variables
    auto xold = new double[xLength];
    for (size_t i = 0; i < xLength; i++)
        // Copy last values to xold
        xold[i] = *param_ptr[i];

    ///////////////////////////////////////////////////////
    /// Start of line search
    ///////////////////////////////////////////////////////

    // Make the initial position alpha1
    alpha1 = 0;
    f1 = f0;

    // Take a step of alpha=1 as alpha2
    alpha2 = 1;
    for (size_t i = 0; i < xLength; i++)
        // calculate the new x
        *param_ptr[i] = xold[i] + alpha2 * s[i];

    f2 = calc(cons, consLength);
    ftimes++;

    // Take a step of alpha 3 that is 2*alpha2
    alpha3 = alpha * 2;
    for (size_t i = 0; i < xLength; i++)
        // calculate the new x
        *param_ptr[i] = xold[i] + alpha3 * s[i];

    f3 = calc(cons, consLength);
    ftimes++;

    // Now reduce or lengthen alpha2 and alpha3 until the minimum is
    // Bracketed by the triplet f1>f2<f3
    while (f2 > f1 || f2 > f3)
        if (f2 > f1) {
            // If f2 is greater than f1 then we shorten alpha2 and alpha3 closer
            // to f1 Effectively both are shortened by a factor of two.
            alpha3 = alpha2;
            f3 = f2;
            alpha2 = alpha2 / 2;
            for (size_t i = 0; i < xLength; i++)
                // calculate the new x
                *param_ptr[i] = xold[i] + alpha2 * s[i];
            f2 = calc(cons, consLength);
            ftimes++;
        } else if (f2 > f3) {
            // If f2 is greater than f3 then we length alpah2 and alpha3 closer
            // to f1 Effectively both are lengthened by a factor of two.
            alpha2 = alpha3;
            f2 = f3;
            alpha3 = alpha3 * 2;
            for (size_t i = 0; i < xLength; i++)
                // calculate the new x
                *param_ptr[i] = xold[i] + alpha3 * s[i];
            f3 = calc(cons, consLength);
            ftimes++;
        }
    // get the alpha for the minimum f of the quadratic approximation
    alphaStar =
        alpha2 + (alpha2 - alpha1) * (f1 - f3) / (3 * (f1 - 2 * f2 + f3));

    // Guarantee that the new alphaStar is within the bracket
    if (alphaStar > alpha3 || alphaStar < alpha1)
        alphaStar = alpha2;
    if (alphaStar != alphaStar)
        // Fix nan problem
        alphaStar = 0.001;

    /// Set the values to alphaStar
    for (size_t i = 0; i < xLength; i++)
        // calculate the new x
        *param_ptr[i] = xold[i] + alphaStar * s[i];
    double fnew = calc(cons, consLength);
    ftimes++;

    /////////////////////////////////////
    /// end of line search
    /////////////////////////////////////

    auto deltaX = new double[xLength];
    auto gradnew = new double[xLength];
    auto gamma = new double[xLength];
    auto gammatDotN = new double[xLength];
    auto FirstSecond = new double *[xLength];
    auto deltaXDotGammatDotN = new double *[xLength];
    auto gammatDotDeltaXt = new double *[xLength];
    auto NDotGammaDotDeltaXt = new double *[xLength];
    for (size_t i = 0; i < xLength; i++) {
        FirstSecond[i] = new double[xLength];
        deltaXDotGammatDotN[i] = new double[xLength];
        gammatDotDeltaXt[i] = new double[xLength];
        NDotGammaDotDeltaXt[i] = new double[xLength];
    }

    // Calculate deltaX
    for (size_t i = 0; i < xLength; i++)
        // Calculate the difference in x for the Hessian update
        deltaX[i] = *param_ptr[i] - xold[i];

    size_t iterations = 1;
    double maxIterNumber = MaxIterations * xLength;
    double deltaXnorm = 1;
    double deltaXtDotGamma, gammatDotNDotGamma, firstTerm, bottom;
    while (deltaXnorm > convergence && fnew > SmallF &&
           iterations < maxIterNumber) {
        //////////////////////////////////////////////////////////////////////
        /// Start of main loop!!!!
        //////////////////////////////////////////////////////////////////////
        bottom = 0;
        deltaXtDotGamma = 0;
        pert = fnew * PertMag;
        if (pert < PertMin)
            pert = PertMin;
        for (size_t i = 0; i < xLength; i++) {
            // Calculate the new gradient vector
            double tmp = *param_ptr[i];
            *param_ptr[i] = tmp - pert;
            double first = calc(cons, consLength);
            *param_ptr[i] = tmp + pert;
            double second = calc(cons, consLength);
            gradnew[i] = 0.5 * (second - first) / pert;
            ftimes++;
            *param_ptr[i] = tmp;
            // Calculate the change in the gradient
            gamma[i] = gradnew[i] - grad[i];
            bottom += deltaX[i] * gamma[i];

            deltaXtDotGamma += deltaX[i] * gamma[i];
        }

        // make sure that bottom is never 0
        if (bottom == 0.)
            bottom = 0.0000000001;

        // calculate all (1xn).(nxn)
        for (size_t i = 0; i < xLength; i++) {
            gammatDotN[i] = 0;
            for (size_t j = 0; j < xLength; j++)
                // This is gammatDotN transpose
                gammatDotN[i] += gamma[j] * N[i][j];
        }

        // calculate all (1xn).(nx1)
        gammatDotNDotGamma = 0;
        for (size_t i = 0; i < xLength; i++)
            gammatDotNDotGamma += gammatDotN[i] * gamma[i];

        // Calculate the first term
        firstTerm = 1 + gammatDotNDotGamma / bottom;

        // Calculate all (nx1).(1xn) matrices
        for (size_t i = 0; i < xLength; i++)
            for (size_t j = 0; j < xLength; j++) {
                FirstSecond[i][j] = deltaX[j] * deltaX[i] / bottom * firstTerm;
                deltaXDotGammatDotN[i][j] = deltaX[i] * gammatDotN[j];
                gammatDotDeltaXt[i][j] = gamma[i] * deltaX[j];
            }

        // Calculate all (nxn).(nxn) matrices
        for (size_t i = 0; i < xLength; i++)
            for (size_t j = 0; j < xLength; j++) {
                NDotGammaDotDeltaXt[i][j] = 0;
                for (size_t k = 0; k < xLength; k++)
                    NDotGammaDotDeltaXt[i][j] +=
                        N[i][k] * gammatDotDeltaXt[k][j];
            }

        // Now calculate the BFGS update on N
        for (size_t i = 0; i < xLength; i++)
            for (size_t j = 0; j < xLength; j++)
                N[i][j] =
                    N[i][j] + FirstSecond[i][j] -
                    (deltaXDotGammatDotN[i][j] + NDotGammaDotDeltaXt[i][j]) /
                        bottom;

        // Calculates
        for (size_t i = 0; i < xLength; i++) {
            s[i] = 0;
            for (size_t j = 0; j < xLength; j++)
                s[i] += -N[i][j] * gradnew[j];
        }

        // copy newest values to the xold
        for (size_t i = 0; i < xLength; i++)
            xold[i] = *param_ptr[i];  // Copy last values to xold

        ///////////////////////////////////////////////////////
        /// Start of line search
        ///////////////////////////////////////////////////////

        // Make the initial position alpha1
        alpha1 = 0;
        f1 = fnew;

        // Take a step of alpha=1 as alpha2
        alpha2 = 1;
        for (size_t i = 0; i < xLength; i++)
            // calculate the new x
            *param_ptr[i] = xold[i] + alpha2 * s[i];
        f2 = calc(cons, consLength);
        ftimes++;

        // Take a step of alpha 3 that is 2*alpha2
        alpha3 = alpha2 * 2;
        for (size_t i = 0; i < xLength; i++)
            // calculate the new x
            *param_ptr[i] = xold[i] + alpha3 * s[i];
        f3 = calc(cons, consLength);
        ftimes++;

        // Now reduce or lengthen alpha2 and alpha3 until the minimum is
        // Bracketed by the triplet f1>f2<f3
        while (f2 > f1 || f2 > f3)
            if (f2 > f1) {
                // If f2 is greater than f1 then we shorten alpha2 and alpha3
                // closer to f1 Effectively both are shortened by a factor of
                // two.
                alpha3 = alpha2;
                f3 = f2;
                alpha2 /= 2;
                for (size_t i = 0; i < xLength; i++)
                    // calculate the new x
                    *param_ptr[i] = xold[i] + alpha2 * s[i];
                f2 = calc(cons, consLength);
                ftimes++;
            } else if (f2 > f3) {
                // If f2 is greater than f3 then we length alpah2 and alpha3
                // closer to f1 Effectively both are lengthened by a factor of
                // two.
                alpha2 = alpha3;
                f2 = f3;
                alpha3 *= 2;
                for (size_t i = 0; i < xLength; i++)
                    // calculate the new x
                    *param_ptr[i] = xold[i] + alpha3 * s[i];
                f3 = calc(cons, consLength);
                ftimes++;
            }

        // get the alpha for the minimum f of the quadratic approximation
        alphaStar =
            alpha2 + (alpha2 - alpha1) * (f1 - f3) / (3 * (f1 - 2 * f2 + f3));

        // Guarantee that the new alphaStar is within the bracket
        if (alphaStar >= alpha3 || alphaStar <= alpha1)
            alphaStar = alpha2;
        // Avoid NaN
        if (alphaStar != alphaStar)
            alphaStar = 0;

        // Set the values to alphaStar
        for (size_t i = 0; i < xLength; i++)
            // calculate the new x
            *param_ptr[i] = xold[i] + alphaStar * s[i];
        fnew = calc(cons, consLength);
        ftimes++;

        /////////////////////////////////////
        /// end of line search
        ////////////////////////////////////

        deltaXnorm = 0;
        for (size_t i = 0; i < xLength; i++) {
            // Calculate the difference in x for the hessian update
            deltaX[i] = *param_ptr[i] - xold[i];
            deltaXnorm += deltaX[i] * deltaX[i];
            grad[i] = gradnew[i];
        }
        deltaXnorm = sqrt(deltaXnorm);
        iterations++;
        /////////////////////////////////////////////////////////////
        /// End of Main loop
        /////////////////////////////////////////////////////////////
    }

#ifdef DEBUG
    cout << "Fnew: " << fnew << endl;
    cout << "Number of Iterations: " << iterations << endl;
    cout << "Number of function calls: " << ftimes << endl;
#endif

    delete[] grad;
    delete[] s;
    for (size_t i = 0; i < xLength; i++) {
        delete[] N[i];
        delete[] FirstSecond[i];
        delete[] deltaXDotGammatDotN[i];
        delete[] gammatDotDeltaXt[i];
        delete[] NDotGammaDotDeltaXt[i];
    }
    delete[] N;
    delete[] xold;
    delete[] deltaX;
    delete[] gradnew;
    delete[] gamma;
    delete[] gammatDotN;
    delete[] FirstSecond;
    delete[] deltaXDotGammatDotN;
    delete[] gammatDotDeltaXt;
    delete[] NDotGammaDotDeltaXt;

    // End of function
    if (fnew < (isFine ? ValidSolutionFine : ValidSoltuionRough)) {
        delete[] origSolution;
        return Success;
    }
    // Replace the bad numbers with the last result.
    for (size_t i = 0; i < xLength; i++)
        *param_ptr[i] = origSolution[i];
    delete[] origSolution;
    return NoSolution;
}
