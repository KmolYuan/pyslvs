/*
 * solve.cpp
 *
 *  Created on: May 4, 2009
 *  Author: Jonathan George
 *  Contributor: KmolYuan
 */

#include "solve.h"
#include <cmath>
#include <vector>
#include "calc.h"

typedef std::vector<double> array1d;
typedef std::vector<array1d> array2d;

double cal_grad(double *param, Constraint *cons, size_t cons_len, double pert) {
    double tmp = *param;
    *param = tmp - pert;
    double first = calc(cons, cons_len);
    *param = tmp + pert;
    double second = calc(cons, cons_len);
    *param = tmp;
    return 0.5 * (second - first) / pert;
}

int solve(double **param, size_t param_len, Constraint *cons, size_t cons_len,
          bool is_fine) {
    // Save the original parameters for later
    array1d param_origin(param_len);
    for (size_t i = 0; i < param_len; i++)
        param_origin[i] = *param[i];
    // Calculate Function at the starting point:
    double f0 = calc(cons, cons_len);
    if (f0 < SmallF)
        return Success;

    // Calculate the gradient at the starting point:
    // The gradient vector (1xn)
    array1d grad(param_len);
    double pert = f0 * PertMag;
    // The norm of the gradient vector
    for (size_t j = 0; j < param_len; j++)
        grad[j] = cal_grad(param[j], cons, cons_len, pert);

    // Estimate the norm of n
    // Initialize n and calculate s
    // The current search direction
    array1d s(param_len);
    // The estimate of the Hessian inverse
    array2d n(param_len, array1d(param_len));
    for (size_t i = 0; i < param_len; i++)
        for (size_t j = 0; j < param_len; j++) {
            if (i == j) {
                n[i][j] = 1;
                // Calculate the initial search vector
                s[i] = -grad[i];
                continue;
            }
            n[i][j] = 0;
        }

    // Initial search vector multiplier
    double alpha = 1;
    // Storage for the previous design variables
    array1d x_old(param_len);
    for (size_t i = 0; i < param_len; i++)
        // Copy last values to x_old
        x_old[i] = *param[i];

    ///////////////////////////////////////////////////////
    /// Start of line search
    ///////////////////////////////////////////////////////
    // Make the initial position alpha1
    double alpha1 = 0;
    double f1 = f0;
    // Take a step of alpha=1 as alpha2
    double alpha2 = 1;
    for (size_t i = 0; i < param_len; i++)
        // calculate the new x
        *param[i] = x_old[i] + alpha2 * s[i];
    double f2 = calc(cons, cons_len);

    // Take a step of alpha 3 that is 2*alpha2
    double alpha3 = alpha * 2;
    for (size_t i = 0; i < param_len; i++)
        // calculate the new x
        *param[i] = x_old[i] + alpha3 * s[i];
    double f3 = calc(cons, cons_len);

    // Now reduce or lengthen alpha2 and alpha3 until the minimum is
    // Bracketed by the triplet f1>f2<f3
    while (f2 > f1 || f2 > f3)
        if (f2 > f1) {
            // If f2 is greater than f1 then we shorten alpha2 and alpha3 closer
            // to f1 Effectively both are shortened by a factor of two
            alpha3 = alpha2;
            f3 = f2;
            alpha2 *= 0.5;
            for (size_t i = 0; i < param_len; i++)
                // calculate the new x
                *param[i] = x_old[i] + alpha2 * s[i];
            f2 = calc(cons, cons_len);
        } else if (f2 > f3) {
            // If f2 is greater than f3 then we length alpah2 and alpha3 closer
            // to f1 Effectively both are lengthened by a factor of two
            alpha2 = alpha3;
            f2 = f3;
            alpha3 *= 2;
            for (size_t i = 0; i < param_len; i++)
                // calculate the new x
                *param[i] = x_old[i] + alpha3 * s[i];
            f3 = calc(cons, cons_len);
        }
    // get the alpha for the minimum f of the quadratic approximation
    double alpha_s =
        alpha2 + (alpha2 - alpha1) * (f1 - f3) / (3 * (f1 - 2 * f2 + f3));

    // Guarantee that the new alpha_s is within the bracket
    if (alpha_s > alpha3 || alpha_s < alpha1)
        alpha_s = alpha2;
    if (alpha_s != alpha_s)
        // Fix nan problem
        alpha_s = 0.001;

    /// Set the values to alpha_s
    for (size_t i = 0; i < param_len; i++)
        // calculate the new x
        *param[i] = x_old[i] + alpha_s * s[i];
    double fnew = calc(cons, cons_len);

    /////////////////////////////////////
    /// end of line search
    /////////////////////////////////////
    array1d delta_x(param_len);
    array1d grad_new(param_len);
    array1d gamma(param_len);
    array1d gn(param_len);  // gamma' dot n
    array2d first_second(param_len, array1d(param_len));
    array2d dgn(param_len, array1d(param_len));  // delta_x dot gamma' dot n
    array2d gd(param_len, array1d(param_len));   // gamma' dot delta_x
    array2d ngd(param_len, array1d(param_len));  // n dot gamma dot delta_x
    // Calculate delta_x
    for (size_t i = 0; i < param_len; i++)
        // Calculate the difference in x for the Hessian update
        delta_x[i] = *param[i] - x_old[i];

    size_t iterations = 1;
    double delta_x_norm = 1;
    while (delta_x_norm > (is_fine ? XConvergenceFine : XConvergenceRough) &&
           fnew > SmallF && iterations < MaxIterations * param_len) {
        //////////////////////////////////////////////////////////////////////
        /// Start of main loop!!!!
        //////////////////////////////////////////////////////////////////////
        double bottom = 0;
        pert = fnew * PertMag;
        if (pert < PertMin)
            pert = PertMin;
        for (size_t i = 0; i < param_len; i++) {
            grad_new[i] = cal_grad(param[i], cons, cons_len, pert);
            // Calculate the change in the gradient
            gamma[i] = grad_new[i] - grad[i];
            bottom += delta_x[i] * gamma[i];
        }
        // make sure that bottom is never 0
        if (bottom == 0.)
            bottom = 1e-10;
        // calculate all (1xn)dot(nxn)
        for (size_t i = 0; i < param_len; i++) {
            gn[i] = 0;
            for (size_t j = 0; j < param_len; j++)
                // This is gn transpose
                gn[i] += gamma[j] * n[i][j];
        }
        // calculate all (1xn)dot(nx1)
        double gng = 0;  // gamma' dot n dot gamma
        for (size_t i = 0; i < param_len; i++)
            gng += gn[i] * gamma[i];
        // Calculate the first term
        double first_term = 1 + gng / bottom;
        // Calculate all (nx1)dot(1xn) matrices
        for (size_t i = 0; i < param_len; i++)
            for (size_t j = 0; j < param_len; j++) {
                first_second[i][j] =
                    delta_x[j] * delta_x[i] / bottom * first_term;
                dgn[i][j] = delta_x[i] * gn[j];
                gd[i][j] = gamma[i] * delta_x[j];
            }
        // Calculate all (nxn)dot(nxn) matrices
        for (size_t i = 0; i < param_len; i++)
            for (size_t j = 0; j < param_len; j++) {
                ngd[i][j] = 0;
                for (size_t k = 0; k < param_len; k++)
                    ngd[i][j] += n[i][k] * gd[k][j];
            }
        // Now calculate the BFGS update on n
        for (size_t i = 0; i < param_len; i++)
            for (size_t j = 0; j < param_len; j++)
                n[i][j] +=
                    first_second[i][j] - (dgn[i][j] + ngd[i][j]) / bottom;
        // Calculates
        for (size_t i = 0; i < param_len; i++) {
            s[i] = 0;
            for (size_t j = 0; j < param_len; j++)
                s[i] += -n[i][j] * grad_new[j];
        }
        // copy newest values to the x_old
        for (size_t i = 0; i < param_len; i++)
            x_old[i] = *param[i];  // Copy last values to x_old

        ///////////////////////////////////////////////////////
        /// Start of line search
        ///////////////////////////////////////////////////////
        // Make the initial position alpha1
        alpha1 = 0;
        f1 = fnew;

        // Take a step of alpha=1 as alpha2
        alpha2 = 1;
        for (size_t i = 0; i < param_len; i++)
            // calculate the new x
            *param[i] = x_old[i] + alpha2 * s[i];
        f2 = calc(cons, cons_len);

        // Take a step of alpha 3 that is 2*alpha2
        alpha3 = alpha2 * 2;
        for (size_t i = 0; i < param_len; i++)
            // calculate the new x
            *param[i] = x_old[i] + alpha3 * s[i];
        f3 = calc(cons, cons_len);

        // Now reduce or lengthen alpha2 and alpha3 until the minimum is
        // Bracketed by the triplet f1>f2<f3
        while (f2 > f1 || f2 > f3)
            if (f2 > f1) {
                // If f2 is greater than f1 then we shorten alpha2 and alpha3
                // closer to f1 Effectively both are shortened by a factor of
                // two
                alpha3 = alpha2;
                f3 = f2;
                alpha2 /= 2;
                for (size_t i = 0; i < param_len; i++)
                    // calculate the new x
                    *param[i] = x_old[i] + alpha2 * s[i];
                f2 = calc(cons, cons_len);
            } else if (f2 > f3) {
                // If f2 is greater than f3 then we length alpah2 and alpha3
                // closer to f1 Effectively both are lengthened by a factor of
                // two
                alpha2 = alpha3;
                f2 = f3;
                alpha3 *= 2;
                for (size_t i = 0; i < param_len; i++)
                    // calculate the new x
                    *param[i] = x_old[i] + alpha3 * s[i];
                f3 = calc(cons, cons_len);
            }

        // get the alpha for the minimum f of the quadratic approximation
        alpha_s =
            alpha2 + (alpha2 - alpha1) * (f1 - f3) / (3 * (f1 - 2 * f2 + f3));

        // Guarantee that the new alpha_s is within the bracket
        if (alpha_s >= alpha3 || alpha_s <= alpha1)
            alpha_s = alpha2;
        // Avoid NaN
        if (alpha_s != alpha_s)
            alpha_s = 0;

        // Set the values to alpha_s
        for (size_t i = 0; i < param_len; i++)
            // calculate the new x
            *param[i] = x_old[i] + alpha_s * s[i];
        fnew = calc(cons, cons_len);

        /////////////////////////////////////
        /// end of line search
        ////////////////////////////////////
        delta_x_norm = 0;
        for (size_t i = 0; i < param_len; i++) {
            // Calculate the difference in x for the hessian update
            delta_x[i] = *param[i] - x_old[i];
            delta_x_norm += delta_x[i] * delta_x[i];
            grad[i] = grad_new[i];
        }
        delta_x_norm = sqrt(delta_x_norm);
        iterations++;
        /////////////////////////////////////////////////////////////
        /// End of Main loop
        /////////////////////////////////////////////////////////////
    }
    // End of function
    if (fnew < (is_fine ? ValidSolutionFine : ValidSoltuionRough))
        return Success;
    // Replace the bad numbers with the last result
    for (size_t i = 0; i < param_len; i++)
        *param[i] = param_origin[i];
    return NoSolution;
}
