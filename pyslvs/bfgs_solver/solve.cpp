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

#define EPS 1e-20
#define MAX_ITER 50
#define PERT_MAG 1e-6
#define PERT_MIN 1e-10
#define CONVERGENCE_ROUGH 1e-8
#define CONVERGENCE_FINE 1e-10
#define VALID_SOLUTION_FINE 1e-12
#define VALID_SOLUTION_ROUGH 1e-4

// Dynamic arrays
using array1d = std::vector<double>;
using array2d = std::vector<array1d>;

namespace {
auto cal_grad(double *param, Constraint *cons, size_t cons_len, double pert)
    -> double {
    auto tmp = *param;
    *param = tmp - pert;
    auto first = calc(cons, cons_len);
    *param = tmp + pert;
    auto second = calc(cons, cons_len);
    *param = tmp;
    return 0.5 * (second - first) / pert;
}
}  // namespace

auto solve(double **param, size_t param_len, Constraint *cons, size_t cons_len,
           bool is_fine) -> bool {
    // Save the original parameters for later
    auto param_origin = array1d(param_len);
    for (auto i = 0u; i < param_len; i++)
        param_origin[i] = *param[i];
    // Calculate Function at the starting point:
    auto f0 = calc(cons, cons_len);
    if (f0 < EPS)
        return true;

    // Calculate the gradient at the starting point:
    // The gradient vector (1xn)
    auto grad = array1d(param_len);
    auto pert = f0 * PERT_MAG;
    // The norm of the gradient vector
    for (auto j = 0u; j < param_len; j++)
        grad[j] = cal_grad(param[j], cons, cons_len, pert);

    // Estimate the norm of n
    // Initialize n and calculate s
    // The current search direction
    auto s = array1d(param_len);
    // The estimate of the Hessian inverse
    auto n = array2d(param_len, array1d(param_len));
    for (auto i = 0u; i < param_len; i++)
        for (auto j = 0u; j < param_len; j++) {
            if (i == j) {
                n[i][j] = 1;
                // Calculate the initial search vector
                s[i] = -grad[i];
                continue;
            }
            n[i][j] = 0;
        }

    // Initial search vector multiplier
    auto alpha = 1.;
    // Storage for the previous design variables
    auto x_old = array1d(param_len);
    for (auto i = 0u; i < param_len; i++)
        // Copy last values to x_old
        x_old[i] = *param[i];

    ///////////////////////////////////////////////////////
    /// Start of line search
    ///////////////////////////////////////////////////////
    // Make the initial position alpha1
    auto alpha1 = 0.;
    auto f1 = f0;
    // Take a step of alpha=1 as alpha2
    auto alpha2 = 1.;
    for (auto i = 0u; i < param_len; i++)
        // calculate the new x
        *param[i] = x_old[i] + alpha2 * s[i];
    auto f2 = calc(cons, cons_len);

    // Take a step of alpha 3 that is 2*alpha2
    auto alpha3 = alpha * 2;
    for (auto i = 0u; i < param_len; i++)
        // calculate the new x
        *param[i] = x_old[i] + alpha3 * s[i];
    auto f3 = calc(cons, cons_len);

    // Now reduce or lengthen alpha2 and alpha3 until the minimum is
    // Bracketed by the triplet f1>f2<f3
    while (f2 > f1 || f2 > f3)
        if (f2 > f1) {
            // If f2 is greater than f1 then we shorten alpha2 and alpha3 closer
            // to f1 Effectively both are shortened by a factor of two
            alpha3 = alpha2;
            f3 = f2;
            alpha2 *= 0.5;
            for (auto i = 0u; i < param_len; i++)
                // calculate the new x
                *param[i] = x_old[i] + alpha2 * s[i];
            f2 = calc(cons, cons_len);
        } else if (f2 > f3) {
            // If f2 is greater than f3 then we length alpah2 and alpha3 closer
            // to f1 Effectively both are lengthened by a factor of two
            alpha2 = alpha3;
            f2 = f3;
            alpha3 *= 2;
            for (auto i = 0u; i < param_len; i++)
                // calculate the new x
                *param[i] = x_old[i] + alpha3 * s[i];
            f3 = calc(cons, cons_len);
        }
    // get the alpha for the minimum f of the quadratic approximation
    auto alpha_s =
        alpha2 + (alpha2 - alpha1) * (f1 - f3) / (3 * (f1 - 2 * f2 + f3));

    // Guarantee that the new alpha_s is within the bracket
    if (alpha_s > alpha3 || alpha_s < alpha1)
        alpha_s = alpha2;
    if (std::isnan(alpha_s))
        // Fix nan problem
        alpha_s = 0.001;

    /// Set the values to alpha_s
    for (auto i = 0u; i < param_len; i++)
        // calculate the new x
        *param[i] = x_old[i] + alpha_s * s[i];
    auto fnew = calc(cons, cons_len);

    /////////////////////////////////////
    /// end of line search
    /////////////////////////////////////
    auto delta_x = array1d(param_len);
    auto grad_new = array1d(param_len);
    auto gamma = array1d(param_len);
    auto gn = array1d(param_len);  // gamma' @ n
    auto first_second = array2d(param_len, array1d(param_len));
    // delta_x @ gamma' @ n
    auto dgn = array2d(param_len, array1d(param_len));
    // gamma' @ delta_x
    auto gd = array2d(param_len, array1d(param_len));
    // n @ gamma @ delta_x
    auto ngd = array2d(param_len, array1d(param_len));
    // Calculate delta_x
    for (auto i = 0u; i < param_len; i++)
        // Calculate the difference in x for the Hessian update
        delta_x[i] = *param[i] - x_old[i];

    auto iterations = 1u;
    auto delta_x_norm = 1.;
    while (delta_x_norm > (is_fine ? CONVERGENCE_FINE : CONVERGENCE_ROUGH)
           && fnew > EPS && iterations < MAX_ITER * param_len) {
        //////////////////////////////////////////////////////////////////////
        /// Start of main loop!!!!
        //////////////////////////////////////////////////////////////////////
        auto bottom = 0.;
        pert = fnew * PERT_MAG;
        if (pert < PERT_MIN)
            pert = PERT_MIN;
        for (auto i = 0u; i < param_len; i++) {
            grad_new[i] = cal_grad(param[i], cons, cons_len, pert);
            // Calculate the change in the gradient
            gamma[i] = grad_new[i] - grad[i];
            bottom += delta_x[i] * gamma[i];
        }
        // make sure that bottom is never 0
        if (bottom == 0.)
            bottom = 1e-10;
        // calculate all (1xn)@(nxn)
        for (auto i = 0u; i < param_len; i++) {
            gn[i] = 0;
            for (auto j = 0u; j < param_len; j++)
                // This is gn transpose
                gn[i] += gamma[j] * n[i][j];
        }
        // calculate all (1xn)@(nx1)
        auto gng = 0.;  // gamma' @ n @ gamma
        for (auto i = 0u; i < param_len; i++)
            gng += gn[i] * gamma[i];
        // Calculate the first term
        auto first_term = 1 + gng / bottom;
        // Calculate all (nx1)@(1xn) matrices
        for (auto i = 0u; i < param_len; i++)
            for (auto j = 0u; j < param_len; j++) {
                first_second[i][j] =
                    delta_x[j] * delta_x[i] / bottom * first_term;
                dgn[i][j] = delta_x[i] * gn[j];
                gd[i][j] = gamma[i] * delta_x[j];
            }
        // Calculate all (nxn)@(nxn) matrices
        for (auto i = 0u; i < param_len; i++)
            for (auto j = 0u; j < param_len; j++) {
                ngd[i][j] = 0;
                for (auto k = 0u; k < param_len; k++)
                    ngd[i][j] += n[i][k] * gd[k][j];
            }
        // Now calculate the BFGS update on n
        for (auto i = 0u; i < param_len; i++)
            for (auto j = 0u; j < param_len; j++)
                n[i][j] +=
                    first_second[i][j] - (dgn[i][j] + ngd[i][j]) / bottom;
        // Calculates
        for (auto i = 0u; i < param_len; i++) {
            s[i] = 0;
            for (auto j = 0u; j < param_len; j++)
                s[i] += -n[i][j] * grad_new[j];
        }
        // copy newest values to the x_old
        for (auto i = 0u; i < param_len; i++)
            x_old[i] = *param[i];  // Copy last values to x_old

        ///////////////////////////////////////////////////////
        /// Start of line search
        ///////////////////////////////////////////////////////
        // Make the initial position alpha1
        alpha1 = 0;
        f1 = fnew;

        // Take a step of alpha=1 as alpha2
        alpha2 = 1;
        for (auto i = 0u; i < param_len; i++)
            // calculate the new x
            *param[i] = x_old[i] + alpha2 * s[i];
        f2 = calc(cons, cons_len);

        // Take a step of alpha 3 that is 2*alpha2
        alpha3 = alpha2 * 2;
        for (auto i = 0u; i < param_len; i++)
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
                alpha2 *= 0.5;
                for (auto i = 0u; i < param_len; i++)
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
                for (auto i = 0u; i < param_len; i++)
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

        if (std::isnan(alpha_s))
            alpha_s = 0;

        // Set the values to alpha_s
        for (auto i = 0u; i < param_len; i++)
            // calculate the new x
            *param[i] = x_old[i] + alpha_s * s[i];
        fnew = calc(cons, cons_len);

        /////////////////////////////////////
        /// end of line search
        ////////////////////////////////////
        delta_x_norm = 0;
        for (auto i = 0u; i < param_len; i++) {
            // Calculate the difference in x for the hessian update
            delta_x[i] = *param[i] - x_old[i];
            delta_x_norm += delta_x[i] * delta_x[i];
            grad[i] = grad_new[i];
        }
        delta_x_norm = sqrt(delta_x_norm);
        iterations += 1;
        /////////////////////////////////////////////////////////////
        /// End of Main loop
        /////////////////////////////////////////////////////////////
    }
    // End of function
    if (fnew < (is_fine ? VALID_SOLUTION_FINE : VALID_SOLUTION_ROUGH))
        return true;
    // Replace the bad numbers with the last result
    for (auto i = 0u; i < param_len; i++)
        *param[i] = param_origin[i];
    return false;
}
