#pragma once

/*
 *  Created on: July 24, 2018
 *  Author: KmolYuan
 */

#include <cmath>
#include "solve.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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


double calc(Constraint *, const int);
