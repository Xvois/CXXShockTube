//
// Created by Sonny Parker on 17/04/2026.
//

#ifndef FLUIDSOLVER_ANALYTICAL_H
#define FLUIDSOLVER_ANALYTICAL_H
#include <cmath>

#include "constants.h"
#include "physics.h"

double calculateSoundSpeed(const Primitive& w);
double pressureFunction(double p_star, double p_side, double rho_side, double cs_side);
double findStarPressure(const Primitive& L, const Primitive& R);
Primitive sampleExactSolution(const Primitive& L, const Primitive& R, double x_initial, double x_current, double t);

#endif //FLUIDSOLVER_ANALYTICAL_H