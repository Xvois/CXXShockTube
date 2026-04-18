//
// Created by Sonny Parker on 17/04/2026.
//

#ifndef FLUIDSOLVER_CONSTANTS_H
#define FLUIDSOLVER_CONSTANTS_H

const double GAMMA = 1.4;
const int N_ZONES = 100;
const double X_MIN = 0.0;
const double X_MAX = 1.0;
const double DX = (X_MAX - X_MIN) / N_ZONES;
const double CFL_NUMBER = 0.5;

#endif //FLUIDSOLVER_CONSTANTS_H