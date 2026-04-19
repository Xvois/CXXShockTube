//
// constants.h
//
// Global simulation parameters for the Euler equation solver.
//
#ifndef FLUIDSOLVER_CONSTANTS_H
#define FLUIDSOLVER_CONSTANTS_H

const double GAMMA  = 1.4;    // Adiabatic exponent
const int    N_ZONES = 100;   // Number of computational zones
const double X_MIN   = 0.0;   // Domain lower bound
const double X_MAX   = 1.0;   // Domain upper bound
const double DX      = (X_MAX - X_MIN) / N_ZONES;  // Zone width
const double CFL_NUMBER = 0.5;  // CFL stability number

#endif // FLUIDSOLVER_CONSTANTS_H
