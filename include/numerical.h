//
// Created by Sonny Parker on 17/04/2026.
//

#ifndef FLUIDSOLVER_NUMERICAL_H
#define FLUIDSOLVER_NUMERICAL_H
#include <vector>
#include "physics.h"
#include "constants.h"

double calculateTimeStep(const std::vector<Conserved>& grid);
std::vector<Conserved> updateLaxFriedrichs(const std::vector<Conserved>& grid, double dt);

#endif //FLUIDSOLVER_NUMERICAL_H