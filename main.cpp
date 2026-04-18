#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>

#include "include/constants.h"
#include "include/physics.h"
#include "include/analytical.h"
#include "include/numerical.h"



// --- Initialization ---

std::vector<Conserved> initialiseProblemA() {
    std::vector<Conserved> grid(N_ZONES);
    for (int z = 0; z < N_ZONES; ++z) {
        double x = X_MIN + (z + 0.5) * DX;
        Primitive w{};
        if (x < 0.3) {
            w.rho = 1.0;
            w.p = 1.0;
            w.v = 0.75;
        } else {
            w.rho = 0.125;
            w.p = 0.1;
            w.v = 0.0;
        }
        grid[z] = primToCons(w);
    }
    return grid;
}


// --- Output ---

void outputCSV(const std::vector<Conserved>& grid, double currentTime, const std::string& filename) {
    std::ofstream file;
    if (currentTime == 0.0) {
        file.open(filename, std::ios::out); // Overwrite/Create new
        file << "time,x,density,velocity,pressure\n";
    } else {
        file.open(filename, std::ios::app); // Append
    }

    for (int i = 0; i < N_ZONES; ++i) {
        double x = X_MIN + (i + 0.5) * DX;
        Primitive w = consToPrim(grid[i]);
        file << std::fixed << std::setprecision(6)
             << currentTime << "," << x << "," << w.rho << "," << w.v << "," << w.p << "\n";
    }
    file.close();
}

// --- Main ---

int main() {
    std::vector<Conserved> grid = initialiseProblemA();

    double t = 0.0;
    const double T_FINAL = 0.2;
    const double OUTPUT_INTERVAL = 0.05;
    double next_output_t = 0.0;
    int step = 0;

    std::string numericalFN = "12345_problemA_results.csv";
    std::string analyticalFN = "exact_solution.csv";

    std::cout << "Starting simulation..." << std::endl;

    while (t <= T_FINAL) {
        // Output snapshot
        if (t >= next_output_t - 1e-10) {
            outputCSV(grid, t, numericalFN);
            next_output_t += OUTPUT_INTERVAL;
        }

        if (t >= T_FINAL) break;

        double dt = calculateTimeStep(grid);

        // Clip dt to hit output intervals exactly
        if (t + dt > next_output_t && next_output_t <= T_FINAL) dt = next_output_t - t;
        if (t + dt > T_FINAL) dt = T_FINAL - t;

        grid = updateLaxFriedrichs(grid, dt);
        t += dt;
        step++;
    }

    std::cout << "Simulation finished. Data saved to " << numericalFN << std::endl;

    std::vector<Conserved> exact_grid(N_ZONES);
    Primitive L = {1.0, 0.75, 1.0}; // Values for Problem A
    Primitive R = {0.125, 0.0, 0.1};

    for(int i=0; i < N_ZONES; ++i) {
        double x = X_MIN + (i + 0.5) * DX;
        Primitive w_exact = sampleExactSolution(L, R, 0.3, x, 0.2);
        exact_grid[i] = primToCons(w_exact);
    }
    outputCSV(exact_grid, 0.2, analyticalFN);

    std::cout << "Simulation finished. Data saved to " << analyticalFN << std::endl;
    return 0;
}