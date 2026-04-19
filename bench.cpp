//
// bench.cpp
//
// Performance benchmark: sequential vs parallel execution.
// Compares wall-clock times for Problem A and Problem B solvers
// using the serial (build_serial) and parallel (build_parallel) builds.
//
// Compile standalone:
//   g++ -std=c++20 -O2 bench.cpp -o bench
//

#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <numeric>
#include <string>
#include <future>
#include <algorithm>

#include "constants.h"
#include "solver.h"

// =========================================================================
// Timing utilities
// =========================================================================

struct BenchResult {
    double mean_ms{0};
    double median_ms{0};
    double min_ms{1e30};
    double max_ms{-1};
    double stddev_ms{0};
    double cv_pct{0};
    int iterations{0};
};

/** Run the benchmark N times and compute statistics. */
BenchResult run_benchmark(void (*func)(const std::string&), int n_iters = 5) {
    std::vector<double> times;
    times.reserve(n_iters);
    
    // Warm-up run
    func("bench_warmup.csv");
    
    // Timed runs
    for (int i = 0; i < n_iters; ++i) {
        auto start = std::chrono::steady_clock::now();
        func(("bench_iter_" + std::to_string(i) + ".csv").c_str());
        auto end = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(end - start).count();
        times.push_back(ms);
    }
    
    // Compute statistics
    BenchResult result;
    result.iterations = n_iters;
    result.mean_ms = std::accumulate(times.begin(), times.end(), 0.0) / n_iters;
    
    double min_val = times[0], max_val = times[0];
    for (double t : times) {
        min_val = std::min(min_val, t);
        max_val = std::max(max_val, t);
    }
    result.min_ms = min_val;
    result.max_ms = max_val;
    
    // Median
    std::vector<double> sorted_times = times;
    std::sort(sorted_times.begin(), sorted_times.end());
    result.median_ms = (n_iters % 2 == 0) 
        ? (sorted_times[n_iters/2 - 1] + sorted_times[n_iters/2]) / 2.0
        : sorted_times[n_iters/2];
    
    // Std dev and CV
    double sq_sum = 0;
    for (double t : times) {
        sq_sum += (t - result.mean_ms) * (t - result.mean_ms);
    }
    result.stddev_ms = std::sqrt(sq_sum / n_iters);
    result.cv_pct = (result.mean_ms > 0) ? (result.stddev_ms / result.mean_ms) * 100.0 : 0;
    
    return result;
}

// =========================================================================
// Main benchmarking logic
// =========================================================================

int main() {
    const int N_ITERS = 5;
    
    std::cout << "============================================\n";
    std::cout << "  Performance Benchmark: Sequential vs Parallel\n";
    std::cout << "============================================\n\n";
    
    std::cout << "Grid: N = " << N_ZONES << ", Gamma = " << GAMMA << "\n";
    std::cout << "Iterations per benchmark: " << N_ITERS << "\n\n";
    
    // Benchmark Problem A and B sequentially
    std::cout << "--- Sequential Mode ---\n";
    BenchResult a_seq = run_benchmark(solveProblemA, N_ITERS);
    BenchResult b_seq = run_benchmark(solveProblemB, N_ITERS);
    std::cout << "\n";
    
    // Benchmark Problem A and B in parallel (using std::async like main.cpp)
    std::cout << "--- Parallel Mode (std::async) ---\n";
    
    // Parallel benchmark for A+B together
    std::vector<double> parallel_times;
    parallel_times.reserve(N_ITERS);
    
    // Warm-up
    std::future<void> futA = std::async(std::launch::async, solveProblemA, "bench_warmup.csv");
    std::future<void> futB = std::async(std::launch::async, solveProblemB, "bench_warmup.csv");
    futA.get();
    futB.get();
    
    // Timed runs
    for (int i = 0; i < N_ITERS; ++i) {
        auto start = std::chrono::steady_clock::now();
        futA = std::async(std::launch::async, solveProblemA, ("bench_iter_A_" + std::to_string(i) + ".csv").c_str());
        futB = std::async(std::launch::async, solveProblemB, ("bench_iter_B_" + std::to_string(i) + ".csv").c_str());
        futA.get();
        futB.get();
        auto end = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(end - start).count();
        parallel_times.push_back(ms);
    }
    
    BenchResult par_ab;
    par_ab.iterations = N_ITERS;
    par_ab.mean_ms = std::accumulate(parallel_times.begin(), parallel_times.end(), 0.0) / N_ITERS;
    double min_val = parallel_times[0], max_val = parallel_times[0];
    for (double t : parallel_times) {
        min_val = std::min(min_val, t);
        max_val = std::max(max_val, t);
    }
    par_ab.min_ms = min_val;
    par_ab.max_ms = max_val;
    
    std::vector<double> sorted_pt = parallel_times;
    std::sort(sorted_pt.begin(), sorted_pt.end());
    par_ab.median_ms = (N_ITERS % 2 == 0) 
        ? (sorted_pt[N_ITERS/2 - 1] + sorted_pt[N_ITERS/2]) / 2.0
        : sorted_pt[N_ITERS/2];
    
    double sq_sum = 0;
    for (double t : parallel_times) {
        sq_sum += (t - par_ab.mean_ms) * (t - par_ab.mean_ms);
    }
    par_ab.stddev_ms = std::sqrt(sq_sum / N_ITERS);
    par_ab.cv_pct = (par_ab.mean_ms > 0) ? (par_ab.stddev_ms / par_ab.mean_ms) * 100.0 : 0;
    
    // Print results
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Sequential:\n";
    std::cout << "  Problem A: mean=" << a_seq.mean_ms << "ms, median=" << a_seq.median_ms 
              << "ms, min=" << a_seq.min_ms << "ms, max=" << a_seq.max_ms << "ms, cv=" 
              << a_seq.cv_pct << "%\n";
    std::cout << "  Problem B: mean=" << b_seq.mean_ms << "ms, median=" << b_seq.median_ms 
              << "ms, min=" << b_seq.min_ms << "ms, max=" << b_seq.max_ms << "ms, cv=" 
              << b_seq.cv_pct << "%\n\n";
    
    std::cout << "Parallel (A+B together):\n";
    std::cout << "  Combined: mean=" << par_ab.mean_ms << "ms, median=" << par_ab.median_ms 
              << "ms, min=" << par_ab.min_ms << "ms, max=" << par_ab.max_ms << "ms, cv=" 
              << par_ab.cv_pct << "%\n\n";
    
    // Calculate speedup
    double total_seq = a_seq.mean_ms + b_seq.mean_ms;
    double speedup_combined = (par_ab.mean_ms > 0) ? total_seq / par_ab.mean_ms : 0;
    
    std::cout << "Summary:\n";
    std::cout << "  Sequential total (A+B): " << total_seq << "ms\n";
    std::cout << "  Parallel combined (A+B): " << par_ab.mean_ms << "ms\n";
    std::cout << "  Speedup (combined): " << std::fixed << std::setprecision(2) << speedup_combined << "x\n\n";
    
    // Calculate theoretical individual parallel times
    double a_par = a_seq.mean_ms;  // Single problem doesn't benefit from async
    double b_par = b_seq.mean_ms;
    double theoretical_par = std::max(a_seq.mean_ms, b_seq.mean_ms);
    double theoretical_speedup = (theoretical_par > 0) ? total_seq / theoretical_par : 0;
    
    std::cout << "Theoretical max speedup (single problem): " << theoretical_speedup << "x\n";
    std::cout << "Actual combined speedup: " << speedup_combined << "x\n";
    std::cout << "Efficiency: " << std::fixed << std::setprecision(1) 
              << (speedup_combined / theoretical_speedup * 100.0) << "%\n";
    
    // Write to file
    std::ofstream out("performance_metrics.txt");
    out << "============================================\n";
    out << "  Performance Metrics: Sequential vs Parallel\n";
    out << "============================================\n\n";
    out << "Configuration:\n";
    out << "  Grid size: N = " << N_ZONES << "\n";
    out << "  Gamma: " << GAMMA << "\n";
    out << "  Benchmark iterations: " << N_ITERS << "\n\n";
    
    out << "Sequential Mode:\n";
    out << "  Problem A: mean=" << std::fixed << std::setprecision(2) << a_seq.mean_ms 
        << "ms, median=" << a_seq.median_ms << "ms, min=" << a_seq.min_ms 
        << "ms, max=" << a_seq.max_ms << "ms, stddev=" << a_seq.stddev_ms 
        << "ms, cv=" << a_seq.cv_pct << "%\n";
    out << "  Problem B: mean=" << b_seq.mean_ms << "ms, median=" << b_seq.median_ms 
        << "ms, min=" << b_seq.min_ms << "ms, max=" << b_seq.max_ms 
        << "ms, stddev=" << b_seq.stddev_ms << "ms, cv=" << b_seq.cv_pct << "%\n\n";
    
    out << "Parallel Mode (std::async, A+B together):\n";
    out << "  Combined: mean=" << par_ab.mean_ms << "ms, median=" << par_ab.median_ms 
        << "ms, min=" << par_ab.min_ms << "ms, max=" << par_ab.max_ms 
        << "ms, stddev=" << par_ab.stddev_ms << "ms, cv=" << par_ab.cv_pct << "%\n\n";
    
    out << "Summary:\n";
    out << "  Sequential total (A+B): " << total_seq << "ms\n";
    out << "  Parallel combined (A+B): " << par_ab.mean_ms << "ms\n";
    out << "  Speedup: " << speedup_combined << "x\n";
    out << "  Theoretical max speedup: " << theoretical_speedup << "x\n";
    out << "  Efficiency: " << (speedup_combined / theoretical_speedup * 100.0) << "%\n\n";
    
    out << "Key Observations:\n";
    out << "  - Parallel execution overlaps A and B computation\n";
    out << "  - Speedup is limited by the slower problem\n";
    out << "  - Thread overhead reduces actual speedup below theoretical\n";
    out << "  - CV < 5% indicates consistent timing across iterations\n";
    out.close();
    
    std::cout << "\nPerformance metrics written to performance_metrics.txt\n";
    
    return 0;
}
