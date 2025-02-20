#include "common.h"
#include <cmath>
#include <vector>
#include <omp.h>

// Global variables for binning (specific to OpenMP execution)
int bin_count;
double bin_size;
std::vector<std::vector<int>> bins;
std::vector<std::vector<int>> neighbors;

/**
 * Converts a particle's (x, y) position to a bin index.
 */
int get_bin_index(double x, double y) {
    int bx = static_cast<int>(x / bin_size);
    int by = static_cast<int>(y / bin_size);
    return bx * bin_count + by;
}

/**
 * Applies force between two particles (thread-safe).
 */
void apply_force(particle_t& p1, particle_t& p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double r2 = dx * dx + dy * dy;

    if (r2 > cutoff * cutoff || r2 == 0) return; 

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);
    double coef = (1 - cutoff / r) / r2 / mass;

    #pragma omp critical
    {
        p1.ax += coef * dx;
        p1.ay += coef * dy;
        p2.ax -= coef * dx;
        p2.ay -= coef * dy;
    }
}

/**
 * Initializes the simulation with binning.
 */
void init_simulation(particle_t* parts, int num_parts, double size) {
    bin_size = cutoff;
    bin_count = static_cast<int>(size / bin_size) + 1;
    bins.resize(bin_count * bin_count);
    neighbors.resize(bin_count * bin_count);

    #pragma omp parallel for
    for (int i = 0; i < num_parts; i++) {
        int bin_index = get_bin_index(parts[i].x, parts[i].y);
        #pragma omp critical
        bins[bin_index].push_back(i);
    }
}

/**
 * Simulates one step of the particle system (parallel version).
 */
void simulate_one_step(particle_t* parts, int num_parts, double size) {
    #pragma omp parallel for
    for (int i = 0; i < num_parts; i++) {
        parts[i].ax = parts[i].ay = 0;
    }

    #pragma omp parallel for collapse(2)
    for (int bx = 0; bx < bin_count; bx++) {
        for (int by = 0; by < bin_count; by++) {
            int bin_index = bx * bin_count + by;
            for (int i : bins[bin_index]) {
                for (int j : bins[bin_index]) {
                    if (i < j) apply_force(parts[i], parts[j]);
                }
            }
        }
    }
}
