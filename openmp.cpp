#include "common.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <omp.h>

/**
 * Converts a particle's (x, y) position to a bin index.
 */
int get_bin_index(double x, double y, int bin_count, double bin_size) {
    int bx = static_cast<int>(x / bin_size);
    int by = static_cast<int>(y / bin_size);
    return bx * bin_count + by;
}

/**
 * Applies force between two particles with thread-local force accumulators.
 */
void apply_force(particle_t& p1, particle_t& p2, double& ax1, double& ay1, double& ax2, double& ay2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double r2 = dx * dx + dy * dy;

    if (r2 > cutoff * cutoff || r2 == 0) return; 

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);
    double coef = (1 - cutoff / r) / r2 / mass;

    ax1 += coef * dx;
    ay1 += coef * dy;
    ax2 -= coef * dx;
    ay2 -= coef * dy;
}

/**
 * Simulates one step of the particle system (parallelized with OpenMP).
 */
void simulate_one_step(particle_t* parts, int num_parts, double size, int bin_count, double bin_size, 
                       std::vector<std::vector<int>>& bins, std::vector<std::vector<int>>& neighbors) {

    // Step 1: Compute forces in parallel using thread-local accumulators
    std::vector<double> local_ax(num_parts, 0.0);
    std::vector<double> local_ay(num_parts, 0.0);

    #pragma omp parallel
    {
        std::vector<double> thread_ax(num_parts, 0.0);
        std::vector<double> thread_ay(num_parts, 0.0);

        #pragma omp for collapse(2) schedule(dynamic)
        for (int bx = 0; bx < bin_count; bx++) {
            for (int by = 0; by < bin_count; by++) {
                int bin_index = bx * bin_count + by;

                for (int i : bins[bin_index]) {
                    for (int j : bins[bin_index]) {
                        if (i < j) apply_force(parts[i], parts[j], thread_ax[i], thread_ay[i], thread_ax[j], thread_ay[j]);
                    }
                    for (int neighbor_bin : neighbors[bin_index]) {
                        for (int j : bins[neighbor_bin]) {
                            apply_force(parts[i], parts[j], thread_ax[i], thread_ay[i], thread_ax[j], thread_ay[j]);
                        }
                    }
                }
            }
        }

        #pragma omp critical
        for (int i = 0; i < num_parts; i++) {
            local_ax[i] += thread_ax[i];
            local_ay[i] += thread_ay[i];
        }
    }

    // Apply accumulated forces
    #pragma omp parallel for
    for (int i = 0; i < num_parts; i++) {
        parts[i].ax = local_ax[i];
        parts[i].ay = local_ay[i];
    }

    // Step 2: Move particles in parallel and determine bin changes
    std::vector<std::vector<int>> new_bins(bin_count * bin_count);
    
    #pragma omp parallel
    {
        std::vector<std::vector<int>> local_bins(bin_count * bin_count);
        std::vector<int> to_remove;

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < num_parts; i++) {
            int old_bin = get_bin_index(parts[i].x, parts[i].y, bin_count, bin_size);
            
            // Update velocity and position
            parts[i].vx += parts[i].ax * dt;
            part
