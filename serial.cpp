#include "common.h"
#include <cmath>
#include <vector>
#include <iostream>

// Global variables
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
 * Applies force between two particles.
 */
void apply_force(particle_t& p1, particle_t& p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double r2 = dx * dx + dy * dy;

    if (r2 > cutoff * cutoff || r2 == 0) return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);
    double coef = (1 - cutoff / r) / r2 / mass;

    p1.ax += coef * dx;
    p1.ay += coef * dy;
    p2.ax -= coef * dx;
    p2.ay -= coef * dy;
}

/**
 * Initializes the simulation with binning.
 */
void init_simulation(particle_t* parts, int num_parts, double size) {
    std::cout << "Init Simulation (Serial): num_parts = " << num_parts 
              << ", size = " << size << std::endl;

    bin_size = cutoff;
    bin_count = static_cast<int>(size / bin_size) + 1;
    bins.resize(bin_count * bin_count);
    neighbors.resize(bin_count * bin_count);

    // Precompute valid neighbor bins
    for (int bx = 0; bx < bin_count; bx++) {
        for (int by = 0; by < bin_count; by++) {
            int bin_index = bx * bin_count + by;
            if (bx > 0) neighbors[bin_index].push_back((bx - 1) * bin_count + by);
            if (by < bin_count - 1) neighbors[bin_index].push_back(bx * bin_count + (by + 1));
        }
    }

    // Assign particles to bins
    for (int i = 0; i < num_parts; i++) {
        int bin_index = get_bin_index(parts[i].x, parts[i].y);
        bins[bin_index].push_back(i);
    }
}

/**
 * Simulates one step of the particle system (serial version).
 */
void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // Reset forces
    for (int i = 0; i < num_parts; i++) {
        parts[i].ax = parts[i].ay = 0;
    }

    // Compute forces using bins
    for (int bx = 0; bx < bin_count; bx++) {
        for (int by = 0; by < bin_count; by++) {
            int bin_index = bx * bin_count + by;
            for (int i : bins[bin_index]) {
                for (int j : bins[bin_index]) {
                    if (i < j) apply_force(parts[i], parts[j]);
                }
                for (int neighbor_bin : neighbors[bin_index]) {
                    for (int j : bins[neighbor_bin]) {
                        apply_force(parts[i], parts[j]);
                    }
                }
            }
        }
    }

    // Move particles and update bins
    std::vector<std::vector<int>> new_bins(bin_count * bin_count);
    for (int i = 0; i < num_parts; i++) {
        int old_bin = get_bin_index(parts[i].x, parts[i].y);

        // Update velocity and position
        parts[i].vx += parts[i].ax * dt;
        parts[i].vy += parts[i].ay * dt;
        parts[i].x += parts[i].vx * dt;
        parts[i].y += parts[i].vy * dt;

        // Reflect off boundaries
        while (parts[i].x < 0 || parts[i].x > size) {
            parts[i].x = parts[i].x < 0 ? -parts[i].x : 2 * size - parts[i].x;
            parts[i].vx = -parts[i].vx;
        }
        while (parts[i].y < 0 || parts[i].y > size) {
            parts[i].y = parts[i].y < 0 ? -parts[i].y : 2 * size - parts[i].y;
            parts[i].vy = -parts[i].vy;
        }

        int new_bin = get_bin_index(parts[i].x, parts[i].y);
        new_bins[new_bin].push_back(i);
    }

    bins = new_bins;  // Update bins
}
