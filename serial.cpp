#include "common.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <unordered_set> // Used for efficient cross-boundary tracking

// Global variables for binning
int bin_count;              // Number of bins per row/column
double bin_size;            // Size of each bin (equal to cutoff distance)
std::vector<std::vector<int>> bins;  // Bins storing particle indices
std::vector<std::vector<int>> neighbors;  // Precomputed valid neighbor bins for each bin

/**
 * Converts a particle's (x, y) position to a bin index.
 */
int get_bin_index(double x, double y) {
    int bx = static_cast<int>(x / bin_size);
    int by = static_cast<int>(y / bin_size);
    return bx * bin_count + by;
}

/**
 * Computes the force exerted between two particles and updates both.
 * Uses Newton's Third Law: F(A on B) = -F(B on A).
 */
void apply_force(particle_t& particle, particle_t& neighbor) {
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    if (r2 > cutoff * cutoff || r2 == 0) return; // Ignore if too far or same particle

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Compute repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    
    // Update accelerations for both particles (Newton’s Third Law)
    double ax = coef * dx;
    double ay = coef * dy;

    particle.ax += ax;
    particle.ay += ay;
    neighbor.ax -= ax; // Opposite direction
    neighbor.ay -= ay;
}

/**
 * Updates the position and velocity of a particle using the Velocity Verlet method.
 * Returns `true` if the particle **crossed a bin boundary**.
 */
bool move(particle_t& p, double size, int& old_bin) {
    int new_bin = get_bin_index(p.x, p.y); // Compute bin before move

    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Reflect particles off the walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }
    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }

    old_bin = new_bin; // Store original bin
    return (new_bin != get_bin_index(p.x, p.y)); // Return true if crossed a bin
}

/**
 * Initializes the binning structure and precomputes neighbor bins.
 */
void init_simulation(particle_t* parts, int num_parts, double size) {
    bin_size = cutoff;
    bin_count = (int)ceil(size / bin_size);
    bins.resize(bin_count * bin_count);
    neighbors.resize(bin_count * bin_count);

    // Precompute valid neighbor bins for each bin  // trivial to parallelize
    for (int bx = 0; bx < bin_count; bx++) {
        for (int by = 0; by < bin_count; by++) {
            int bin_index = bx * bin_count + by;

            // Compute the pre-selected 4 neighbors: left, top, top-left, top-right
            if (bx > 0) neighbors[bin_index].push_back((bx - 1) * bin_count + by);     // Left
            if (by < bin_count - 1) neighbors[bin_index].push_back(bx * bin_count + (by + 1)); // Top
            if (bx > 0 && by < bin_count - 1) neighbors[bin_index].push_back((bx - 1) * bin_count + (by + 1)); // Top-left
            if (bx < bin_count - 1 && by < bin_count - 1) neighbors[bin_index].push_back((bx + 1) * bin_count + (by + 1)); // Top-right
        }
    }

    // Assign each particle to a bin  // trivial to parallelize
    for (int i = 0; i < num_parts; i++) {
        int bin_index = get_bin_index(parts[i].x, parts[i].y);
        bins[bin_index].push_back(i);
    }
}

/**
 * Simulates a single time step of the particle system.
 */
void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // Step 1: Compute forces using precomputed neighbors & avoiding redundancy  // NOT SURE HOW WE CAN AVOID RACE CONDITIONS HERE
    for (int bx = 0; bx < bin_count; bx++) {
        for (int by = 0; by < bin_count; by++) {
            int bin_index = bx * bin_count + by;

            for (int i : bins[bin_index]) {
                parts[i].ax = parts[i].ay = 0; // Reset acceleration
            }

            // Apply forces within current bin
            for (int i : bins[bin_index]) {   // Use critical section here
                for (int j : bins[bin_index]) {
                    if (i < j) apply_force(parts[i], parts[j]); // Avoid duplicate calculations
                }

                // Apply forces with precomputed neighbor bins
                for (int neighbor_bin : neighbors[bin_index]) { // use critical section here
                    for (int j : bins[neighbor_bin]) {
                        apply_force(parts[i], parts[j]);
                    }
                }
            }
        }
    }

    // Use barrier here

    // Step 2: Move particles and update bin assignments **only for crossing particles**  // Particle movement is trivial to parallelize. WILL THE BOUNDARY CROSSING PART CAUSE RACE CONDITION?
    for (int i = 0; i < num_parts; i++) {
        int old_bin;
        if (move(parts[i], size, old_bin)) { // If the particle moved to a different bin
            int new_bin = get_bin_index(parts[i].x, parts[i].y);

            bins[old_bin].erase(remove(bins[old_bin].begin(), bins[old_bin].end(), i), bins[old_bin].end());

            
            bins[new_bin].push_back(i);

        }
    }
}
