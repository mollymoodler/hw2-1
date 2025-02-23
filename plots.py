import numpy as np
import matplotlib.pyplot as plt

# Load data from file
data = np.loadtxt("serial_results.txt")
num_particles = data[:, 0]
times = data[:, 1]

# Log-log plot
plt.figure(figsize=(8,6))
plt.loglog(num_particles, times, marker='o', linestyle='-', label="Serial Runtime")

# Fit a line to check O(n) complexity
fit = np.polyfit(np.log(num_particles), np.log(times), 1)
plt.loglog(num_particles, np.exp(fit[1]) * num_particles**fit[0], '--', label=f'O(n) fit (slope={fit[0]:.2f})')

plt.xlabel("Number of Particles (log scale)")
plt.ylabel("Execution Time (log scale)")
plt.title("Serial Simulation Performance")
plt.legend()
plt.grid(True)
plt.show()
