#!/bin/bash

# Define particle sizes to test
particle_counts=(100 500 1000 5000 10000 50000 100000)

# Output file for results
output_file="serial_results.txt"
rm -f $output_file  # Clear previous results

# Run serial simulation for each particle count
for n in "${particle_counts[@]}"; do
    result=$(./serial -n $n | grep "Simulation Time" | awk '{print $5, $8}')
    echo "$result" >> $output_file
done

echo "Results saved to $output_file"
