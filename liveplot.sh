#!/bin/bash
# This script continuously runs the link_distance.py script and updates the plot.

# Run an infinite loop
while true; do
    # Wait for 2 seconds
    sleep 2
    # Clear the terminal
    clear
    # Run the link_distance.py script
    python3 link_dist.py
    # Print the current date and time
    echo "Updating plot at $(date +'%m/%d/%Y %I:%M:%S %p')"
done