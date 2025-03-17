#!/usr/bin/env python3
"""
Simple visualization script for Energy Correlator outputs.
This script reads the histogram data produced by the example_enc_calculator
and creates plots to visualize the results.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def load_histogram(filename):
    """Load histogram data from the given file."""
    # Initialize variables
    data = {}

    # Read the file
    with open(filename, 'r') as f:
        content = f.read()

    # Extract variables by running the content as Python code
    exec(content, globals(), data)

    return data


def plot_penc(hist_data, filename):
    # Plotting
    plt.plot(hist_data['theta1_centers'], hist_data['hist'])

    # TODO: process

    # Saving
    plt.savefig(f"{filename.replace('.py', '')}.png")
    plt.close()

def plot_re3c(hist_data, filename):
    # Plotting
    # TODO: implement
    plt.plot(hist_data['theta1_centers'], hist_data['hist'])

    # TODO: process

    # Saving
    plt.savefig(f"{filename.replace('.py', '')}.png")
    plt.close()



def main():
    """Main function to create visualizations from histogram data."""
    # Find the output files
    output_files = [f for f in os.listdir(os.getcwd()) if
                    (f.startswith("penc_example") or
                     f.startswith('re3c_example'))
                    and f.endswith(".py")]

    if not output_files:
        print("Error: No output files found in the output directory.")
        return 1

    print(f"Found {len(output_files)} output files.")

    # Process each file
    for filename in output_files:
        print(f"Processing file: {filename}")

        # Load the histogram data
        try:
            hist_data = load_histogram(filename)
        except Exception as e:
            print(f"Error loading histogram from {filename}: {e}")
            continue

        if filename.startswith("penc_example"):
            plot_penc(hist_data, filename)
        else:
            plot_re3c(hist_data, filename)
        print(f"Created visualizations for {filename}")

    print(f"All visualizations saved")
    return 0

if __name__ == "__main__":
    main()
