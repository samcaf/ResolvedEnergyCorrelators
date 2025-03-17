# Resolved Energy Correlator Example

This directory contains a minimal working example to calculate Projected Energy Correlators (PENC) and Resolved Energy Correlators (R3NC) as presented in our paper.

## Requirements

- C++ compiler (GCC or Clang)
- FastJet (version 3.0+)
- Pythia 8
- Python 3.6+ with matplotlib and numpy (for visualization)

## Build Instructions

To build the example, run:

```bash
make
```

This will compile the example energy correlator calculator.

## Running the Example

To run a complete example calculation, use:

```bash
./run_example.sh
```

This script will:
1. Run the energy correlator calculation on a sample dataset
2. Generate histograms in the `output` directory
3. Produce visualization plots

## Manual Execution

If you prefer to run the example manually:

1. Calculate correlators:
   ```
   ./example_enc_calculator --weights 1 0 --n_events 10000 --file_prefix example_run
   ```

2. Visualize results:
   ```
   python visualize_example.py
   ```

## Input Parameters

The main executable accepts the following parameters:
- `--weights`: Energy weighting parameters (pairs of values)
- `--n_events`: Number of events to process
- `--file_prefix`: Prefix for output files
- `--use_deltaR`: Use ΔR instead of angle (default: true for pp collisions)
- `--use_pt`: Use pT instead of energy (default: true for pp collisions)
- `--nbins`: Number of bins in the histogram (default: 100)
- `--minbin`: Minimum bin in log scale (default: -8)
- `--maxbin`: Maximum bin in log scale (default: 0.05)

## Output Format

The program produces histogram files in the `output` directory with the following format:
- Each file contains bin edges and centers for θ₁, θ₂/θ₁, and ϕ dimensions
- The main histogram data is stored as a 3D array
- Files are saved in a Python-readable format
