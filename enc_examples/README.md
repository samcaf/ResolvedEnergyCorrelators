# Resolved Energy Correlator Example

This directory contains a minimal working example to calculate Projected Energy Correlators (PENC) and Resolved Energy Correlators (R3NC) as presented in [PAPER].

## Requirements

- C++ compiler (GCC or Clang)
- FastJet (version 3.0+)
- Python 3.6+ with matplotlib and numpy (for visualization)

## Build Instructions

To download the CMS Open Data dataset used for the examples in this folder, return to the main directory and run

```bash
make get_cms_od
```

Then, return to this directory and run
```bash
chmod +x visualize.py
```


## Running the Example (TODO: IMPLEMENT)

To run a complete example calculation, use:

```bash
./run_examples.sh
```

This script will:
1. Run the energy correlator calculation on a sample dataset
2. Generate histograms in the `output` directory
3. Produce visualization plots


## Manual Execution

If you prefer to run the example manually:

1. Calculate correlators:

    * PENCs:
    ```
    ./write_penc DESCRIPTION
    ```

    - `--weights`: Energy weighting parameters (pairs of values)
    - `--n_events`: Number of events to process
    - `--nbins`: Number of bins in the histogram (default: 100)
    - `--minbin`: Minimum bin in log scale (default: -8)
    - `--maxbin`: Maximum bin in log scale (default: 0.05)

    - Each file contains bin edges and centers for θ₁
    - Files are saved in a Python-readable format
    - The main histogram data is stored as a 1D array


    * RE3C:
    ```
    ./write_re3c
    ```

    Input parameters are the same except weights are not an option anymore (wording)


    - Each file contains bin edges and centers for θ₁, θ₂/θ₁, and ϕ
    - Files are saved in a Python-readable format
    - The main histogram data is stored as a 3D array


2. Visualize results:
   ```
   ./visualize.py
   ```

   - DISCUSSION of output


## Output Format

The program produces histogram files of the form ``penc_example_n-[VALUE].py`` and ``re3c_example.py``, as well as plots of the form ``penc_example_n-[VALUE].png`` and ``re3c_example.png``.
