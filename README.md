# DelphiScribe

DelphiScribe a project for obtaining energy-weighted correlations in Pythia and on CMS Open Data.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Introduction

DelphiScribe is a project dedicated to computing energy-weighted correlation in particle physics datasets.
It facilitates analysis on both simulated data from Pythia and real data from the CMS Open Data.


## Features

- **New Angles on Energy Correlators**:

  Introduces a new parameterization for N-point Energy Correlators (ENCs) as detailed in [2410.xxxx].

## Installation

To install DelphiScribe, follow these steps:

1. **Clone the Repository**

   ```
   bash git clone https://github.com/yourusername/DelphiScribe.git
   ```

2. **Navigate to the Directory**

    ```
    cd DelphiScribe
    ```

3. **Configure the Makefile**

    Before compiling the code, you'll need to edit `Makefile.inc` to set up the necessary directories:
    * Open `Makefile.inc` in your preferred text editor.
    * Modify the following variables to match the installation paths and versions on your system:
          + `SOFTWARE_DIR`
          + `PYTHIA_DIR`
          + `FASTJET_DIR`

4. **Configure Open Data File Path**

    If you intend to use information from  the CMS Open Data 2011A Jet Primary Dataset:
    * Open `write/include/opendata_utils.h`.
    * Change the `cms_jets_file` variable to point to the location of the CMS Open Data file on your machine.

5. **Compile the Code**

    Run `make` to compile all relevant C++ code and set up a Python virtual environment
    ```
    make
    ```

6. **Create Output Directories**

    Create the necessary directories for storing output data and figures:
    ```
    mkdir -p output/new_encs
    mkdir -p output/new_enc_figures
    ```


## Usage

After installation, you can start computing Energy Correlators using the following commands.


### New Angles on Energy Correlators

To generate files containing N-Point Energy Correlators (ENCs) with the keyword `opendata_test` in the directory `./output/new_encs/`, try running the following commands.
Additional examples, including examples for computing ENCs in Pythia, can be found in `./bin/`.

* Projected ENCs (PENCs)

Generate PENCs by running:

```
./write/new_enc/2particle --use_opendata true --use_deltaR --use_pt --weights 1.0 [or any list] --n_events 100000 --nbins 200 --file_prefix opendata_test
```

* Resolved 3-Point ENCs (RE3Cs)

Generate RE3Cs with:

```
./write/new_enc/3particle --use_opendata true --use_deltaR --use_pt --weights 1.0 1.0 [or any set of pairs] --n_events 100000 --nbins 150 --file_prefix opendata_test
```

Resolved 4-Point ENCs (RE4Cs)

Generate RE4Cs using:
```
./write/new_enc/4particle --use_opendata true --use_deltaR --use_pt --weights 1.0 1.0 1.0 [or any set of triples] --n_events 100000 --nbins 150 --file_prefix opendata_test
```

Each of these commands will produce associated files in output/new_encs. You can use the plotting tools in plot/encs, which can be modified to produce your own versions of the plots from [2410.xxxx].

## Project Structure

The project is organized into several directories, each serving a specific purpose:

write/: Contains the main C++ source code for computing Energy Correlators.

* new_enc/: Executables for computing Projected ENCs, Resolved 3-Point ENCs, and Resolved 4-Point ENCs.

* src/: Core C++ source files.

* utils/: Utility functions and classes for data processing.

* include/: Header files defining interfaces and data structures.


plot/: Houses Python scripts for data visualization.

* histogram.py: Contains a histogram class which is useful for plotting ENCs, testing their normalization, integrating over variables, finding sub-histograms, etc..

* plotter.py: Plotter class inspired by the MIT Open Data plot format.

* utils/: Python utility modules.

* encs/: Project-specific plotting code for ENCs.


bin/: Compiled binaries and executable files.


output/: Directory for storing output data and generated figures.

* new_encs/: Output data files from ENC computations.

* new_enc_figures/: Figures and plots generated from the data.

Makefile.inc: Configuration file for compilation settings.

README.md: Project documentation.

## Contributing

We welcome contributions from the community! To contribute:

1. Fork the Repository

    Click on the 'Fork' button at the top right corner of the repository page.

2. Create a New Branch

    ```
    git checkout -b feature/YourFeature
    ```

3. Commit Your Changes

    ```
    git commit -m "Add your feature"
    ```

4. Push to Your Branch

    ```
    git push origin feature/YourFeature
    ```

5. Open a Pull Request

    Navigate to the original repository and click on 'New Pull Request'.


## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Contact
For any questions or suggestions:

* **Email**: samuelaf@mit.edu

* **GitHub Issues**: Issue Tracker
