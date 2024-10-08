# ECScribe

A project for obtaining energy-weighted correlations in Pythia and on CMS Open Data.

<img src="output/display/w_combined_1d.png" width="200"> <img src="output/display/qcd_3particle_bullseye.png" width="200"> <img src="output/display/od_4particle_bullseye.png" width="200"> <img src="output/display/od_nonpert_density.png" width="200">

## Table of Contents


- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)


# Introduction

ECScribe is a project dedicated to computing energy-weighted correlations in particle physics datasets.
It facilitates analysis on both simulated data from Pythia and real data from the CMS Open Data.


# Features

- **Jet Properties**:

  Basic properties of jets -- useful for testing the structure of the C++ code and for seeing how python histograms can be plotted.

- **New Angles on Energy Correlators**:

  Introduces a new parameterization for N-point Energy Correlators (ENCs) as detailed in [2410.xxxx].
            <details>
                <summary>*Pseudocode*</summary>
                <img src="output/display/enc_pseudocode.png" width="600">
            </details>



- **Energy Weighted Observable Correlations (EWOCs)**:

  Not yet included: will introduce a new type of energy correlator on non-angular observables as detailed in [24yy.xxxx].

- **Python Classes**:

  * `./plot/plotter.py` contains a python plotting class inspired by the MIT Open Data plotting format;
  * `./plot/histogram.py` contains a python histogram class with several functionalities, including:
      - integrating out variables;
      - finding sub-histograms;
      - plotting 1-, 2-, and 3-dimensional data;
  * Examples of the usage for these classes can be found, for example, in `./plot/jet_properties`.


## Project Structure
The project is organized into several directories, each serving a specific purpose:

`Makefile.inc`: Configuration file for compilation settings.

`README.md`: Project documentation.

<details>
<summary>
`write/`: Contains the main C++ source code for computing Energy Correlators:
</summary>
* `new_enc/`: Executables for computing Projected ENCs, Resolved 3-Point ENCs, and Resolved 4-Point ENCs;
* `src/`: Core C++ source files;
* `utils/`: Utility functions and classes for data processing;
* `include/`: Header files defining interfaces and data structures;
* `data/`: Houses datasets, including the CMS 2011A Jet Primary Dataset.
</details>

<details>
<summary>
`plot/`: Houses python tools for data visualization:
</summary>
* `plotter.py`: Plotter class inspired by the MIT Open Data plot format.
* `histogram.py`: Contains a histogram class which is useful for plotting ENCs, testing their normalization, integrating over variables, finding sub-histograms, etc..
* `utils/`: Python utility modules.
* `jet_properties/`: Example plotting code for plotting properties of jets from different samples.
* `encs/`: Project-specific plotting code for ENCs.
</details>

<details>
<summary>
`output/`: Directory for storing output data and generated figures.
</summary>
* `new_encs/`: Output data files from ENC computations.
* `new_enc_figures/`: Figures and plots generated from the data.
</details>

`bin/`: Example command-line code for generating events and output.


# Installation

To install ECScribe, follow these steps:

1. **Clone the Repository**
   ```
   git clone https://github.com/samcaf/ECScribe.git
   ```

2. **Navigate to the Directory**

    ```
    cd ECScribe
    ```

3. **Configure the Makefile**

    Before compiling the code, you'll need to edit `Makefile.inc` to set up the necessary directories:
    * Open `Makefile.inc` in your preferred text editor.
    * Modify the following variables to match the installation paths and versions on your system:
      - `SOFTWARE_DIR`
      - `PYTHIA_VERSION`
      - `FASTJET_VERSION`
      - **Note:** If you don't have Pythia or FastJet installed, you can try running the following after setting up the directories/version numbers above:
         ```
         make download pythia
         make install pythia
         ```
         and
         ```
         make download fastjet
         make install fastjet
         ```
         respectively.

4. **Compile the Code**

    Run `make` to compile all relevant C++ code, prepare output directories, download a text file containing the CMS 2011A Jet Primary Dataset (from [FastEEC](https://github.com/abudhraj/FastEEC/releases/tag/0.1)) to `write/data`, and set up a Python virtual environment
    ```
    make
    ```

5. **Configure Open Data File Path**

    If you intend to use information from  the CMS Open Data 2011A Jet Primary Dataset:
    * Open `write/include/opendata_utils.h`.
    * Change the `cms_jets_file` variable to point to the location of the CMS Open Data file on your machine.



# Usage

After installation, you are ready to start computing Energy-Weighted Correlations!


## Jet Properties

To generate histogram files containing jet properties (mass, transverse momentum, pseudorapidity, etc.) with the keyword `opendata_test` in the directory `./output/jet_properties/`, try running the following command:
```
./write/jet_properties --use_opendata true --n_events 100000 --nbins 100 --file_prefix opendata_test
```
You can look at/use/modify the plotting tools in `./plot/jet_properties` for some example plots of jet masses.
Additional examples, including examples for computing ENCs in Pythia, can be found in `./bin/`.


## New Angles on Energy Correlators

To generate files containing N-Point Energy Correlators (ENCs) with the keyword `opendata_test` in the directory `./output/new_encs/`, try running one of the commands below.

You can use the plotting tools in `./plot/encs`, which can be modified to produce your own versions of the plots from [2410.xxxx].
Additional examples for computing ENCs, including examples for computing ENCs in Pythia, can be found in `./bin/`.

### Projected ENCs (PENCs)

<img src="output/display/od_highN_1d.png" width="200"> <img src="output/display/w_combined_1d.png" width="200">

Generate PENCs by running:

```
./write/new_enc/2particle --use_opendata true --use_deltaR --use_pt --weights 1.0 --n_events 100000 --nbins 200 --file_prefix opendata_test
```
The weight 1.0 indicates the energy weight associated with a particle in the jet -- or the value of N-1 for the ENC. It can be replaced by any list of weights (any list of the desired values for N-1);

### Resolved 3-Point ENCs (RE3Cs)

<img src="output/display/qcd_3particle_bullseye.png" width="200"> <img src="output/display/od_newdef_density.png" width="200">

Generate RE3Cs with:

```
./write/new_enc/3particle --use_opendata true --use_deltaR --use_pt --weights 1.0 1.0 --n_events 100000 --nbins 150 --file_prefix opendata_test
```
The weights (1.0, 1.0) indicate the energy weights associated with a pair of resolved particles, and can be changed to any pair or list of pairs;

### Resolved 4-Point ENCs (RE4Cs)

<img src="output/display/od_4particle_bullseye.png" width="200"> <img src="output/display/od_4particle_bullseye-1.png" width="200">

Generate RE4Cs using:
```
./write/new_enc/4particle --use_opendata true --use_deltaR --use_pt --weights 1.0 1.0 1.0 --n_events 100000 --nbins 150 --file_prefix opendata_test
```
The weights (1.0, 1.0, 1.0) can be changed to any list of triples.



## Contributing

We welcome contributions from the community! To [contribute](https://github.com/actions/checkout/blob/main/CONTRIBUTING.md):

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

* **GitHub Issues**: [Issue Tracker](https://github.com/samcaf/ECScribe/issues)
