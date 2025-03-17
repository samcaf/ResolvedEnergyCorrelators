// ---------------------------------
// Basic imports
// ---------------------------------
#include <iostream>
#include <cmath>
#include <locale>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include <map>
#include <utility>
#include <stdexcept>

#include <chrono>
using namespace std::chrono;

// for including infinity as an overflow bin
#include <limits>

// ---------------------------------
// HEP imports
// ---------------------------------
#include "fastjet/PseudoJet.hh"

// Local imports:
#include "../../include/general_utils.h"
#include "../../include/jet_utils.h"
#include "../../include/cmdln.h"
#include "../../include/pythia_cmdln.h"

#include "../../include/enc_utils.h"

#include "../../include/opendata_utils.h"


// Type definition for histograms
typedef std::vector<double> Hist;
// (we are only differential in a single angle)


// ####################################
// Main
// ####################################
/**
* @brief: Generates events with Pythia and creates EWOC histograms.
*
* @return: int
*/
int main (int argc, char* argv[]) {
    // Starting timer
    auto start = high_resolution_clock::now();
    const int n_events = cmdln_int("n_events",
                                   argc, argv,
                                   100000);

    // Getting the list of energy weights from command line
    // (We are calculating the projected `PENC`)
    std::vector <double> n_weights;
    for(int iarg=0; iarg<argc; ++iarg) {
        if(str_eq(argv[iarg], "--weights"))
            while (iarg+1 < argc and
                    // next arg doesn't start with '--'
                   std::string(argv[iarg+1]).find("--") == std::string::npos) {
                ++iarg;
                n_weights.emplace_back(atof(argv[iarg]));
            }
    }

    if (n_weights.size() == 0)
        throw std::invalid_argument(
            "Must be given at least 1 weight (via `--weights`).");

    // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
    // Histogram Settings
    // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
    const int    nbins   = cmdln_int("nbins", argc, argv, 100,
                                     false);
    // NOTE: using logarithmically spaced bins,
    // NOTE:   with minbin and maxbin in base 10

    const double minbin  = cmdln_double("minbin", argc, argv,
                                  -8, false);
    const double maxbin  = cmdln_double("maxbin", argc, argv,
                                  1, false);
    const bool   uflow = true, oflow = true;


    // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    // Initializing bin edges and centers
    // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    // - - - - - - - - - - - - - - -
    // For theta1
    // - - - - - - - - - - - - - - -
    // Setting up under/overflow
    int nbins_finite = nbins, bins_finite_start = 0;
    if (uflow) {nbins_finite -= 1; bins_finite_start += 1;}
    if (oflow) nbins_finite -= 1;

    // (logarithmic, from 10^minbin to 10^maxbin)
    const std::vector<double> bin_edges   = get_bin_edges(
                                            minbin, maxbin, nbins,
                                            uflow, oflow);
    const std::vector<double> bin_centers = get_bin_centers(
                                            minbin, maxbin, nbins,
                                            uflow, oflow);


    // =====================================
    // Output Setup
    // =====================================
    // Set up histograms
    std::vector<Hist> enc_hists;
    // Set up histogram output files
    std::vector<std::string> enc_outfiles;

    for (auto n : n_weights){
        // Setting up histograms
        enc_hists.emplace_back(Hist (nbins));

        // Setting up output files,
        std::string filename = "penc_example_n-" +
                str_round(n, 2);

        filename = periods_to_hyphens(filename);
        filename += ".py";
        // writing a header with relevant information
        // for a PENC projected to depend on only a single angle,
        write_enc_header(filename, argc, argv,
                         std::vector<double> {n},
                         true);
        // and adding them to the dict of output files
        enc_outfiles.push_back(filename);
    }

    // ---------------------------------
    // CMS Open Data
    // ---------------------------------
    od::EventReader cms_jet_reader(od::cms_jets_file);


    // ---------------------------------
    // =====================================
    // Analyzing events
    // =====================================
    // ---------------------------------
    // Initializing good_jets, sorted angles and weights
    std::vector<PseudoJet> good_jets;
    std::vector<std::pair<double, double>> sorted_angsweights;

    // Reserving memory
    good_jets.reserve(5);
    sorted_angsweights.reserve(50);

    // =====================================
    // Looping over events
    // =====================================
    for (int iev = 0; iev < n_events; ++iev){
        progressbar(static_cast<double>(iev+1)/
                    double(n_events));
        // -----------------------------------------
        // Reading jet from CMS Opendata
        // -----------------------------------------
        PseudoJet jet;
        cms_jet_reader.read_jet(jet);

        // Storing jet constituents
        const std::vector<PseudoJet>& constituents = jet.constituents();
        double weight_tot = 0;
        for (const auto& particle : constituents)
            weight_tot += particle.pt();

        // ---------------------------------
        // Loop on "special" particle
        for (const auto& part_sp : constituents) {
            // Energy-weighting factor for "special" particle
            double weight_sp = part_sp.pt() / weight_tot;
            // Initializing sum of weights
            // within an angle of 1st particle
            double sum_weight1 = weight_sp;
            // At particle j within the loop below,
            // sum_weight1 = \sum_{thetak < thetaj} weight1_k

            // Preparing contact terms for each PENC weight:
            for (size_t iweight = 0; iweight < n_weights.size(); ++iweight) {
                enc_hists[iweight][0] += std::pow(sum_weight1,
                                          1+n_weights[iweight]);
            }

            // (Sorting:
            //   * [theta1]: angle relative to special particle
            //   * [weight1]: either E2/Ejet or pt2/ptjet
            //  by theta1)
            sorted_angsweights.clear();

            // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
            // Loop on particles
            for (const auto& part1 : constituents) {
                // Energy-weighting factor for particle 1
                double weight1 = part1.pt() / weight_tot;

                // Angle relative to "special" particle
                double theta1 = part_sp.delta_R(part1);

                sorted_angsweights.emplace_back(theta1, weight1);
            } // end second particle loop
            // Sorting angles/weights by angle as promised :)
            std::sort(sorted_angsweights.begin(),
                      sorted_angsweights.end());

            // Loop on second particle
            // (calculating change in cumulative PENC)
            for (size_t jpart=1; jpart<sorted_angsweights.size(); ++jpart){
                double theta1 = sorted_angsweights[jpart].first;
                double weight1 = sorted_angsweights[jpart].second;

                // Calculating the theta1 bin in the histogram
                int bin = bin_position(theta1, minbin, maxbin,
                                       nbins, "log",
                                       uflow, oflow);

                // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
                // Looping on PENC weights [`n's]
                for (size_t iweight = 0; iweight < n_weights.size(); ++iweight) {
                    double hist_weight = weight_sp*(
                                std::pow(sum_weight1+weight1,
                                         n_weights[iweight])
                                -
                                std::pow(sum_weight1,
                                          n_weights[iweight])
                            );

                    enc_hists[iweight][bin] += hist_weight;
                }
                // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

                // Preparing for the next term in the loop!
                sum_weight1 += weight1;
            } // end calculation/particle loop
            // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        } // end "special particle" loop
        // ---------------------------------
    } // end event loop
    // =====================================


    // -----------------------------------
    // =====================================
    // Writing output files
    // =====================================
    // -----------------------------------
    for (size_t iweight = 0; iweight < n_weights.size(); ++iweight) {
        // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        // Output setup
        // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-
        // Opening histogram output file
        // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-
        std::string filename = enc_outfiles[iweight];
        std::fstream outfile;
        outfile.open(filename, std::ios_base::in |
                               std::ios_base::out |
                               std::ios_base::app);

        // Checking for existence
        if (!outfile.is_open()) {
            std::stringstream errMsg;
            errMsg << "File for EnC output was expected "
                   << "to be open, but was not open.\n\n"
                   << "It is possible the file was unable to "
                   << "be created at the desired location:\n\n\t"
                   << "filename = " << filename << "\n\n"
                   << "Is the filename an absolute path? If not, "
                   << "that might be the problem.";
            throw std::runtime_error(errMsg.str().c_str());
        }

        // =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
        // Writing bins to files
        // =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
        // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        // theta1s
        // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        // -:-:-:-:-:-:-:-:-:-:-:-:
        // bin edges
        // -:-:-:-:-:-:-:-:-:-:-:-:
        outfile << "theta1_edges = [\n\t";

        // nbins+1 bin edges:
        //   include -infty and infty for under/overflow
        for (int ibin = 0; ibin < nbins; ++ibin)
            outfile << std::pow(10, bin_edges[ibin]) << ", ";
        if (std::isinf(bin_edges[nbins]))
            outfile << "np.inf\n";
        else
            outfile << std::pow(10, bin_edges[nbins]) << "\n";

        outfile << "]\n\n";

        // -:-:-:-:-:-:-:-:-:-:-:-:
        // bin centers
        // -:-:-:-:-:-:-:-:-:-:-:-:
        outfile << "theta1_centers = [\n\t";

        for (int ibin = 0; ibin < nbins-1; ++ibin)
            outfile << std::pow(10, bin_centers[ibin]) << ", ";
        if (std::isinf(bin_centers[nbins-1]))
            outfile << "np.inf\n";
        else
            outfile << std::pow(10, bin_centers[nbins-1]) << "\n";

        outfile << "]\n\n";


        // =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
        // Processing/writing histogram
        // =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
        // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-
        // Normalizing histogram
        // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-
        // Currently, hist contains
        //   hist[ibin] = N_jets * d Sigma[theta1]
        // Now, changing all finite bins:
        //   hist[ibin] -> (dSigma/dtheta1)
        // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-
        double total_sum = 0.0;

        // Looping over all bins
        for (int bin=0; bin < nbins; ++bin) {
            // Dealing with expectation value over N jets
            enc_hists[iweight][bin] /= n_events;
            total_sum += enc_hists[iweight][bin];

            // Not normalizing outflow bins further
            if (bin < bins_finite_start or bin >= nbins_finite)
                continue;

            // Otherwise, getting differential "volume" element
            double dlogtheta1 = (bin_edges[bin+1] - bin_edges[bin]);
            if (dlogtheta1 == 0) {
                throw std::runtime_error("Found invalid bin "
                                         "width dlogtheta1=0.");
            }
            enc_hists[iweight][bin] /= dlogtheta1;

            // NOTE: This is theta1 times the
            // NOTE:    linearly normed distribution
        }

        // Printing out weight information
        float total_integral = 0;
        // Looping over all bins
        for (int bin=0; bin < nbins; ++bin) {
            // Not normalizing outflow bins further
            if (bin < bins_finite_start or bin >= nbins_finite) {
                total_integral += enc_hists[iweight][bin];
                continue;
            }

            // Otherwise, getting differential "volume" element
            double dlogtheta1 = (bin_edges[bin+1] - bin_edges[bin]);
            if (dlogtheta1 == 0) {
                throw std::runtime_error("Found invalid bin "
                                         "width dlogtheta1=0.");
            }

            total_integral += enc_hists[iweight][bin]*dlogtheta1;
        }

        // Printing normalization
        double n = n_weights[iweight];
        std::cout << "\nTotal weight for n=" << n << ": "
                  << total_sum;
        std::cout << "\nIntegrated weight for n=" << n << ": "
                  << total_integral;

        // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-
        // Writing finalized histogram
        // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-
        outfile << "hist = [\n\t";

        // loop over theta1s
        for (int bin = 0; bin < nbins; ++bin) {
            outfile << std::setprecision(10)
                    << enc_hists[iweight][bin];
            if (bin != nbins-1)
                outfile << ", ";
            else
                outfile << "\n";
        }
        outfile << "]";

        // Closing the file
        outfile.close();
    }


    // ---------------------------------
    // Verifying successful run
    // ---------------------------------
    std::cout << "\nComplete!\n";
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop-start);
    std::cout << "Analyzed and saved data from "
              << std::to_string(n_events)
              << " events in "
              << std::to_string(float(duration.count())/std::pow(10, 6))
              << " seconds.\n";

    return 0;
}
