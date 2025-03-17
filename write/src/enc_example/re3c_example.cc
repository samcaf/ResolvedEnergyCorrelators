/**
 * @file    new_enc_3particle.cc
 *
 * @brief   Code for generating histograms for "new angles on"
 *          n-point projected energy correlators (ENCs);
 *          in this file, we consider three particles, and are
 *          differential in three angles:
 *          theta1, theta2, and the azimuthal angle phi_{2*1}.
 */


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


// =====================================
// Type definitions for histograms
// =====================================
// Multi-dimensional Histograms
typedef std::vector<double> Hist1d;
typedef std::vector< std::vector<double>> Hist2d;
typedef std::vector< std::vector< std::vector<double>>> Hist3d;

// Using pairs of weights to specify the doubly-projected correlator
typedef std::pair<double, double> weight_t;


// =====================================
// Additional Utilities
// =====================================
inline double mod2pi(double phi) {
    while (phi > PI)
        phi -= TWOPI;
    while (phi <= -PI)
        phi += TWOPI;

    return phi;
}

double enc_azimuth(const PseudoJet part1,
                   const PseudoJet part_sp,
                   const PseudoJet part2) {
    // NOTE: I think this doesn't work for very fat jets --
    // NOTE:   roughly because fat jets require a bit more care
    // NOTE:   with mod(2pi) arithmetic
    // - - - - - - - - - - - - - - - - - - -
    // Normalized vectors in eta-phi plane
    // - - - - - - - - - - - - - - - - - - -
    // From part_sp to part1
    double x1 = part1.rap() - part_sp.rap();
    double y1 = mod2pi(part1.phi() - part_sp.phi());
    const double n1  = sqrt(std::pow(x1, 2.) + std::pow(y1, 2.));
    if (n1 == 0)
        return 0;
    x1 /= n1; y1 /= n1;

    // From part_sp to part2
    double x2;
    x2 = part2.rap() - part_sp.rap();
    double y2 = mod2pi(part2.phi() - part_sp.phi());
    const double n2  = sqrt(std::pow(x2, 2.) + std::pow(y2, 2.));
    if (n2 == 0)
        return 0;
    x2 /= n2; y2 /= n2;

    // - - - - - - - - - - - - - - - - - - -
    // Getting the "angle" between these vectors
    // - - - - - - - - - - - - - - - - - - -
    const double dot = x1*x2 + y1*y2;
    const double det = x1*y2 - y1*x2;
    double phi = mod2pi(atan2(det, dot));

    // Setting it to be between -pi and pi
    phi = phi > PI ? phi - TWOPI : phi;
    if (phi < -PI or PI < phi)
        throw std::range_error(
                "Found azimuthal angle not between -pi and pi.");

    return phi;
}


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

    // =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
    // Basic Pythia Settings
    // =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
    // 50k e+ e=:=> hadrons events, by default
    const int         n_events      = cmdln_int("n_events",
                                          argc, argv,
                                          100000);

    // =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
    // Histogram Settings
    // =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
    const int   nbins     = cmdln_int("nbins", argc, argv,
                                      100, false);
    // NOTE: using logarithmically spaced bins,
    // NOTE:   with minbin and maxbin in base 10

    const double minbin   = cmdln_double("minbin", argc, argv,
                                   -8, false);
    const double maxbin   = cmdln_double("maxbin", argc, argv,
                                   0.05, false);
    const bool bin1_uflow = true, bin1_oflow = true;

    // Phi is binned linearly, with same nbins by default
    const int   nphibins  = cmdln_int("nphibins", argc, argv,
                                nbins, false);

    // theta2/theta1 is binned linearly by default, but can be log
    const bool lin_bin2   = cmdln_bool("lin_bin2", argc, argv,
                                 true, false);

    // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    // Initializing bin edges and centers
    // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    // - - - - - - - - - - - - - - -
    // For theta1
    // - - - - - - - - - - - - - - -
    // Setting up under/overflow
    int nbins1_finite = nbins, bin1_finite_start = 0;
    if (bin1_uflow) {nbins1_finite -= 1; bin1_finite_start += 1;}
    if (bin1_oflow) nbins1_finite -= 1;

    // (logarithmic, from 10^minbin to 10^maxbin)
    const std::vector<double> bin1_edges   = get_bin_edges(
                                            minbin, maxbin, nbins,
                                            bin1_uflow, bin1_oflow);
    const std::vector<double> bin1_centers = get_bin_centers(
                                            minbin, maxbin, nbins,
                                            bin1_uflow, bin1_oflow);

    // - - - - - - - - - - - - - - -
    // Bins for theta2/theta
    // - - - - - - - - - - - - - - -
    // (from 0 to 1, with a variable bin-spacing scheme)
    // Options dependent on the binning of theta2/theta1:
    const std::string bin2_scheme = lin_bin2 ? "lin" : "log";
    // t2/t1 in (0, 1) in linear bins, or
    // t2/t1 in (1e-minbin, 1) if logarithmic bins
    const double bin2_min = lin_bin2 ?   0   : minbin;
    const double bin2_max = lin_bin2 ?   1   : 0;
    // Use underflow only if binning logarithmically
    const bool bin2_uflow        = lin_bin2 ? false : true;
    const int  bin2_finite_start = lin_bin2 ?   0   : 1;
    const int  nbins2_finite     = lin_bin2 ? nbins : nbins-1;

    const std::vector<double> bin2_edges = get_bin_edges(
                                        bin2_min, bin2_max,
                                        nbins, bin2_uflow, false);
                            /* underflow = (!bin_lin2) */
                            /* but overflow = false */
    const std::vector<double> bin2_centers = get_bin_centers(
                                        bin2_min, bin2_max,
                                        nbins, bin2_uflow, false);

    // - - - - - - - - - - - - - - -
    // For "azimuthal" angle phi
    // - - - - - - - - - - - - - - -
    // (linear, from -pi to pi)
    const std::vector<double> phi_edges   = get_bin_edges(
                                                -PI, PI, nphibins,
                                                false, false);
    const std::vector<double> phi_centers = get_bin_centers(
                                                -PI, PI, nphibins,
                                                false, false);
    // Bin for phi = 0
    const int phizerobin = bin_position(0, -PI, PI, nphibins,
                                        "linear", false, false);

    // =====================================
    // Output Setup
    // =====================================
    // Set up histograms
    Hist3d enc_hist (nbins, Hist2d(nbins, Hist1d(nphibins)));
    // Set up histogram output files
    std::string enc_outfile = "re3c_example.py";

    // writing a header with relevant information
    write_enc_header(enc_outfile, argc, argv,
                 std::vector<double> {1, 1},
                 true);

    // ---------------------------------
    // CMS Open Data
    // ---------------------------------
    od::EventReader cms_jet_reader(od::cms_jets_file);

    // ---------------------------------
    // =====================================
    // Analyzing events
    // =====================================
    // ---------------------------------
    // Initializing particles, good_jets, sorted angles and weights
    std::vector<PseudoJet> good_jets;
    std::vector<std::pair<double, PseudoJet>> sorted_angs_parts;

    // Reserving memory
    good_jets.reserve(5);
    sorted_angs_parts.reserve(50);

    // =====================================
    // Looping over events
    // =====================================
    for (int iev = 0; iev < n_events; ++iev) {
        progressbar(static_cast<double>(iev+1)/
                    double(n_events));
        // -----------------------------------------
        // CMS Open Data (gives jets from the start)
        // -----------------------------------------
        PseudoJet jet;
        cms_jet_reader.read_jet(jet);
        good_jets.emplace_back(std::move(jet));

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

            enc_hist[0][0][phizerobin] +=
                    std::pow(weight_sp, 3);

            // (Sorting:
            //   * [theta1]: angle relative to special particle
            //   * [weight1]: either E2/Ejet or pt2/ptjet
            //  by theta1)
            sorted_angs_parts.clear();

            // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
            // Loop on particles
            for (const auto& part1 : constituents) {
                // Angle relative to "special" particle
                double theta1 = part_sp.delta_R(part1);

                sorted_angs_parts.emplace_back(theta1, part1);
            } // end second particle loop
            // Sorting angles/weights by angle as promised :)
            std::sort(sorted_angs_parts.begin(),
                      sorted_angs_parts.end(),
                      [](auto& left, auto& right) {
                          return left.first < right.first;
                     });
            // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

            // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
            // Loop on first non-special particle
            // (calculating change in cumulative E^nu C)
            for (size_t jpart=1; jpart<sorted_angs_parts.size(); ++jpart) {
                // Properties of 1st particle
                double theta1    = sorted_angs_parts[jpart].first;
                PseudoJet& part1 = sorted_angs_parts[jpart].second;
                double weight1 = part1.pt() / weight_tot;

                // Calculating the theta1 bin in the histogram
                int bin1 = bin_position(theta1, minbin, maxbin,
                                        nbins, "log",
                                        bin1_uflow, bin1_oflow);

                // Initializing the sum of weights
                // within an angle of the 2nd non-special particle
                std::vector<double> sum_weight2(nphibins);

                sum_weight2[phizerobin] += weight_sp;

                // Preparing contact terms:
                // part2 = part_sp != part_1
                enc_hist[bin1][0][phizerobin] +=
                    2*std::pow(weight_sp, 2)*
                      std::pow(weight1, 1);

                // part2 = part1 != part_sp
                enc_hist[bin1][nbins-1][phizerobin] +=
                            std::pow(weight_sp, 1)*
                            std::pow(weight1, 2);

                // -----------------------------------
                // Loop on second non-special particle
                for (size_t kpart=1; kpart<jpart; ++kpart) {
                    // Getting 2nd particle
                    double theta2    = sorted_angs_parts[kpart].first;
                    PseudoJet& part2 = sorted_angs_parts[kpart].second;
                    double weight2 = part2.pt() / weight_tot;
                    double theta2_over_theta1 =
                        theta1 == 0 ? 0 : theta2/theta1;

                    // Calculating the theta2/theta1 bin position
                    int bin2 = bin_position(theta2_over_theta1,
                                        bin2_min, bin2_max,
                                        nbins, bin2_scheme,
                                        bin2_uflow, false);
                                    /* Variable spacing scheme,
                                     * but with no overflow. */

                    // Getting azimuthal angle
                    // (angle from part1 to part_sp to part2
                    //  in rapidity-azimuth plane)
                    double phi = enc_azimuth(
                            part1, part_sp, part2);

                    // Calculating the phi bin
                    int binphi = bin_position(phi, -PI, PI,
                                          nphibins, "linear",
                                          false, false);

                    // *:*:*:*:*:*:*:*:*:*:*:*:*:*:*:*:*:*:*
                    // Adding to the histogram
                    // *:*:*:*:*:*:*:*:*:*:*:*:*:*:*:*:*:*:*
                    double delta_weight1 = (
                           std::pow(sum_weight1+weight1, 1)
                           -
                           std::pow(sum_weight1, 1)
                         );
                    double delta_weight2 = (
                           std::pow(sum_weight2[binphi]
                                     + weight2, 1)
                           -
                           std::pow(sum_weight2[binphi], 1)
                         );
                    double perm = 2;
                    // imagine a triangle with theta_j < theta_i;
                    // need to count twice to get the full
                    // sum on all pairs (see also contact term)

                    double hist_weight = weight_sp *
                                delta_weight1 *
                                delta_weight2;

                    enc_hist[bin1][bin2][binphi] +=
                            perm*hist_weight;

                    // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
                    // Preparing for the next particle in the loop!
                    sum_weight2[binphi] += weight2;
                } // end calculation/2nd particle loop
                // -----------------------------------

                // Preparing for the particle in the loop!
                sum_weight1 += weight1;
            } // end 1st particle loop
            // -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

        // ---------------------------------
        } // end "special particle" loop
        // ---------------------------------
    } // end event loop
    // =====================================


    // ===================================
    // Writing histograms to output file
    // ===================================
    // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    // Output setup
    // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-
    // Opening histogram output file
    // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-
    std::fstream outfile;
    outfile.open(enc_outfile, std::ios_base::in |
                           std::ios_base::out |
                           std::ios_base::app);

    // Checking for existence
    if (!outfile.is_open()) {
        std::stringstream errMsg;
        errMsg << "File for EnC output was expected "
               << "to be open, but was not open.\n\n"
               << "It is possible the file was unable to "
               << "be created at the desired location:\n\n\t"
               << "filename = " << enc_outfile << "\n\n"
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
        outfile << std::pow(10, bin1_edges[ibin]) << ", ";
    if (std::isinf(bin1_edges[nbins]))
        outfile << "np.inf\n";
    else
        outfile << std::pow(10, bin1_edges[nbins]) << "\n";

    outfile << "]\n\n";

    // -:-:-:-:-:-:-:-:-:-:-:-:
    // bin centers
    // -:-:-:-:-:-:-:-:-:-:-:-:
    outfile << "theta1_centers = [\n\t";

    for (int ibin = 0; ibin < nbins-1; ++ibin)
        outfile << std::pow(10, bin1_centers[ibin]) << ", ";
    if (std::isinf(bin1_centers[nbins-1]))
        outfile << "np.inf\n";
    else
        outfile << std::pow(10, bin1_centers[nbins-1]) << "\n";

    outfile << "]\n\n";


    // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    // theta2_over_theta1s
    // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    // -:-:-:-:-:-:-:-:-:-:-:-:
    // bin edges
    // -:-:-:-:-:-:-:-:-:-:-:-:
    outfile << "theta2_over_theta1_edges = [\n\t";

    // nbins+1 bin edges:
    for (int ibin = 0; ibin < nbins+1; ++ibin) {
        double bin2_edge = lin_bin2 ? bin2_edges[ibin]
                                    : std::pow(10, bin2_edges[ibin]);
        outfile << bin2_edge;
        if (ibin < nbins)
            outfile << ", ";
        else
            outfile << std::endl;
    }

    outfile << "]\n\n";

    // -:-:-:-:-:-:-:-:-:-:-:-:
    // bin centers
    // -:-:-:-:-:-:-:-:-:-:-:-:
    outfile << "theta2_over_theta1_centers = [\n\t";

    for (int ibin = 0; ibin < nbins; ++ibin) {
        double bin2_val = lin_bin2 ? bin2_centers[ibin]
                                   : std::pow(10, bin2_centers[ibin]);
        outfile << bin2_val;
        if (ibin < nbins-1)
            outfile << ", ";
        else
            outfile << std::endl;
    }
    outfile << "]\n\n";


    // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    // phis
    // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    // -:-:-:-:-:-:-:-:-:-:-:-:
    // bin edges
    // -:-:-:-:-:-:-:-:-:-:-:-:
    outfile << "phi_edges = [\n\t";

    // nphibins+1 bin edges:
    for (int ibin = 0; ibin < nphibins; ++ibin)
        outfile << phi_edges[ibin] << ", ";
    outfile << phi_edges[nphibins] << "\n";

    outfile << "]\n\n";

    // -:-:-:-:-:-:-:-:-:-:-:-:
    // bin centers
    // -:-:-:-:-:-:-:-:-:-:-:-:
    outfile << "phi_centers = [\n\t";

    for (int ibin = 0; ibin < nphibins-1; ++ibin)
        outfile << phi_centers[ibin] << ", ";
    outfile << phi_centers[nphibins-1] << "\n";

    outfile << "]\n\n";

    // =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
    // Processing/writing histogram
    // =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
    // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-
    // Normalizing histogram
    // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-
    // Currently, hist contains
    //   hist[ibin] = N_jets * d^3 Sigma[theta1][theta2/theta1][phi]
    // Now, changing all finite bins:
    //   hist[ibin] -> (theta1^2 * d^3Sigma/dtheta1 dtheta2 dphi)
    // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-
    double total_sum = 0.0;

    // Looping over all bins
    for (int bin1=0; bin1<nbins; ++bin1) {
        for (int bin2=0; bin2<nbins; ++bin2) {
            for (int binphi=0; binphi<nphibins; ++binphi) {
                // Dealing with expectation value over N jets
                enc_hist[bin1][bin2][binphi] /= n_events;

                total_sum += enc_hist[bin1][bin2][binphi];

                // Not normalizing outflow bins further
                if (bin1 < bin1_finite_start
                        or bin1 >= nbins1_finite
                        or bin2 < bin2_finite_start
                        or bin2 >= nbins2_finite)
                    continue;

                // Getting differential "volume" element
                double dlogtheta1 = (bin1_edges[bin1+1] - bin1_edges[bin1]);
                double dtheta2_over_theta1 = (bin2_edges[bin2+1] - bin2_edges[bin2]);
                double dphi = (phi_edges[binphi+1] - phi_edges[binphi]);

                double dvol = dlogtheta1 * dtheta2_over_theta1 * dphi;
                enc_hist[bin1][bin2][binphi] /= dvol;

                // NOTE: This is theta1^2 times the
                // NOTE:    linearly normed distribution
            }
        }
    }

    double total_integral = 0.0;
    for (int bin1=0; bin1<nbins; ++bin1) {
        for (int bin2=0; bin2<nbins; ++bin2) {
            for (int binphi=0; binphi<nphibins; ++binphi) {
                if (bin1 < bin1_finite_start
                        or bin1 >= nbins1_finite
                        or bin2 < bin2_finite_start
                        or bin2 >= nbins2_finite) {
                    total_integral += enc_hist[bin1][bin2][binphi];
                    continue;
                }

                // Getting differential "volume" element
                double dlogtheta1 = (bin1_edges[bin1+1] - bin1_edges[bin1]);
                double dtheta2_over_theta1 = (bin2_edges[bin2+1] - bin2_edges[bin2]);
                double dphi = (phi_edges[binphi+1] - phi_edges[binphi]);

                double dvol = dlogtheta1 * dtheta2_over_theta1 * dphi;

                total_integral += enc_hist[bin1][bin2][binphi] * dvol;
            }
        }
    }


    // Printing normalization
    std::cout << "\nTotal weight = " << total_sum;
    std::cout << "\nIntegrated weight = " << total_integral;


    // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-
    // Writing histogram
    // -:-:-:-:-:-:-:-:-:-:-:-:-:-:-
    outfile << "hist = [\n\t";

    // theta1s
    for (int bin1 = 0; bin1 < nbins; ++bin1) {
        outfile << "[\n\t";

        // theta2s
        for (int bin2 = 0; bin2 < nbins; ++bin2) {
            // Phis
            if (nphibins == 1){
                outfile << enc_hist[bin1][bin2][0];
                outfile << (bin2 != nbins-1 ? ", "
                                            : "\n");
            } else {
                outfile << "\t[\n\t\t\t";

                // Loop over phis
                for (int binphi = 0; binphi < nphibins-1; ++binphi) {
                    outfile << std::setprecision(10)
                            << enc_hist[bin1][bin2][binphi] << ", ";
                }
                outfile << enc_hist[bin1][bin2][nphibins-1] << "\n";

                outfile << (bin2 != nbins-1 ? "\t\t],\n\t"
                                            : "\t\t]\n");
            }
        }

        outfile << (bin1 != nbins-1 ? "\t],\n\t"
                                    : "\t]\n");
    }
    outfile << "]";

    outfile.close();

    // ---------------------------------
    // =====================================
    // Verifying successful run
    // =====================================
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
