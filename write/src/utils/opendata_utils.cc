/**
 * @file    opendata_utils.cc
 *
 * @brief   A utility file for reading/processing open data.
 *
 * @author: Samuel Alipour-fard
 */
#include <sstream>
#include "fastjet/PseudoJet.hh"

// Local imports
#include "../../include/jet_utils.h"
#include "../../include/opendata_utils.h"


// =====================================
// Switches, flags, and options
// =====================================
// Cut on jets from the CMS Jet 2011A Dataset
const float CMS_ETA_CUT       = 1.9;
const float CMS_R_JET         = 0.5;
const std::string CMS_JET_ALG = "akt";
const float CMS_PT_MIN        = 500;
const float CMS_PT_MAX        = 550;

namespace od
{

    // =====================================
    // File Reading Utilities
    // =====================================
    /**
    * @brief: Takes in a text file containing processed CMS Open Data.
    *         Reads events into a given vector of PseudoJets when
    *         read_event is called.
    *
    * @param: inputfile   The text file containing CMS Open Data
    */
    EventReader::EventReader(const std::string& inputfile)
                : current_event(0) {
        source.open(inputfile);
        if (!source.good()) {
            throw std::runtime_error("od::EventReader: Error: "
                                     "Input file not found.");
        }
    }

    EventReader::~EventReader() {
        if (source.is_open()) {
            source.close();
        }
    }

    bool EventReader::read_jet(
            fastjet::PseudoJet& jet) {
        std::string line;
        while (std::getline(source, line)) {
            if (line.substr(0, 4) == "#END") {
                return false;
            }
            if (line.substr(0, 1) == "#") {
                continue;
            }

            std::istringstream linestream(line);
            int event_number;
            double pt, eta, phi;
            linestream >> event_number >> pt >> eta >> phi;

            if (linestream.fail()) {
                throw std::runtime_error("od::EventReader: Error:"
                   " Incorrect input format in line: " + line);
            }

            if (event_number > current_event) {
                current_event = event_number;
                return true;
            }

            double px = pt * cos(phi);
            double py = pt * sin(phi);
            double pz = pt * sinh(eta);
            double E = sqrt(px * px + py * py + pz * pz);
            jet = join(jet, fastjet::PseudoJet(px, py, pz, E));
        }
        return false;
    }


    // DEBUG: Old OD method
    void read_events(std::vector< std::vector<PseudoJet> >& events,
                     int nevents, std::string inputfile) {
        std::ifstream source(inputfile);

        if (not source.good())
            throw std::runtime_error("od::read_events:\n\tError: Input file not found.");

        std::string line;
        int event_count = 0;

        while (std::getline(source, line) && event_count < nevents) {
            // Reading the file line-by-line
            std::istringstream linestream(line);

            if (line.substr(0,4) == "#END") {return;}
            if (line.substr(0,1) == "#") {continue;}

            int event_number;
            double pt, eta, phi;
            linestream >> event_number >> pt >> eta >> phi;

            // Check if all variables were successfully read
            if (linestream.fail())
                throw std::runtime_error("od::read_events:\n\tError: "
                                         "Incorrect input format in line: "
                                         + line);

            double px = pt * cos(phi);
            double py = pt * sin(phi);
            double pz = pt * sinh(eta);
            double E = sqrt(px * px + py * py + pz * pz);
            fastjet::PseudoJet particle(px, py, pz, E);
            if (event_number >= nevents) {
                break;
            } else {
                events[event_number].push_back(particle);
            }
        }

        source.close();

        //Check if enough events were provided
        for (int i=0; i<nevents; ++i)
            if (events[i].size()==0)
                throw std::runtime_error("od::read_events:\n\tError: "
                                         "Input file does not have enough events.");
    }
    // END DEBUG

}
