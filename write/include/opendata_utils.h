#ifndef READ_EVENT_H
#define READ_EVENT_H

#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

#include "fastjet/PseudoJet.hh"

// =====================================
// Switches, flags, and options
// =====================================
// Cut on jets from the CMS Jet 2011A Dataset
extern const float CMS_ETA_CUT;
extern const float CMS_R_JET;
extern const std::string CMS_JET_ALG;
extern const float CMS_PT_MIN;
extern const float CMS_PT_MAX;


namespace od
{
    // Opendata Files
    const std::string cms_jets_file = "/home/samaf/Documents/Research/ResolvedEnergyCorrelators/write/data/cms_jet_run2011A.opendata.txt";

    // Class for reading files
    class EventReader {
    private:
        std::ifstream source;
        int current_event;

    public:
        EventReader(const std::string& inputfile);
        ~EventReader();

        bool read_jet(fastjet::PseudoJet& jet);
    };

    // DEBUG: Old OD method
    // Utilities for reading files
    void read_events(
            std::vector<std::vector<fastjet::PseudoJet>>& events,
            int nevents, std::string inputfile=cms_jets_file);
    // END DEBUG
}
#endif
