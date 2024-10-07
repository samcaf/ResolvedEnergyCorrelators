/**
 * @file    pythia_cmdln.h
 *
 * @brief   A command line utility header file for pythia event generation.
 */
#ifndef PYTHIA_CMDLN
#define PYTHIA_CMDLN

#include <string>
#include <string.h>
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>

#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"

// Local imports:
#include "general_utils.h"
#include "jet_utils.h"
#include "cmdln.h"


// =====================================
// Command Line Basics
// =====================================
// Command Line Defaults
extern const int                           _NEVENTS_DEFAULT;
extern const double                        _ENERGY_DEFAULT;
extern const std::string                   _LEVEL_DEFAULT;
extern const std::string                   _OUTSTATE_DEFAULT;
extern const int                           _PID_1_DEFAULT;
extern const int                           _PID_2_DEFAULT;

int showermodel_cmdln(int argc, char* argv[]);


// =====================================
// Setup Utilities
// =====================================
// Verification of paramters
int checkPythiaInputs(int argc, char* argv[]);

// Setting up pythia
void setup_pythia_cmdln(Pythia8::Pythia &pythia, int argc, char* argv[]);

void write_jetproperty_header(std::string filename,
                              int argc, char* argv[],
                              bool python_format);
#endif
