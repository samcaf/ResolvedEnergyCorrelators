#ifndef ENC_HEADER
#define ENC_HEADER

#include <string>
#include <string.h>

#include "cmdln.h"
#include "pythia_cmdln.h"


extern const std::string enc_banner;

void write_enc_header(std::string filename,
                           int argc, char* argv[],
                           std::vector<double> weight,
                           bool python_format);

#endif
