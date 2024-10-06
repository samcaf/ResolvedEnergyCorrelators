#===========================================
# Install directory prefixes.
#===========================================
# Versions
PYTHIA_VERSION=8310
FASTJET_VERSION=3.4.2

# Main Directory
SOFTWARE_DIR=/home/samaf/Software

# Subdirectories
PYTHIA_DIR=$(SOFTWARE_DIR)/pythia$(PYTHIA_VERSION)
FASTJET_DIR=$(SOFTWARE_DIR)/fastjet-$(FASTJET_VERSION)
CONTRIB_DIR=$(SOFTWARE_DIR)/fastjet-install

# LIB_SUFFIX=.dylib
LIB_SUFFIX=.so

# Pythia
PYTHIA_URL=https:\/\/pythia.org\/download\/pythia83\/pythia$(PYTHIA_VERSION).tgz
PYTHIA_BIN=$(PYTHIA_DIR)/bin
PYTHIA_INCLUDE=$(PYTHIA_DIR)/include
PYTHIA_LIB=$(PYTHIA_DIR)/lib
PYTHIA_SHARE=$(PYTHIA_DIR)/share/Pythia8
PYTHIA=$(PYTHIA_LIB)/libpythia8$(LIB_SUFFIX)

# Fastjet
FASTJET_URL=http:\/\/fastjet.fr\/repo\/fastjet-$(FASTJET_VERSION).tar.gz
FASTJET_LIB=$(FASTJET_DIR)/src/.libs
FASTJET_INCLUDE=$(FASTJET_DIR)/include
FASTJET=$(FASTJET_LIB)/libfastjet$(LIB_SUFFIX)

# fastjet-contrib
CONTRIB_LIB=$(CONTRIB_DIR)/lib
CONTRIB_INCLUDE=$(CONTRIB_DIR)/include


#===========================================
# Compilation flags
#===========================================
# Basic compilation
CXX=g++ -g

# CXX_COMMON=-O2 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_COMMON=-O2 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_COMMON:=-I$(PYTHIA_INCLUDE) -I$(FASTJET_INCLUDE) $(CXX_COMMON) $(GZIP_LIB)
CXX_COMMON+= -L$(PYTHIA_LIB) -L$(FASTJET_LIB) -Wl,-rpath,$(PYTHIA_LIB) -Wl,-rpath,$(FASTJET_LIB) -lpythia8 -lfastjet -ldl

CXX_SHARED=-dynamiclib
CXX_SONAME=-Wl,-dylib_install_name,@rpath/
