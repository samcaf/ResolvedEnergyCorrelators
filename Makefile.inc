#===========================================
# Install directory prefixes.
#===========================================
# =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# If you want to install via make
# =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# Versions
PYTHIA_VERSION=8310
FASTJET_VERSION=3.4.2

# Main Directory
SOFTWARE_INSTALL_DIR=/home/samaf/Software
PYTHIA_INSTALL_DIR=$(SOFTWARE_DIR)/pythia$(PYTHIA_VERSION)
FASTJET_INSTALL_DIR=$(SOFTWARE_DIR)/fastjet-$(FASTJET_VERSION)

# =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# Without direct references to dirs
# =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
FJCONFIG=fastjet-config
P8CONFIG=pythia8-config

FASTJET_DIR=$(shell $(FJCONFIG) --prefix)
PYTHIA_DIR=$(shell $(P8CONFIG) --prefix)

# Subdirectories
CONTRIB_INCLUDE=$(FJCONFIG)/include

# LIB_SUFFIX=.dylib
LIB_SUFFIX=.so

# Pythia
PYTHIA_URL=https:\/\/pythia.org\/download\/pythia83\/pythia$(PYTHIA_VERSION).tgz
PYTHIA_BIN=$(PYTHIA_DIR)/bin
PYTHIA_INCLUDE=$(PYTHIA_DIR)/include
PYTHIA_LIB=$(PYTHIA_DIR)/lib
PYTHIA_SHARE=$(PYTHIA_DIR)/share/Pythia8
PYTHIA=$(PYTHIA_LIB)/libpythia8$(LIB_SUFFIX)

# FastJet
FASTJET_URL=http:\/\/fastjet.fr\/repo\/fastjet-$(FASTJET_VERSION).tar.gz
FASTJET_LIB=$(FASTJET_DIR)/src/.libs
FASTJET_INCLUDE=$(FASTJET_DIR)/include
FASTJET=$(FASTJET_LIB)/libfastjet$(LIB_SUFFIX)

# FastJet Contrib
CONTRIB_LIB=$(FASTJET_DIR)/lib
CONTRIB_TOOLS=$(CONTRIB_LIB)/libfastjettools$(LIB_SUFFIX)
CONTRIB=$(CONTRIB_LIB)/libfastjetcontribfragile$(LIB_SUFFIX)


#===========================================
# Compilation flags
#===========================================
# Basic compilation
CXX=g++ -g

CXX_COMMON=-O2 -std=c++17 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_FASTJET:=$(CXX_COMMON) `fastjet-config --libs --cxxflags`
CXX_COMMON:=-I$(PYTHIA_INCLUDE) -I$(CONTRIB_INCLUDE) -I$(FASTJET_INCLUDE) $(CXX_COMMON)
CXX_COMMON+= -L$(PYTHIA_LIB) -L$(CONTRIB_LIB) -L$(FASTJET_LIB) -Wl,-rpath,$(PYTHIA_LIB) -Wl,-rpath,$(CONTRIB_LIB) -Wl,-rpath,$(FASTJET_LIB) -lpythia8 -lfastjet -lfastjetcontribfragile -lfastjettools -ldl
