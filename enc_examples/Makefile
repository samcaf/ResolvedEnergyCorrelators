# Makefile for Energy Correlator Example

CXX = g++
CXXFLAGS = -Wall -O2 -std=c++11

# Add FastJet and Pythia flags
# You may need to adjust these paths for your system
FASTJET_CONFIG = fastjet-config
FASTJET_CXXFLAGS = $(shell $(FASTJET_CONFIG) --cxxflags)
FASTJET_LIBS = $(shell $(FASTJET_CONFIG) --libs)

# Pythia configuration - adjust paths as needed
PYTHIA8_DIR = $(shell pythia8-config --prefix)
PYTHIA8_INCLUDE = -I$(PYTHIA8_DIR)/include
PYTHIA8_LIBS = -L$(PYTHIA8_DIR)/lib -lpythia8

# Include directories
INCLUDE_DIRS = -Iinclude $(PYTHIA8_INCLUDE) $(FASTJET_CXXFLAGS)

# Libraries
LIBS = $(FASTJET_LIBS) $(PYTHIA8_LIBS)

# Target executable
TARGET = example_enc_calculator

# Source files
SOURCES = example_enc_calculator.cc

# Object files
OBJECTS = $(SOURCES:.cc=.o)

# Main rule
all: $(TARGET)

# Rule to make the executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# Rule to make object files
%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean
