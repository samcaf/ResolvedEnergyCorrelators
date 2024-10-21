#!/usr/bin/env bash

./write/new_enc/4particle --use_opendata false --use_deltaR --use_pt --weights 1 1 1 --energy 14000 --isr on --fsr on --mpi on --pid_1 2212 --pid_2 2212 --outstate top --jet_rad 0.8 --pt_min 500 --pt_max 550 --minbin -6 --n_events 1000000 --nbins 25 --recursive_phi true --file_prefix top_1M_25bins
