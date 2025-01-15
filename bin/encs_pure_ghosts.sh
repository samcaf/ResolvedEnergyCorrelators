#!/bin/bash

# Ghost-dominated PENCs
./write/new_enc/2particle \
    --use_opendata false \
    --use_deltaR --use_pt --weights 1 2 3 \
    --energy 100 --pid_1 2212 --pid_2 2212 --outstate qcd \
    --level parton --isr on --fsr on --mpi off \
    --jet_rad 0.8 --pt_min 0 --pt_max 10000 \
    --minbin -6 --nbins 500 \
    --n_events 10 \
    --add_uniform_ghosts true --mean_ghost_pt 100 \
    --file_prefix pure_ghosts_100k_500bins
