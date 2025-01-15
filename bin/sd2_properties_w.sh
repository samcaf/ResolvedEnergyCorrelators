#!/bin/bash

# Start a new tmux session named "my_session"
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux new-session -d -s ewocs_lhc_mass_pt
fi

# SD, beta = 2
# pt_min = 500

# ==============================
# Basic Resolution
# ==============================
# parton, hadron w/o MPI
tmux send-keys './write/groomed_properties --use_opendata false \
    --n_events 1000000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --isr true --fsr true \
    --level parton --mpi false \
    --jet_rad 0.8 --jet_alg akt --jet_scheme E_scheme \
    --beta_SD 2 --zcut 0.1 --R_SD 0.8 \
    --nbins 250 \
    --pt_min 500 --pt_max 100000 --pt_max_bin 3500 \
    --file_prefix lhc_w_parton_ptmin500_AKT8_Escheme_sd2

./write/groomed_properties --use_opendata false \
    --n_events 1000000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --isr true --fsr true \
    --level hadron --mpi false \
    --jet_rad 0.8 --jet_alg akt --jet_scheme E_scheme \
    --beta_SD 2 --zcut 0.1 --R_SD 0.8 \
    --nbins 250 \
    --pt_min 500 --pt_max 100000 --pt_max_bin 3500 \
    --file_prefix lhc_w_nompi_ptmin500_AKT8_Escheme_sd2
' C-m
tmux select-layout even-vertical

# hadron + MPI
tmux split-window -v
tmux send-keys './write/groomed_properties --use_opendata false \
    --n_events 1000000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --isr true --fsr true \
    --level hadron --mpi true \
    --jet_rad 0.8 --jet_alg akt --jet_scheme E_scheme \
    --beta_SD 2 --zcut 0.1 --R_SD 0.8 \
    --nbins 250 \
    --pt_min 500 --pt_max 100000 --pt_max_bin 3500 \
    --file_prefix lhc_w_ptmin500_AKT8_Escheme_sd2
' C-m
tmux select-layout even-vertical


# Attach to the session
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux attach-session -t my_session
fi
