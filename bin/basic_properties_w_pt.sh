#!/bin/bash

# Start a new tmux session named "my_session"
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux new-session -d -s ewocs_lhc_mass_pt
fi

# hadrons + MPI
# pt_min = 400
tmux send-keys './write/jet_properties --use_opendata false \
    --n_events 100000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --level hadron \
    --isr true --fsr true --mpi true \
    --jet_rad 0.8 --jet_alg akt --jet_scheme E_scheme \
    --nbins 500 \
    --pt_min 400 --pt_max 100000 --pt_max_bin 3500 \
    --file_prefix lhc_w_ptmin400_nomax_AKT8
' C-m
tmux select-layout even-vertical

# pt_min = 500
tmux split-window -v
tmux send-keys './write/jet_properties --use_opendata false \
    --n_events 100000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --level hadron \
    --isr true --fsr true --mpi true \
    --jet_rad 0.8 --jet_alg akt --jet_scheme E_scheme \
    --nbins 500 \
    --pt_min 500 --pt_max 100000 --pt_max_bin 3500 \
    --file_prefix lhc_w_ptmin500_nomax_AKT8
' C-m
tmux select-layout even-vertical

# pt_min = 600
tmux split-window -v
tmux send-keys './write/jet_properties --use_opendata false \
    --n_events 100000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --level hadron \
    --isr true --fsr true --mpi true \
    --jet_rad 0.8 --jet_alg akt --jet_scheme E_scheme \
    --nbins 500 \
    --pt_min 600 --pt_max 100000 --pt_max_bin 3500 \
    --file_prefix lhc_w_ptmin600_nomax_AKT8
' C-m
tmux select-layout even-vertical


# Attach to the session
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux attach-session -t my_session
fi
