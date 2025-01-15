#!/bin/bash

# Start a new tmux session named "my_session"
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux new-session -d -s ewoc_lhc_deltaR
fi


# MPI,  pt_min = 500

# ==============================
# Basic Resolution
# ==============================
# # Charged only
# tmux send-keys './write/ewocs \
#     --pair_obs deltaR --n_events 1000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level hadron \
#     --isr true --fsr true --mpi true \
#     --jet_rad 0.8 --sub_rad 0.0 0.01 0.03 0.1 0.3 0.6 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --nbins 500 --minbin -6 --maxbin 1 \
#     --pt_min 500 --pt_max 100000 \
#     --charged_only \
#     --file_prefix lhc_w_ptmin500_nomax_charged
# ' C-m
# tmux select-layout even-vertical

# # Smearing
# tmux split-window -v
# tmux send-keys './write/ewocs \
#     --pair_obs deltaR --n_events 100000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level hadron \
#     --isr true --fsr true --mpi true \
#     --jet_rad 0.8 --sub_rad 0.0 0.01 0.03 0.1 0.3 0.6 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --nbins 500 --minbin -6 --maxbin 1 \
#     --pt_min 500 --pt_max 100000 \
#     --smear_factor 3 \
#     --file_prefix lhc_w_ptmin500_nomax_threesmeared
# ' C-m
# tmux select-layout even-vertical


# ==============================
# High Resolution
# ==============================
# # Charged only
# tmux split-window -v
# tmux send-keys './write/ewocs \
#     --pair_obs deltaR --n_events 5000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level hadron \
#     --isr true --fsr true --mpi true \
#     --jet_rad 0.8 --sub_rad 0.0 0.01 0.03 0.1 0.3 0.6 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --nbins 120 --minbin -0.7 --maxbin -0.3 \
#     --pt_min 500 --pt_max 100000 \
#     --charged_only \
#     --file_prefix lhc_w_highres120bins_charged
# ' C-m
# tmux select-layout even-vertical

# Smearing
tmux split-window -v
tmux send-keys './write/ewocs \
    --pair_obs deltaR --n_events 5000000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --level hadron \
    --isr true --fsr true --mpi true \
    --jet_rad 0.8 --sub_rad 0.0 0.01 0.03 0.1 0.3 0.6 \
    --jet_alg akt --sub_alg kt \
    --jet_scheme WTA_pt --sub_scheme WTA_pt \
    --nbins 120 --minbin -0.7 --maxbin -0.3 \
    --pt_min 500 --pt_max 100000 \
    --smear_factor 1 \
    --file_prefix lhc_w_highres120bins_onesmeared
' C-m
tmux select-layout even-vertical

# Attach to the session
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux attach-session -t my_session
fi
