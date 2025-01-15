#!/bin/bash

# Start a new tmux session named "my_session"
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux new-session -d -s ewoc_lhc_deltaR
fi


# pt_min = 500
# Parton, Hadron, MPI

# ==============================
# Basic Resolution
# ==============================
# # Parton and Hadron
# tmux send-keys './write/ewocs \
#     --pair_obs deltaR --n_events 1000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level parton \
#     --isr true --fsr true --mpi false \
#     --jet_rad 0.8 --sub_rad 0.0 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --nbins 500 --minbin -6 --maxbin 1 \
#     --pt_min 500 --pt_max 100000 \
#     --file_prefix lhc_w_ptmin500_nomax_parton

# ./write/ewocs \
#     --pair_obs deltaR --n_events 1000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level hadron \
#     --isr true --fsr true --mpi false \
#     --jet_rad 0.8 --sub_rad 0.0 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --nbins 500 --minbin -6 --maxbin 1 \
#     --pt_min 500 --pt_max 100000 \
#     --file_prefix lhc_w_ptmin500_nomax_nompi
# ' C-m
# tmux select-layout even-vertical

# # MPI
# tmux split-window -v
# tmux send-keys './write/ewocs \
#     --pair_obs deltaR --n_events 1000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level hadron \
#     --isr true --fsr true --mpi true \
#     --jet_rad 0.8 --sub_rad 0.0 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --nbins 500 --minbin -6 --maxbin 1 \
#     --pt_min 500 --pt_max 100000 \
#     --file_prefix lhc_w_ptmin500_nomax
# ' C-m
# tmux select-layout even-vertical


# # MPI only, pt_min = 400
# tmux split-window -v
# tmux send-keys './write/ewocs \
#     --pair_obs deltaR --n_events 1000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level hadron \
#     --isr true --fsr true --mpi true \
#     --jet_rad 0.8 --sub_rad 0.0 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --nbins 500 --minbin -6 --maxbin 1 \
#     --pt_min 400 --pt_max 100000 \
#     --file_prefix lhc_w_ptmin400_nomax
# ' C-m
# tmux select-layout even-vertical

# # MPI only, pt_min = 600
# tmux split-window -v
# tmux send-keys './write/ewocs \
#     --pair_obs deltaR --n_events 1000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level hadron \
#     --isr true --fsr true --mpi true \
#     --jet_rad 0.8 --sub_rad 0.0 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --nbins 500 --minbin -6 --maxbin 1 \
#     --pt_min 600 --pt_max 100000 \
#     --file_prefix lhc_w_ptmin600_nomax
# ' C-m
# tmux select-layout even-vertical



# ==============================
# High Resolution
# ==============================
# Parton and Hadron
tmux split-window -v
tmux send-keys './write/ewocs \
    --pair_obs deltaR --n_events 5000000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --level parton \
    --isr true --fsr true --mpi false \
    --jet_rad 0.8 --sub_rad 0.0 \
    --jet_alg akt --sub_alg kt \
    --jet_scheme WTA_pt --sub_scheme WTA_pt \
    --nbins 120 --minbin -0.7 --maxbin -0.3 \
    --pt_min 500 --pt_max 100000 \
    --file_prefix lhc_w_highres120bins_parton

./write/ewocs \
    --pair_obs deltaR --n_events 5000000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --level hadron \
    --isr true --fsr true --mpi false \
    --jet_rad 0.8 --sub_rad 0.0 \
    --jet_alg akt --sub_alg kt \
    --jet_scheme WTA_pt --sub_scheme WTA_pt \
    --nbins 120 --minbin -0.7 --maxbin -0.3 \
    --pt_min 500 --pt_max 100000 \
    --file_prefix lhc_w_highres120bins_nompi
' C-m
tmux select-layout even-vertical

# MPI
tmux split-window -v
tmux send-keys './write/ewocs \
    --pair_obs deltaR --n_events 5000000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --level hadron \
    --isr true --fsr true --mpi true \
    --jet_rad 0.8 --sub_rad 0.0 \
    --jet_alg akt --sub_alg kt \
    --jet_scheme WTA_pt --sub_scheme WTA_pt \
    --nbins 120 --minbin -0.7 --maxbin -0.3 \
    --pt_min 500 --pt_max 100000 \
    --file_prefix lhc_w_highres120bins
' C-m
tmux select-layout even-vertical


# Attach to the session
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux attach-session -t my_session
fi
