#!/bin/bash

# Start a new tmux session named "my_session"
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux new-session -d -s ewocs_lhc_mass
fi


# Parton, Hadron, MPI, pt_min = 500

# ==============================
# Basic Resolution
# ==============================
# Parton and Hadron
# tmux send-keys './write/ewocs \
#     --pair_obs mass --n_events 2000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level parton \
#     --isr true --fsr true --mpi false \
#     --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --contact_terms false \
#     --nbins 700 --minbin -2 --maxbin 4 \
#     --pt_min 500 --pt_max 100000 \
#     --file_prefix lhc_w_ptmin500_nomax_parton
# '

# ./write/ewocs \
#     --pair_obs mass --n_events 2000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level hadron \
#     --isr true --fsr true --mpi false \
#     --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --contact_terms true \
#     --nbins 700 --minbin -2 --maxbin 4 \
#     --pt_min 500 --pt_max 100000 \
#     --file_prefix lhc_w_ptmin500_nomax_nompi
# ' C-m
# tmux select-layout even-vertical

# # MPI
# tmux split-window -v
# tmux send-keys './write/ewocs \
#     --pair_obs mass --n_events 2000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level hadron \
#     --isr true --fsr true --mpi true \
#     --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --contact_terms true \
#     --nbins 700 --minbin -2 --maxbin 4 \
#     --pt_min 500 --pt_max 100000 \
#     --file_prefix lhc_w_ptmin500_nomax
# ' C-m
# tmux select-layout even-vertical

# : "
# # MPI only, pt_min = 400
# tmux split-window -v
# tmux send-keys './write/ewocs \
#     --pair_obs mass --n_events 2000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level hadron \
#     --isr true --fsr true --mpi true \
#     --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --contact_terms true \
#     --nbins 700 --minbin -2 --maxbin 4 \
#     --pt_min 400 --pt_max 100000 \
#     --file_prefix lhc_w_ptmin400_nomax
# ' C-m
# tmux select-layout even-vertical


# # MPI only, pt_min = 600
# tmux split-window -v
# tmux send-keys './write/ewocs \
#     --pair_obs mass --n_events 2000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level hadron \
#     --isr true --fsr true --mpi true \
#     --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --contact_terms true \
#     --nbins 700 --minbin -2 --maxbin 4 \
#     --pt_min 600 --pt_max 100000 \
#     --file_prefix lhc_w_ptmin600_nomax
# ' C-m
# tmux select-layout even-vertical
# "


# ==============================
# High Resolution
# ==============================
# Parton and Hadron
tmux split-window -v
tmux send-keys './write/ewocs \
    --pair_obs mass --n_events 5000000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --level parton \
    --isr true --fsr true --mpi false \
    --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
    --jet_alg akt --sub_alg kt \
    --jet_scheme WTA_pt --sub_scheme WTA_pt \
    --contact_terms false \
    --nbins 120 --minbin 1.76 --maxbin 2.00 \
    --pt_min 500 --pt_max 100000 \
    --file_prefix lhc_w_highres120bins_parton

./write/ewocs \
    --pair_obs mass --n_events 5000000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --level hadron \
    --isr true --fsr true --mpi false \
    --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
    --jet_alg akt --sub_alg kt \
    --jet_scheme WTA_pt --sub_scheme WTA_pt \
    --contact_terms true \
    --nbins 120 --minbin 1.76 --maxbin 2.00 \
    --pt_min 500 --pt_max 100000 \
    --file_prefix lhc_w_highres120bins_nompi
' C-m
tmux select-layout even-vertical

# MPI
tmux split-window -v
tmux send-keys './write/ewocs \
    --pair_obs mass --n_events 5000000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --level hadron \
    --isr true --fsr true --mpi true \
    --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
    --jet_alg akt --sub_alg kt \
    --jet_scheme WTA_pt --sub_scheme WTA_pt \
    --contact_terms true \
    --nbins 120 --minbin 1.76 --maxbin 2.00 \
    --pt_min 500 --pt_max 100000 \
    --file_prefix lhc_w_highres120bins
# ' C-m
# tmux select-layout even-vertical


# Attach to the session
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux attach-session -t my_session
fi
