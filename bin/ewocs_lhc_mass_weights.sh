#!/bin/bash

# Start a new tmux session named "my_session"
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux new-session -d -s ewocs_lhc_mass
fi

# -----------------------------------
# Parton, Hadron, and MPI,
# pt_min = 500
# Energy Weight = 2, 3
# -----------------------------------

# =====================================
# Basic resolution
# =====================================
# # Parton, Hadron
# tmux send-keys './write/ewocs \
#     --pair_obs mass --n_events 2000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level parton \
#     --isr true --fsr true --mpi false \
#     --weight 2 \
#     --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --contact_terms false \
#     --nbins 700 --minbin -2 --maxbin 4 \
#     --pt_min 500 --pt_max 100000 \
#     --file_prefix lhc_w_ptmin500_nomax_parton_weight2-0

# ./write/ewocs \
#     --pair_obs mass --n_events 2000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level parton \
#     --isr true --fsr true --mpi false \
#     --weight 3 \
#     --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --contact_terms true \
#     --nbins 700 --minbin -2 --maxbin 4 \
#     --pt_min 500 --pt_max 100000 \
#     --file_prefix lhc_w_ptmin500_nomax_parton_weight3-0

# ./write/ewocs \
#     --pair_obs mass --n_events 2000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level hadron \
#     --isr true --fsr true --mpi false \
#     --weight 2 \
#     --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --contact_terms true \
#     --nbins 700 --minbin -2 --maxbin 4 \
#     --pt_min 500 --pt_max 100000 \
#     --file_prefix lhc_w_ptmin500_nomax_nompi_weight2-0

# ./write/ewocs \
#     --pair_obs mass --n_events 2000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level hadron \
#     --isr true --fsr true --mpi false \
#     --weight 3 \
#     --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --contact_terms true \
#     --nbins 700 --minbin -2 --maxbin 4 \
#     --pt_min 500 --pt_max 100000 \
#     --file_prefix lhc_w_ptmin500_nomax_nompi_weight3-0
# ' C-m
# tmux select-layout even-vertical


# # MPI, Weight = 2
# tmux split-window -v
# tmux send-keys './write/ewocs \
#     --pair_obs mass --n_events 2000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level hadron \
#     --isr true --fsr true --mpi true \
#     --weight 2 \
#     --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --contact_terms true \
#     --nbins 700 --minbin -2 --maxbin 4 \
#     --pt_min 500 --pt_max 100000 \
#     --file_prefix lhc_w_ptmin500_nomax_weight2-0
# ' C-m
# tmux select-layout even-vertical


# # MPI, Weight = 3
# tmux split-window -v
# tmux send-keys './write/ewocs \
#     --pair_obs mass --n_events 2000000 \
#     --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
#     --level hadron \
#     --isr true --fsr true --mpi true \
#     --weight 3 \
#     --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
#     --jet_alg akt --sub_alg kt \
#     --jet_scheme WTA_pt --sub_scheme WTA_pt \
#     --contact_terms true \
#     --nbins 700 --minbin -2 --maxbin 4 \
#     --pt_min 500 --pt_max 100000 \
#     --file_prefix lhc_w_ptmin500_nomax_weight3-0
# ' C-m
# tmux select-layout even-vertical


# =====================================
# High resolution
# =====================================
# Parton, Hadron
tmux split-window -v
tmux send-keys './write/ewocs \
    --pair_obs mass --n_events 5000000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --level parton \
    --isr true --fsr true --mpi false \
    --weight 2 \
    --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
    --jet_alg akt --sub_alg kt \
    --jet_scheme WTA_pt --sub_scheme WTA_pt \
    --contact_terms false \
    --nbins 120 --minbin 1.76 --maxbin 2.00 \
    --pt_min 500 --pt_max 100000 \
    --file_prefix lhc_w_highres120bins_parton_weight2-0

./write/ewocs \
    --pair_obs mass --n_events 5000000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --level parton \
    --isr true --fsr true --mpi false \
    --weight 3 \
    --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
    --jet_alg akt --sub_alg kt \
    --jet_scheme WTA_pt --sub_scheme WTA_pt \
    --contact_terms true \
    --nbins 120 --minbin 1.76 --maxbin 2.00 \
    --pt_min 500 --pt_max 100000 \
    --file_prefix lhc_w_highres120bins_parton_weight3-0

./write/ewocs \
    --pair_obs mass --n_events 5000000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --level hadron \
    --isr true --fsr true --mpi false \
    --weight 2 \
    --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
    --jet_alg akt --sub_alg kt \
    --jet_scheme WTA_pt --sub_scheme WTA_pt \
    --contact_terms true \
    --nbins 120 --minbin 1.76 --maxbin 2.00 \
    --pt_min 500 --pt_max 100000 \
    --file_prefix lhc_w_highres120bins_nompi_weight2-0

./write/ewocs \
    --pair_obs mass --n_events 5000000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --level hadron \
    --isr true --fsr true --mpi false \
    --weight 3 \
    --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
    --jet_alg akt --sub_alg kt \
    --jet_scheme WTA_pt --sub_scheme WTA_pt \
    --contact_terms true \
    --nbins 120 --minbin 1.76 --maxbin 2.00 \
    --pt_min 500 --pt_max 100000 \
    --file_prefix lhc_w_highres120bins_nompi_weight3-0
' C-m
tmux select-layout even-vertical


# MPI, Weight = 2
tmux split-window -v
tmux send-keys './write/ewocs \
    --pair_obs mass --n_events 5000000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --level hadron \
    --isr true --fsr true --mpi true \
    --weight 2 \
    --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
    --jet_alg akt --sub_alg kt \
    --jet_scheme WTA_pt --sub_scheme WTA_pt \
    --contact_terms true \
    --nbins 120 --minbin 1.76 --maxbin 2.00 \
    --pt_min 500 --pt_max 100000 \
    --file_prefix lhc_w_highres120bins_weight2-0
' C-m
tmux select-layout even-vertical


# MPI, Weight = 3
tmux split-window -v
tmux send-keys './write/ewocs \
    --pair_obs mass --n_events 5000000 \
    --pid_1 2212 --pid_2 2212 --energy 14000 --outstate w \
    --level hadron \
    --isr true --fsr true --mpi true \
    --weight 3 \
    --jet_rad 0.8 --sub_rad 0.0 0.03 0.3 0.6 \
    --jet_alg akt --sub_alg kt \
    --jet_scheme WTA_pt --sub_scheme WTA_pt \
    --contact_terms true \
    --nbins 120 --minbin 1.76 --maxbin 2.00 \
    --pt_min 500 --pt_max 100000 \
    --file_prefix lhc_w_highres120bins_weight3-0
' C-m
tmux select-layout even-vertical


# Attach to the session
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux attach-session -t my_session
fi
