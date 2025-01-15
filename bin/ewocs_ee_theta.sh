#!/bin/bash

# Start a new tmux session named "my_session"
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux new-session -d -s ewocs_lhc_theta
fi


# Split the window vertically and run the command for twoangle in the new pane
# WTA E Scheme, CA/ca
tmux send-keys './write/ewocs \
    --pair_obs theta --n_events 1000000 \
    --pid_1 11 --pid_2 -11 --energy 1000 --outstate qcd \
    --level parton \
    --fsr true \
    --jet_rad 0.5 --sub_rad 0.0 \
    --jet_alg ee_ca --sub_alg ee_ca \
    --jet_scheme WTA_modp --sub_scheme WTA_modp \
    --contact_terms true \
    --nbins 500 --minbin -6 --maxbin 0 \
    --pt_min 100 --pt_max 100000 \
    --file_prefix ee_qcd_WTA_parton
' C-m
tmux select-layout even-vertical

tmux split-window -v
tmux send-keys './write/ewocs \
    --pair_obs theta --n_events 1000000 \
    --pid_1 11 --pid_2 -11 --energy 1000 --outstate qcd \
    --level hadron \
    --fsr true \
    --jet_rad 0.5 --sub_rad 0.0 \
    --jet_alg ee_ca --sub_alg ee_ca \
    --jet_scheme WTA_modp --sub_scheme WTA_modp \
    --contact_terms true \
    --nbins 500 --minbin -6 --maxbin 0 \
    --pt_min 100 --pt_max 100000 \
    --file_prefix ee_qcd_WTA_hadron
' C-m
tmux select-layout even-vertical



# Split the window vertically and run the command for twoangle in the new pane
# WTA E Scheme, CA/ca
tmux split-window -v
tmux send-keys './write/ewocs \
    --pair_obs theta --n_events 1000000 \
    --pid_1 11 --pid_2 -11 --energy 1000 --outstate qcd \
    --level parton \
    --fsr true \
    --jet_rad 0.5 --sub_rad 0.0 \
    --jet_alg ee_ca --sub_alg ee_ca \
    --jet_scheme WTA_modp --sub_scheme E_scheme \
    --contact_terms true \
    --nbins 500 --minbin -6 --maxbin 0 \
    --pt_min 100 --pt_max 100000 \
    --file_prefix ee_qcd_Escheme_parton
' C-m
tmux select-layout even-vertical

tmux split-window -v
tmux send-keys './write/ewocs \
    --pair_obs theta --n_events 1000000 \
    --pid_1 11 --pid_2 -11 --energy 1000 --outstate qcd \
    --level hadron \
    --fsr true \
    --jet_rad 0.5 --sub_rad 0.0 \
    --jet_alg ee_ca --sub_alg ee_ca \
    --jet_scheme WTA_modp --sub_scheme E_scheme \
    --contact_terms true \
    --nbins 500 --minbin -6 --maxbin 0 \
    --pt_min 100 --pt_max 100000 \
    --file_prefix ee_qcd_Escheme_hadron
' C-m
tmux select-layout even-vertical


# Attach to the session
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux attach-session -t my_session
fi
