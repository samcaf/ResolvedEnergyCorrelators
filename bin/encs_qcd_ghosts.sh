#!/bin/bash

# Start a new tmux session named "my_session"
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux new-session -d -s new_enc_qcd_100k
fi


# PENC
tmux send-keys '
./write/new_enc/2particle \
    --use_opendata false \
    --use_deltaR --use_pt --weights 1 2 3 \
    --energy 14000 --pid_1 2212 --pid_2 2212 --outstate qcd \
    --level parton --isr on --fsr on --mpi off \
    --jet_rad 0.8 --pt_min 500 --pt_max 550 \
    --minbin -6 --nbins 500 \
    --n_events 100000 \
    --add_uniform_ghosts true --mean_ghost_pt 0.1 \
    --file_prefix qcd_100M_ghosts_parton_100k_500bins

./write/new_enc/2particle \
    --use_opendata false \
    --use_deltaR --use_pt --weights 1 2 3 \
    --energy 14000 --pid_1 2212 --pid_2 2212 --outstate qcd \
    --level hadron --isr on --fsr on --mpi off \
    --jet_rad 0.8 --pt_min 500 --pt_max 550 \
    --minbin -6 --nbins 500 \
    --n_events 100000 \
    --add_uniform_ghosts true --mean_ghost_pt 0.1 \
    --file_prefix qcd_100M_ghosts_nompi_100k_500bins
' C-m

# RE3C
# Parton
tmux split-window -v
tmux send-keys '
./write/new_enc/3particle \
    --use_opendata false \
    --use_deltaR --use_pt --weights 1 1 2 2 3 3 \
    --energy 14000 --pid_1 2212 --pid_2 2212 --outstate qcd \
    --level parton --isr on --fsr on --mpi off \
    --jet_rad 0.8 --pt_min 500 --pt_max 550 \
    --minbin -6 --nbins 150 \
    --n_events 100000 \
    --add_uniform_ghosts true --mean_ghost_pt 0.1 \
    --file_prefix qcd_100M_ghosts_parton_100k_150bins
'
# Hadron
tmux split-window -v
tmux send-keys '
./write/new_enc/3particle \
    --use_opendata false \
    --use_deltaR --use_pt --weights 1 1 2 2 3 3 \
    --energy 14000 --pid_1 2212 --pid_2 2212 --outstate qcd \
    --level hadron --isr on --fsr on --mpi off \
    --jet_rad 0.8 --pt_min 500 --pt_max 550 \
    --minbin -6 --nbins 150 \
    --n_events 100000 \
    --add_uniform_ghosts true --mean_ghost_pt 0.1 \
    --file_prefix qcd_100M_ghosts_nompi_100k_150bins
' C-m
tmux select-layout even-vertical

# Old E3C
# Parton
tmux split-window -v
tmux send-keys '
./write/new_enc/old_3particle \
    --use_opendata false \
    --use_deltaR --use_pt --weights 1 1 \
    --energy 14000 --pid_1 2212 --pid_2 2212 --outstate qcd \
    --level parton --isr on --fsr on --mpi off \
    --jet_rad 0.8 --pt_min 500 --pt_max 550 \
    --minbin -6 --nbins 150 \
    --n_events 100000 \
    --add_uniform_ghosts true --mean_ghost_pt 0.1 \
    --file_prefix qcd_100M_ghosts_parton_100k_150bins
' C-m
tmux select-layout even-vertical
# Hadron
tmux split-window -v
tmux send-keys '
./write/new_enc/old_3particle \
    --use_opendata false \
    --use_deltaR --use_pt --weights 1 1 \
    --energy 14000 --pid_1 2212 --pid_2 2212 --outstate qcd \
    --level hadron --isr on --fsr on --mpi off \
    --jet_rad 0.8 --pt_min 500 --pt_max 550 \
    --minbin -6 --nbins 150 \
    --n_events 100000 \
    --add_uniform_ghosts true --mean_ghost_pt 0.1 \
    --file_prefix qcd_100M_ghosts_nompi_100k_150bins
' C-m
tmux select-layout even-vertical

# Attach to the session
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux attach-session -t my_session
fi
