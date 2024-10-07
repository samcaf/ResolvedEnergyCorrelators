#!/bin/bash

# Start a new tmux session named "my_session"
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux new-session -d -s test_w_properties_10k
fi


# Run the command for AKT12 jet properties in the first pane
tmux send-keys './write/jet_properties --use_opendata false --energy 14000 --isr on --fsr on --mpi on --pid_1 2212 --pid_2 2212 --outstate w --jet_rad 1.2 --pt_min 500 --pt_max 550 --minbin -6 --n_events 10000 --nbins 200 --file_prefix w_AKT12_10k_200bins' C-m

# Split the window vertically and run the command for AKT10 in new pane
tmux split-window -v
tmux send-keys './write/jet_properties --use_opendata false --energy 14000 --isr on --fsr on --mpi on --pid_1 2212 --pid_2 2212 --outstate w --jet_rad 1.0 --pt_min 500 --pt_max 550 --minbin -6 --n_events 10000 --nbins 200 --file_prefix w_AKT10_10k_200bins' C-m

# Split the window vertically and run the command for AKT8 in new pane
tmux split-window -v
tmux send-keys './write/jet_properties --use_opendata false --energy 14000 --isr on --fsr on --mpi on --pid_1 2212 --pid_2 2212 --outstate w --jet_rad 0.8 --pt_min 500 --pt_max 550 --minbin -6 --n_events 10000 --nbins 200 --file_prefix w_AKT8_10k_200bins' C-m

# Split the window vertically and run the command for AKT6 in new pane
tmux split-window -v
tmux send-keys './write/jet_properties --use_opendata false --energy 14000 --isr on --fsr on --mpi on --pid_1 2212 --pid_2 2212 --outstate w --jet_rad 0.6 --pt_min 500 --pt_max 550 --minbin -6 --n_events 10000 --nbins 200 --file_prefix w_AKT6_10k_200bins' C-m


# Attach to the session
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux attach-session -t my_session
fi
