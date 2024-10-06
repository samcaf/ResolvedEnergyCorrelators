#!/bin/bash

# Start a new tmux session named "my_session"
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux new-session -d -s new_enc_w_1M
fi

# Run the command for 2particle in the first pane
tmux send-keys './write/new_enc/2particle --use_opendata false --use_deltaR --use_pt --weights 0.5 1 1.5 2 2.5 2.71828 3 3.14159 3.5 4 --energy 14000 --isr on --fsr on --mpi on --pid_1 2212 --pid_2 2212 --outstate w --jet_rad 1.0 --pt_min 500 --pt_max 550 --minbin -6 --n_events 1000000 --nbins 500 --file_prefix w_1M_500bins' C-m

# Split the window vertically and run the command for 3particle in the new pane
tmux split-window -v
tmux send-keys './write/new_enc/3particle --use_opendata false --use_deltaR --use_pt --weights 1 1 --energy 14000 --isr on --fsr on --mpi on --pid_1 2212 --pid_2 2212 --outstate w --jet_rad 1.0 --pt_min 500 --pt_max 550 --minbin -6 --n_events 1000000 --nbins 150 --file_prefix w_1M_150bins' C-m

# Split the window vertically and run the command for 3particle in the new pane
tmux split-window -v
tmux send-keys './write/new_enc/old_3particle --use_opendata false --use_deltaR --use_pt --weights 1 1 --energy 14000 --isr on --fsr on --mpi on --pid_1 2212 --pid_2 2212 --outstate w --jet_rad 1.0 --pt_min 500 --pt_max 550 --minbin -6 --n_events 1000000 --nbins 150 --file_prefix w_1M_150bins' C-m

# Split the window vertically and run the command for 4particle in the new pane
tmux split-window -v
tmux send-keys './write/new_enc/4particle --use_opendata false --use_deltaR --use_pt --weights 1 1 1 --energy 14000 --isr on --fsr on --mpi on --pid_1 2212 --pid_2 2212 --outstate w --jet_rad 1.0 --pt_min 500 --pt_max 550 --minbin -6 --n_events 1000000 --nbins 20 --recursive_phi true --file_prefix w_rec_1M_20bins' C-m


# Attach to the session
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux attach-session -t my_session
fi
