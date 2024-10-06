# Start a new tmux session named "my_session"
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux new-session -d -s new_enc_od_100k
fi

# Run the command for oneangle in the first pane
tmux send-keys './write/new_enc/2particle --use_opendata true --use_deltaR --use_pt --weights 0.5 1 2 2.5 2.71828 3 3.14159 3.5 4 --minbin -6 --n_events 100000 --nbins 200 --file_prefix od_100k_200bins' C-m

# Split the window vertically and run the command for twoangle in the new pane
tmux split-window -v
vmux send-keys './write/new_enc/3particle --use_opendata true --use_deltaR --use_pt --weights 1 1 --minbin -6 --n_events 100000 --nbins 150 --file_prefix od_100k_150bins' C-m

# Split the window vertically and run the command for twoangle in the new pane
tmux split-window -v
tmux send-keys './write/new_enc/old_3particle --use_opendata true --use_deltaR --use_pt --weights 1 1 --minbin -6 --n_events 100000 --nbins 150 --file_prefix od_100k_150bins' C-m

# # Split the window vertically and run the command for threeangle in the new pane
tmux split-window -v
tmux send-keys './write/new_enc/4particle --use_opendata true --use_deltaR --use_pt --weights 1 1 1 --minbin -6 --n_events 100000 --nbins 20 --recursive_phi true --file_prefix od_rec_100k_20bins' C-m

# Attach to the session
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux attach-session -t my_session
fi
