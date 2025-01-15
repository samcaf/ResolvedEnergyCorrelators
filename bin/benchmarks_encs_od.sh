# Start a new tmux session named "my_session"
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux new-session -d -s new_enc_od_100k
fi

# Run the command for oneangle in the first pane
tmux send-keys './write/benchmark/new_enc_2particle --use_opendata true --use_deltaR --use_pt --minbin -6 --n_events 100000 --nbins 200' C-m
tmux send-keys './write/benchmark/new_enc_3particle --use_opendata true --use_deltaR --use_pt --minbin -6 --n_events 100000 --nbins 150' C-m
tmux send-keys './write/benchmark/new_enc_4particle --use_opendata true --use_deltaR --use_pt --minbin -6 --n_events 100000 --nbins 50' C-m

tmux select-layout even-vertical

# Attach to the session
if ! [ "$TERM_PROGRAM" = tmux ]; then
tmux attach-session -t my_session
fi
