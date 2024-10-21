#!/usr/bin/env bash

./write/new_enc/4particle --use_opendata true --use_deltaR --use_pt --weights 1 1 1 --minbin -6 --n_events 100000 --nbins 25 --recursive_phi true --file_prefix od_100k_25bins
