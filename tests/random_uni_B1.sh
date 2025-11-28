#!/bin/sh
FOLDER="/home/wsljan/MP-Workspace/Skyscraper-Invariant/example_files/indecomps_at/random"
PROGRAM="/home/wsljan/MP-Workspace/Skyscraper-Invariant/build/random_uni_B1"

for i in 2 3 4 5 6 7 8 9 10; do
    for j in 1 2 3 4 5; do
        OUTPUT="$FOLDER/random_uni_B1_${i}_$j.scc"
        echo "Generating $OUTPUT"
        "$PROGRAM" "$i" "$OUTPUT"
    done
done

