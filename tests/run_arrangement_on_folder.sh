#!/bin/sh
FOLDER="${1:-.}"
PROGRAM="/home/wsljan/MP-Workspace/Skyscraper-Invariant/build/arrangement_test"
# Process each file matching the pattern axb.scc (where a,b are integers)
for file in "$FOLDER"/*.scc; do
    "$PROGRAM" "$file" 
done
