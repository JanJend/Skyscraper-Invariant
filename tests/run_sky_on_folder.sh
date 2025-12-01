#!/bin/sh

FOLDER="${1:-.}"  # Use provided folder or current directory
PROGRAM="/home/wsljan/MP-Workspace/Skyscraper-Invariant/build/hnf_main"
OUTPUT="$FOLDER/sky_experiments.md"

# Clear/create the output file
echo "# Experiment Results" > "$OUTPUT"
echo "" >> "$OUTPUT"
echo "Generated on: $(date)" >> "$OUTPUT"
echo "" >> "$OUTPUT"

# Process each file
for file in "$FOLDER"/*.scc "$FOLDER"/*.firep; do
    [ -e "$file" ] || continue
    echo "## $(basename "$file")" >> "$OUTPUT"
    echo "" >> "$OUTPUT"
    
    start=$(date +%s.%N)
    "$PROGRAM" "$file" -p >> "$OUTPUT" 2>&1
    end=$(date +%s.%N)
    duration=$(echo "$end - $start" | bc)
    echo "" >> "$OUTPUT"
    echo "**Execution time:** ${duration}s" >> "$OUTPUT"
done

echo "Results saved to: $OUTPUT"