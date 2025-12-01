#!/bin/sh
FOLDER="${1:-.}"
PROGRAM="/home/wsljan/MP-Workspace/Skyscraper-Invariant/build/hnf_at_origin"
OUTPUT="$FOLDER/origin_opt_experiments.md"
echo "# Experiment Results" > "$OUTPUT"
echo "" >> "$OUTPUT"
echo "Generated on: $(date)" >> "$OUTPUT"
echo "" >> "$OUTPUT"
# Process each file matching the pattern axb.scc (where a,b are integers)
for file in "$FOLDER"/*.scc; do
    [ -e "$file" ] || continue
    basename=$(basename "$file")
    echo "## $basename" >> "$OUTPUT"
    echo "" >> "$OUTPUT"
    start=$(date +%s%N)
    echo '```' >> "$OUTPUT"
    "$PROGRAM" "$file" "1.6,0.1" "true" >> "$OUTPUT" 2>&1
    echo '```' >> "$OUTPUT"
    end=$(date +%s%N)
    duration=$(awk "BEGIN {printf \"%.3f\", ($end - $start) / 1000000}")
    echo "" >> "$OUTPUT"
    echo "**Execution time:** ${duration}ms" >> "$OUTPUT"
    echo "" >> "$OUTPUT"
done
echo "Results saved to: $OUTPUT"