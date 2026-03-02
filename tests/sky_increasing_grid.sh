#!/bin/sh

FILE="${1}"  # Input file (required)
PROGRAM="/home/wsljan/MP-Workspace/Skyscraper-Invariant/build/hnf_main"
OUTPUT="${FILE%.scc}_full_experiments.md"  # Output file based on input filename

# Check if file exists
if [ ! -f "$FILE" ]; then
    echo "Error: File '$FILE' not found"
    exit 1
fi

# Grid sizes to test
GRID_SIZES="20 40 80 160"

# Clear/create the output file
echo "# Experiment Results for $(basename "$FILE")" > "$OUTPUT"
echo "" >> "$OUTPUT"
echo "Generated on: $(date)" >> "$OUTPUT"
echo "" >> "$OUTPUT"

# Process each grid size
for size in $GRID_SIZES; do
    echo "## Grid Size: ${size}x${size}" >> "$OUTPUT"
    echo "" >> "$OUTPUT"
    start=$(date +%s.%N)
    timeout 60m "$PROGRAM" "$FILE" -p -y -r${size},${size} >> "$OUTPUT" 2>&1
    exit_code=$?
    end=$(date +%s.%N)
    duration=$(echo "$end - $start" | bc)
    echo "" >> "$OUTPUT"
    if [ $exit_code -eq 124 ]; then
        echo "**Status:** TIMEOUT (60 minutes exceeded)" >> "$OUTPUT"
    fi
    echo "**Execution time:** ${duration}s" >> "$OUTPUT"
    echo "" >> "$OUTPUT"
done

echo "Results saved to: $OUTPUT"