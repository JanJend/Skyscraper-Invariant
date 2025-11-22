#!/bin/sh

FOLDER="${1:-.}"  # Use provided folder or current directory
PROGRAM="/home/wsljan/MP-Workspace/Persistence-Algebra/build/hnf_main"
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
    
    # Capture start time
    start=$(date +%s)
    
    echo '```' >> "$OUTPUT"
    "$PROGRAM" "$file" "$file" >> "$OUTPUT" 2>&1
    echo '```' >> "$OUTPUT"
    
    # Capture end time and calculate duration
    end=$(date +%s)
    duration=$((end - start))
    
    echo "" >> "$OUTPUT"
    echo "**Execution time:** ${duration}s" >> "$OUTPUT"
    echo "" >> "$OUTPUT"
done

echo "Results saved to: $OUTPUT"