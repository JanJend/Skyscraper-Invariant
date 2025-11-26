#!/bin/sh
FOLDER="${1:-.}"
PROGRAM="/home/wsljan/MP-Workspace/Persistence-Algebra/build/pres_to_quiver"
OUTPUT="$FOLDER/quiver_experiments.md"
echo "# Experiment Results" > "$OUTPUT"
echo "" >> "$OUTPUT"
echo "Generated on: $(date)" >> "$OUTPUT"
echo "" >> "$OUTPUT"
  
for file in "$FOLDER"/*.scc; do
    [ -e "$file" ] || continue
    basename=$(basename "$file")
    
    echo "## $basename" >> "$OUTPUT"
    echo "" >> "$OUTPUT"
    
    for i in 2 3 4 5 6; do
        echo "### i=$i" >> "$OUTPUT"s
        echo "" >> "$OUTPUT"
        start=$(date +%s%N)
        echo '```' >> "$OUTPUT"
        "$PROGRAM" "$file" "$i" >> "$OUTPUT" 2>&1
        echo '```' >> "$OUTPUT"
        end=$(date +%s%N)
        duration=$(awk "BEGIN {printf \"%.3f\", ($end - $start) / 1000000}")
        echo "" >> "$OUTPUT"
        echo "**Execution time:** ${duration}ms" >> "$OUTPUT"
        echo "" >> "$OUTPUT"
    done
done
  
echo "Results saved to: $OUTPUT"