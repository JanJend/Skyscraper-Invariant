#!/bin/sh
  
FOLDER="${1:-.}"
PROGRAM="/home/wsljan/MP-Workspace/Persistence-Algebra/build/hnf_at_origin"
OUTPUT="$FOLDER/origin_experiments.md"
  
echo "# Experiment Results" > "$OUTPUT"
echo "" >> "$OUTPUT"
echo "Generated on: $(date)" >> "$OUTPUT"
echo "" >> "$OUTPUT"
  
# Process each file matching the pattern axb.scc (where a,b are integers)
for file in "$FOLDER"/*x*.scc; do
    [ -e "$file" ] || continue
    
    basename=$(basename "$file")
    case "$basename" in
        [0-9]*x[0-9]*.scc)
            ;;
        *)
            continue
            ;;
    esac
    
    echo "## $basename" >> "$OUTPUT"
    echo "" >> "$OUTPUT"
    start=$(date +%s%N)
    echo '```' >> "$OUTPUT"
    "$PROGRAM" "$file" "$file" >> "$OUTPUT" 2>&1
    echo '```' >> "$OUTPUT"
    end=$(date +%s%N)
    duration=$(( (end - start) / 1000000 ))
    echo "" >> "$OUTPUT"
    echo "**Execution time:** ${duration}ms" >> "$OUTPUT"
    echo "" >> "$OUTPUT"
done
  
echo "Results saved to: $OUTPUT"