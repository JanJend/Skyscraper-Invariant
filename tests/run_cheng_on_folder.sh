#!/bin/sh
FOLDER="${1:-.}"
VENV="/home/wsljan/myenv"
PROGRAM="/home/wsljan/MP-Workspace/sky-inv-quiv/scripts/cheng_alg.py"
OUTPUT="$FOLDER/cheng_alg_Q_experiments.md"

echo "# Experiment Results" > "$OUTPUT"
echo "" >> "$OUTPUT"
echo "Generated on: $(date)" >> "$OUTPUT"
echo "" >> "$OUTPUT"


# Process each file matching the pattern axb.scc (where a,b are integers)
for file in "$FOLDER"/*x*.txt; do
    [ -e "$file" ] || continue
    basename=$(basename "$file")
    case "$basename" in
        [0-9]*x[0-9]*.txt)
            ;;
        *)
            continue
            ;;
    esac
    
    # Skip 32x32 and 16x16 files
    case "$basename" in
        *32x32*.txt|*16x16*.txt)
            echo "Skipping $basename (32x32 or 16x16)"
            continue
            ;;
    esac

    echo "## $basename" >> "$OUTPUT"
    echo "" >> "$OUTPUT"
    
    start=$(date +%s%N)
    echo '```' >> "$OUTPUT"
    timeout 600s "$VENV/bin/python" "$PROGRAM" "$file" "Q" >> "$OUTPUT" 2>&1
    exit_code=$?
    echo '```' >> "$OUTPUT"
    end=$(date +%s%N)

    if [ $exit_code -eq 124 ]; then
        echo "" >> "$OUTPUT"
        echo "**Status:** TIMEOUT (exceeded 600 seconds)" >> "$OUTPUT"
        echo "" >> "$OUTPUT"
    else
        duration=$(awk "BEGIN {printf \"%.3f\", ($end - $start) / 1000000}")
        echo "" >> "$OUTPUT"
        echo "**Execution time:** ${duration}ms" >> "$OUTPUT"
        echo "" >> "$OUTPUT"
fi
done

echo "Results saved to: $OUTPUT"