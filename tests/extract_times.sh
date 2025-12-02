#!/bin/bash

input_file="$1"
output_file="${input_file%.md}_extracted.txt"

awk '
/^## / {
    sub(/^## /, "", $0)
    filename = $0
    next
}
/\*\*Execution time:/ {
    sub(/\*\*Execution time:\*\* /, "", $0)
    print filename " | " $0
}
' "$input_file" > "$output_file"

echo "Extracted to: $output_file"