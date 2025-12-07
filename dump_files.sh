#!/bin/bash

# Configuration
OUTPUT_FILE="full_directory_dump.txt"
CURRENT_DIR="/global/project/hpcg6049/somatic_pipeline/scripts/_ordered_scripts_20251116_004802/"
# Extensions to ignore (space separated)
IGNORE_EXTS="gz tbi fai fq bed dict"

# Clear or create the output file
> "$OUTPUT_FILE"

echo "Starting dump to $OUTPUT_FILE..."

# ---------------------------------------------------------
# Part 1: Directory Structure
# ---------------------------------------------------------
echo "Writing directory structure..."
echo "### DIRECTORY STRUCTURE ###" >> "$OUTPUT_FILE"

# Check if 'tree' is installed, otherwise use 'find' as a fallback
if command -v tree &> /dev/null; then
    # -a: All files, -I: Ignore .git and the output file
    # We purposely show the full tree structure so you know files exist,
    # even if we skip their content below.
    tree -a -I ".git|node_modules|$OUTPUT_FILE" "$CURRENT_DIR" >> "$OUTPUT_FILE"
else
    find "$CURRENT_DIR" -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g' >> "$OUTPUT_FILE"
fi

echo -e "\n\n### FILE CONTENTS ###\n" >> "$OUTPUT_FILE"

# ---------------------------------------------------------
# Part 2: File Contents
# ---------------------------------------------------------
# Find all regular files (-type f)
# -not -path '*/.*': Ignore hidden files/folders like .git
# -not -name "$OUTPUT_FILE": IMPORTANT - Do not read the file we are writing to
find "$CURRENT_DIR" -type f -not -path '*/.*' -not -name "$OUTPUT_FILE" -print0 | while IFS= read -r -d '' file; do

    # 1. Extension Filtering
    # Check if file ends with any of the ignored extensions
    skip_file=false
    for ext in $IGNORE_EXTS; do
        if [[ "$file" == *."$ext" ]]; then
            echo "Skipping excluded extension (.$ext): $file"
            skip_file=true
            break
        fi
    done

    if [ "$skip_file" = true ]; then
        continue
    fi

    # 2. Binary Filtering
    # Check if the file is a text file (not binary) to avoid printing garbage
    if grep -Iq . "$file"; then
        echo "Adding: $file"
        
        echo "==============================================================================" >> "$OUTPUT_FILE"
        echo "FILE PATH: $file" >> "$OUTPUT_FILE"
        echo "==============================================================================" >> "$OUTPUT_FILE"
        
        cat "$file" >> "$OUTPUT_FILE"
        
        echo -e "\n\n" >> "$OUTPUT_FILE"
    else
        echo "Skipping binary file: $file"
    fi

done

echo "Done! All contents saved to $OUTPUT_FILE"
