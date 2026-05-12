#!/bin/bash
# Enable pipefail to ensure the pipeline fails if any command within it fails
set -o pipefail

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_metaphlan_dir> <output_tsv_dir>"
    echo "Example: $0 ../03.result/05.metaphlan ../03.result/05.metaphlan"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2

# Create the output directory if it does not exist
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
    echo "Created output directory: $OUTPUT_DIR"
fi

OUTPUT_FILE="${OUTPUT_DIR}/taxonomy.tsv"

echo "========================================="
echo "Input directory : $INPUT_DIR"
echo "Output file     : $OUTPUT_FILE"
echo "========================================="

# Activate the conda environment properly to avoid ArgumentError
echo "Activating conda environment (humann4_mpa4.1.1)..."
source ~/xiaohan/Software/miniconda3/bin/activate
conda activate humann4_mpa4.1.1

# Locate MetaPhlAn bugs list tables using the correct file extension
echo "Locating MetaPhlAn bugs list tables..."
FILES_TO_MERGE=$(find "$INPUT_DIR" -type f -name "*_metaphlan_bugs_list.tsv")

# Check if files were actually found
if [ -z "$FILES_TO_MERGE" ]; then
    echo "Error: No files matching '*_metaphlan_bugs_list.tsv' were found in $INPUT_DIR or its subdirectories."
    echo "Please double-check your input directory."
    exit 1
fi

# Execute the merging and formatting operations
echo "Merging MetaPhlAn tables..."

# Merge tables and clean up headers
# Changed 's/_1_metaphlan//g' to 's/_metaphlan_bugs_list//g' to match the actual filenames
merge_metaphlan_tables.py $FILES_TO_MERGE \
    | sed 's/_metaphlan_bugs_list//g' \
    | tail -n+2 \
    | sed '1 s/clade_name/ID/' \
    | sed '2i #metaphlan4' \
    > "$OUTPUT_FILE"

# Check if the pipeline executed successfully
if [ $? -eq 0 ]; then
    echo "Success! Merged table saved to $OUTPUT_FILE"
else
    echo "Error: Failed to merge tables."
    exit 1
fi
