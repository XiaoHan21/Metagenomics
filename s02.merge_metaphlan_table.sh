#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_metaphlan_dir> <output_tsv_dir>"
    echo "Example: $0 /path/to/05.metaphlan /path/to/output_dir"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2

# Check if the input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory '$INPUT_DIR' does not exist."
    exit 1
fi

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

# Activate the conda environment
echo "Activating conda environment (humann4_mpa4.1.1)..."
source ~/xiaohan/Software/miniconda3/bin/activate
conda activate humann4_mpa4.1.1

# Execute the merging and formatting operations
echo "Merging MetaPhlAn tables..."

# Note: The original path matching pattern /*/*_metaphlan_profile.tsv is preserved here
merge_metaphlan_tables.py "${INPUT_DIR}"/*/*_metaphlan_profile.tsv \
    | sed 's/_1_metaphlan//g' \
    | tail -n+2 \
    | sed '1 s/clade_name/ID/' \
    | sed '2i #metaphlan4' \
    > "$OUTPUT_FILE"

# Check if the output file was generated successfully
if [ $? -eq 0 ]; then
    echo "Success! Merged table saved to $OUTPUT_FILE"
else
    echo "Error: Failed to merge tables."
    exit 1
fi
