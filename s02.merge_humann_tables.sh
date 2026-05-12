#!/bin/bash
# Enable pipefail to ensure the pipeline fails if any command within it fails
set -o pipefail

# Support either 1 or 2 arguments. If 1 argument is provided, output dir = input dir.
if [ "$#" -eq 1 ]; then
    INPUT_DIR=$1
    OUTPUT_DIR=$1
elif [ "$#" -eq 2 ]; then
    INPUT_DIR=$1
    OUTPUT_DIR=$2
else
    echo "Usage: $0 <humann_dir> [output_tsv_dir]"
    echo "Example: $0 ./03.result/06.humann"
    exit 1
fi

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

echo "========================================="
echo "Input directory : $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "========================================="

# Activate the conda environment properly
echo "Activating conda environment (humann4_mpa4.1.1)..."
source ~/xiaohan/Software/miniconda3/bin/activate
conda activate humann4_mpa4.1.1

# Define the standard HUMAnN output table types
TABLE_TYPES=("pathabundance" "genefamilies" "reactions")

# 1. Create a tmp directory inside the output directory
TMP_DIR="${OUTPUT_DIR}/tmp_merge"
mkdir -p "$TMP_DIR"
echo "Created temporary directory: $TMP_DIR"

# 2. Copy files to the tmp directory (Safely ignoring the tmp directory itself)
echo "Copying files from subdirectories to the temporary directory..."
for type in "${TABLE_TYPES[@]}"; do
    # The '-path "$TMP_DIR" -prune' part ensures find skips the tmp directory
    find "$INPUT_DIR" -path "$TMP_DIR" -prune -o -type f -name "*_${type}.tsv" -exec cp {} "$TMP_DIR/" \;
done

# Verify that the temporary directory actually contains files
if [ -z "$(ls -A "$TMP_DIR")" ]; then
    echo "Error: No HUMAnN .tsv files were found to copy."
    echo "Cleaning up and exiting..."
    rm -rf "$TMP_DIR"
    exit 1
fi

echo "Starting to merge and process HUMAnN tables..."

# 3. Loop through table types to merge, normalize, and split
for type in "${TABLE_TYPES[@]}"; do
    echo "-----------------------------------------"
    echo "Processing: $type"
    OUTPUT_FILE="${OUTPUT_DIR}/${type}.tsv"
    RELAB_FILE="${OUTPUT_DIR}/${type}_relab.tsv"

    # Step A: Merge tables from the flat temporary directory
    echo "  -> Merging tables..."
    humann_join_tables \
        --input "$TMP_DIR" \
        --file_name "$type" \
        --output "$OUTPUT_FILE"

    if [ $? -eq 0 ]; then
        # Format sampleID by removing the '_Abundance' suffix
        sed -i 's/_Abundance//g' "$OUTPUT_FILE"
    else
        echo "Error: Failed to merge $type tables."
        rm -rf "$TMP_DIR"
        exit 1
    fi

    # Step B: Normalize to relative abundance (relab)
    echo "  -> Normalizing to relative abundance (relab)..."
    humann_renorm_table \
        --input "$OUTPUT_FILE" \
        --units relab \
        --output "$RELAB_FILE"

    if [ $? -ne 0 ]; then
        echo "Error: Failed to normalize $type."
        rm -rf "$TMP_DIR"
        exit 1
    fi

    # Optional: Preview the first 5 lines of the relab file in the log
    # head -n 5 "$RELAB_FILE"

    # Step C: Stratify into function-related species and function-only
    echo "  -> Splitting into stratified and unstratified tables..."
    humann_split_stratified_table \
        --input "$RELAB_FILE" \
        --output "$OUTPUT_DIR"

    if [ $? -eq 0 ]; then
        echo "Success! Processing completed for $type"
    else
        echo "Error: Failed to split $type tables."
        rm -rf "$TMP_DIR"
        exit 1
    fi
done

# 4. Clean up temporary directory
echo "-----------------------------------------"
echo "Cleaning up temporary directory..."
rm -rf "$TMP_DIR"

echo "========================================="
echo "All operations completed successfully."
