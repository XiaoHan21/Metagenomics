#!/usr/bin/env bash
# ==============================================================
# merge_bracken_counts.sh
# Merge individual sample *.count files into combined OTU tables
# and remove human/Chordata contamination.
#
# Usage:
#   bash merge_bracken_counts.sh <bracken_dir> <prop>
#
# Example:
#   bash merge_bracken_counts.sh /path/to/04.bracken 0.001
#
# Dependencies:
#   - Rscript environment containing filter_feature_table.R
# ==============================================================

BRACKEN_DIR=$1      # e.g. /home/han_xiao/xiaohan/Work/02-SpatialGF/03.result/21.Stool_metagenomics/04.bracken 
PROP=$2             # e.g. 0.1 (prevalance filter threshold for R script)
RSCRIPT="/home/han_xiao/xiaohan/Software/miniconda3/envs/R4.4/bin/Rscript"
FILTER_SCRIPT="/home/han_xiao/xiaohan/Software/EasyMicrobiome-1.20/script/filter_feature_table.R"

# Check input arguments
if [[ -z "$BRACKEN_DIR" || -z "$PROP" ]]; then
    echo "Usage: $0 <bracken_dir> <prop>"
    exit 1
fi


# Loop through taxonomy levels: Phylum (P), Genus (G), Species (S)
for tax in P G S; do
    TAX_DIR="${BRACKEN_DIR}/${tax}"
    echo ">>> Processing level: ${tax}"

    # 1. Collect sample names (sorted)
    samples=($(ls -v "$TAX_DIR"/*.count | xargs -n1 basename | sed 's/.count//'))

    if [[ ${#samples[@]} -eq 0 ]]; then
        echo "No .count files found in ${TAX_DIR}, skipping..."
        continue
    fi

    # 2. Create header (replace first line with 'Taxonomy')
    tail -n+2 "$TAX_DIR/${samples[0]}.brk" | LC_ALL=C sort | cut -f1 | sed "1 s/^/Taxonomy\n/" > "$TAX_DIR/header.txt"

    # 3. Merge all sample count files into one OTU table
    paste "$TAX_DIR/header.txt" "$TAX_DIR"/*.count > "$BRACKEN_DIR/bracken.${tax}.txt"
	
    # 4. Filter low-abundance taxa using R script
    "$RSCRIPT" "$FILTER_SCRIPT" \
        -i "$BRACKEN_DIR/bracken.${tax}.txt" \
        -p "$PROP" \
        -o "$BRACKEN_DIR/bracken.${tax}.${PROP}.txt"
done

# 5. Remove human-related taxa
grep -v 'Chordata'     "$BRACKEN_DIR/bracken.P.${PROP}.txt" > "$BRACKEN_DIR/bracken.P.${PROP}-H.txt"
grep -v "Homo"         "$BRACKEN_DIR/bracken.G.${PROP}.txt" > "$BRACKEN_DIR/bracken.G.${PROP}-H.txt"
grep -v 'Homo sapiens' "$BRACKEN_DIR/bracken.S.${PROP}.txt" > "$BRACKEN_DIR/bracken.S.${PROP}-H.txt"

echo ">>> All taxonomy levels processed successfully."

