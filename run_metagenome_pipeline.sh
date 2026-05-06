#!/usr/bin/env bash

# ==================== STRICT ERROR HANDLING ====================
# Exit immediately on error to prevent generating incomplete or corrupted results
set -e
set -o pipefail

# ==================== INPUT ARGUMENTS ====================
SAMPLE_DIR=$1     # Path to the input sample directory (containing raw fastq.gz)
MODE=$2           # Profiling mode: 'kraken', 'metaphlan', or 'humann'
MAIN_DIR=$3       # Main output directory for pipeline results
THREADS=$4        # Number of CPU threads to use

CLEAN_INTERMEDIATE=true     # Flag to clean up large intermediate files from host removal/QC
CLEAN_FQ=false              # Flag to clean up final Kneaddata clean fastq files (depends on available disk space)

if [[ -z "$SAMPLE_DIR" || -z "$MODE" || -z "$MAIN_DIR" || -z "$THREADS" ]]; then
    echo "Usage: $0 <sample_dir> <mode: kraken|metaphlan|humann> <output_dir> <Threads>"
    exit 1
fi

SAMPLE=$(basename "$SAMPLE_DIR")
echo ">>> Running sample: $SAMPLE (mode=$MODE)"

# ==================== PATHS & SETTINGS ====================
RAW_MERGED_DIR="${MAIN_DIR}/01.fq"
KNEADDATA_DIR="${MAIN_DIR}/02.kneaddata"
KRAKEN_DIR="${MAIN_DIR}/03.kraken2"
BRACKEN_DIR="${MAIN_DIR}/04.bracken"
META_DIR="${MAIN_DIR}/05.metaphlan"
HUMANN_DIR="${MAIN_DIR}/06.humann"

# Database paths (Ensure these paths exist on your server)
DB_DIR="/home/han_xiao/xiaohan/Database"
DB_KNEADDATA_HUMAN="${DB_DIR}/kneaddata/hg_39"
DB_KNEADDATA_MOUSE="${DB_DIR}/kneaddata/mouse_C57BL_6NJ"
DB_KNEADDATA_SILVA="${DB_DIR}/kneaddata/SILVA_128_LSUParc_SSUParc_ribosomal_RNA"
DB_KRAKEN2="${DB_DIR}/kraken2/k2_pluspf_20250402"
DB_METAPHLAN="${DB_DIR}/mpa422"
METAPHLAN_INDEX="mpa_vJan25_CHOCOPhlAnSGB_202503"

HUMAN_DB_METAPHLAN="${DB_DIR}/mpa411"
HUMAN_METAPHLAN_INDEX="mpa_vOct22_CHOCOPhlAnSGB_202403"

READ_LEN=150
TRIMMOMATIC="/home/han_xiao/xiaohan/Software/Trimmomatic-0.39"
CONDA_ACTIVATE="/home/han_xiao/xiaohan/Software/miniconda3/bin/activate"

# ==================== Step 1: Merge raw fastq files ====================
mkdir -p "$RAW_MERGED_DIR/$SAMPLE"

fq1=($(ls "$SAMPLE_DIR"/*_1.fq.gz 2>/dev/null | sort))
fq2=($(ls "$SAMPLE_DIR"/*_2.fq.gz 2>/dev/null | sort))

if [[ ${#fq1[@]} -eq 0 || ${#fq2[@]} -eq 0 ]]; then
    echo "ERROR: No fq.gz found for $SAMPLE"
    exit 1
fi

echo ">>> Merging lanes for $SAMPLE"
cat "${fq1[@]}" > "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fq.gz"
cat "${fq2[@]}" > "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fq.gz"

# Fix paired-end sequence IDs to ensure compatibility with downstream tools
zcat "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fq.gz" | \
    sed '1~4 s/ 1:/.1:/;1~4 s/$/\/1/' | gzip > "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fixed.fq.gz"
zcat "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fq.gz" | \
    sed '1~4 s/ 2:/.1:/;1~4 s/$/\/2/' | gzip > "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fixed.fq.gz"

rm -rf "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fq.gz" "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fq.gz"
mv "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fixed.fq.gz" "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fq.gz"
mv "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fixed.fq.gz" "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fq.gz"


# ==================== Step 2: Quality Control & Host Removal ====================
mkdir -p "$KNEADDATA_DIR/$SAMPLE"
source "$CONDA_ACTIVATE" kneaddate

echo ">>> Running Kneaddata for $SAMPLE"
kneaddata \
    -i1 "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fq.gz" \
    -i2 "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fq.gz" \
    -o "$KNEADDATA_DIR/$SAMPLE" \
    --reference-db "$DB_KNEADDATA_HUMAN" \
    --reference-db "$DB_KNEADDATA_MOUSE" \
    --reference-db "$DB_KNEADDATA_SILVA" \
    --output-prefix "$SAMPLE" \
    -t "$THREADS" \
    --trimmomatic "$TRIMMOMATIC" \
    --bowtie2-options "--very-sensitive --dovetail -p $THREADS" \
    --remove-intermediate-output \
    --max-memory 70g \
    --reorder

if [[ $CLEAN_INTERMEDIATE == true ]]; then
    rm -f "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fq.gz" "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fq.gz"
    rm -f "$KNEADDATA_DIR/$SAMPLE/"*contam* "$KNEADDATA_DIR/$SAMPLE/"*unmatched*
fi

# ==================== Step 3: Taxonomic & Functional Profiling ====================

if [[ "$MODE" == "kraken" ]]; then
    echo ">>> Running Kraken2 + Bracken for $SAMPLE"
    source "$CONDA_ACTIVATE" kraken2.1.5
    mkdir -p "$KRAKEN_DIR" "$BRACKEN_DIR"

    kraken2 \
        --db "$DB_KRAKEN2" \
        --paired "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_1.fastq" "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_2.fastq" \
        --threads "$THREADS" \
        --use-names --report-zero-counts \
        --report "$KRAKEN_DIR/${SAMPLE}.report" \
        --output "$KRAKEN_DIR/${SAMPLE}.output"

    kreport2mpa.py -r "$KRAKEN_DIR/${SAMPLE}.report" --display-header -o "$KRAKEN_DIR/${SAMPLE}.mpa"
    rm -rf "$KRAKEN_DIR/${SAMPLE}.output"

    if [[ $CLEAN_FQ == true ]]; then
        rm -rf "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_1.fastq" "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_2.fastq"
    fi

    for tax in P G S; do
        mkdir -p "$BRACKEN_DIR/$tax"
        bracken \
            -d "$DB_KRAKEN2" \
            -i "$KRAKEN_DIR/${SAMPLE}.report" \
            -r "$READ_LEN" -l "$tax" -t 0 \
            -o "$BRACKEN_DIR/$tax/${SAMPLE}.brk"

        tail -n+2 "$BRACKEN_DIR/$tax/${SAMPLE}.brk" | LC_ALL=C sort | cut -f6 | sed "1 s/^/$SAMPLE\n/" > "$BRACKEN_DIR/$tax/${SAMPLE}.count"
    done

elif [[ "$MODE" == "metaphlan" ]]; then
    echo ">>> Running MetaPhlAn4 for $SAMPLE"
    source "$CONDA_ACTIVATE" mpa4.2.2
    
    # Create sample-specific subdirectory for MetaPhlAn
    mkdir -p "$META_DIR/$SAMPLE"

    metaphlan \
        "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_1.fastq","$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_2.fastq" \
        --mapout "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}.bowtie2.bz2" \
        --input_type fastq \
        --db_dir "$DB_METAPHLAN" \
        --index "$METAPHLAN_INDEX" \
        --nproc "$THREADS" \
        --offline \
        -o "$META_DIR/$SAMPLE/${SAMPLE}_metaphlan_bugs_list.tsv"

    sgb_to_gtdb_profile.py -i "$META_DIR/$SAMPLE/${SAMPLE}_metaphlan_bugs_list.tsv" -o "$META_DIR/$SAMPLE/${SAMPLE}_gtdb_bugs_list.tsv"

    if [[ $CLEAN_FQ == true ]]; then
        rm -rf "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_1.fastq" "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_2.fastq"
    fi

elif [[ "$MODE" == "humann" ]]; then
    echo ">>> Running MetaPhlAn4 + HUMAnN4 for $SAMPLE"
    source "$CONDA_ACTIVATE" humann4_mpa4.1.1
    
    # Create sample-specific subdirectories for both HUMAnN and MetaPhlAn
    mkdir -p "$HUMANN_DIR/$SAMPLE" "$META_DIR/$SAMPLE"

    # Concatenate paired-end reads as input for HUMAnN4
    CONCAT_FQ="$KNEADDATA_DIR/$SAMPLE/${SAMPLE}.fastq"
    echo ">>> Concatenating paired reads..."
    cat "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_1.fastq" \
        "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_2.fastq" > "$CONCAT_FQ"

    # Run HUMAnN4 (calls MetaPhlAn4 under the hood)
    # Output is directed to the sample-specific directory
    humann \
        --input "$CONCAT_FQ" \
        --threads "$THREADS" \
        --metaphlan-options "--input_type fastq --bowtie2db $HUMAN_DB_METAPHLAN --index $HUMAN_METAPHLAN_INDEX --offline -t rel_ab_w_read_stats --nproc $THREADS" \
        --output "$HUMANN_DIR/$SAMPLE"

    # ---------------------------------------------------------
    # Move and standardize the MetaPhlAn profile output
    # HUMAnN generated profiles usually look like: P0_D1_1_metaphlan_profile.tsv
    # We use a wildcard to capture it dynamically.
    # ---------------------------------------------------------
    MPA_PROFILE=$(ls "$HUMANN_DIR/$SAMPLE"/${SAMPLE}_*metaphlan_profile.tsv 2>/dev/null | head -n 1)
    
    if [[ -n "$MPA_PROFILE" && -f "$MPA_PROFILE" ]]; then
        echo ">>> Moving and converting MetaPhlAn profile to $META_DIR/$SAMPLE..."
        
        # Move the profile to the 05.metaphlan sample folder and rename it for consistency
        mv "$MPA_PROFILE" "$META_DIR/$SAMPLE/${SAMPLE}_metaphlan_bugs_list.tsv"
        
        # Convert the SGB profile to GTDB taxonomy format
        sgb_to_gtdb_profile.py -i "$META_DIR/$SAMPLE/${SAMPLE}_metaphlan_bugs_list.tsv" \
                               -o "$META_DIR/$SAMPLE/${SAMPLE}_gtdb_bugs_list.tsv"
    else
        echo "WARNING: MetaPhlAn profile (*_metaphlan_profile.tsv) not found in HUMAnN dir for $SAMPLE"
    fi

    # ---------------------------------------------------------
    # Clean up intermediate and temporary files
    # ---------------------------------------------------------
    echo ">>> Cleaning up HUMAnN intermediate temp directory..."
    # Remove the large temp folder containing intermediate sam, fasta, and bt2 files
    rm -rf "$HUMANN_DIR/$SAMPLE/${SAMPLE}_humann_temp"

    if [[ $CLEAN_FQ == true ]]; then
        echo ">>> Cleaning up KneadData fastq files..."
        rm -f "$CONCAT_FQ"
        rm -rf "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_1.fastq" "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_2.fastq"
    fi

else
    echo "ERROR: Unknown mode $MODE"
    exit 1
fi