#!/usr/bin/env bash

# ==================== STRICT ERROR HANDLING ====================
# Exit immediately if a pipeline or command exits with a non-zero status (fails).
# This prevents the script from pretending to succeed when a crucial middle step crashes.
set -e
set -o pipefail

# ==================== INPUT ARGUMENTS ====================
SAMPLE_DIR=$1     # Path to the input sample directory containing raw fastq.gz files
MODE=$2           # Profiling mode: 'kraken' (Kraken2+Bracken) or 'metaphlan' (MetaPhlAn4)
MAIN_DIR=$3       # Main output directory for the pipeline results
THREADS=$4        # Number of CPU threads to use for parallel processing

CLEAN_INTERMEDIATE=true     # Flag to delete large intermediate fastq files to save disk space
CLEAN_FQ=false

# Check if all 4 required arguments are provided; exit with usage instructions if not
if [[ -z "$SAMPLE_DIR" || -z "$MODE" || -z "$MAIN_DIR" || -z "$THREADS" ]]; then
    echo "Usage: $0 <sample_dir> <mode: kraken|metaphlan> <output_dir> <Threads>"
    exit 1
fi

# Extract the base name of the sample directory to use as the unique sample ID
SAMPLE=$(basename "$SAMPLE_DIR")
echo ">>> Running sample: $SAMPLE (mode=$MODE)"

# ==================== PATHS & SETTINGS ====================
# Define output subdirectories for each step of the pipeline
RAW_MERGED_DIR="${MAIN_DIR}/01.fq"
KNEADDATA_DIR="${MAIN_DIR}/02.kneaddata"
KRAKEN_DIR="${MAIN_DIR}/03.kraken2"
BRACKEN_DIR="${MAIN_DIR}/04.bracken"
META_DIR="${MAIN_DIR}/05.metaphlan"

# Define paths to reference databases
DB_DIR="/home/han_xiao/xiaohan/Database"
DB_KNEADDATA_HUMAN="${DB_DIR}/kneaddata/hg_39"                                    # Human genome for decontamination
DB_KNEADDATA_MOUSE="${DB_DIR}/kneaddata/mouse_C57BL_6NJ"                          # Mouse genome for decontamination
DB_KNEADDATA_SILVA="${DB_DIR}/kneaddata/SILVA_128_LSUParc_SSUParc_ribosomal_RNA"  # rRNA database for decontamination
DB_KRAKEN2="${DB_DIR}/kraken2/k2_pluspf_20250402"                                 # Kraken2 standard/plusPF database
DB_METAPHLAN="${DB_DIR}/mpa422"                                                   # MetaPhlAn4 database directory
METAPHLAN_INDEX="mpa_vJan25_CHOCOPhlAnSGB_202503"                                 # MetaPhlAn4 specific index name

READ_LEN=150 # Expected read length, used for Bracken abundance estimation
TRIMMOMATIC="/home/han_xiao/xiaohan/Software/Trimmomatic-0.39"
CONDA_ACTIVATE="/home/han_xiao/xiaohan/Software/miniconda3/bin/activate"

# ==================== Step 1: Merge raw fastq files ====================
# This step merges multiple sequencing lanes into single R1 and R2 files per sample
mkdir -p "$RAW_MERGED_DIR/$SAMPLE"

# Find all R1 and R2 fastq.gz files in the sample directory and sort them alphabetically
fq1=($(ls "$SAMPLE_DIR"/*_1.fq.gz 2>/dev/null | sort))
fq2=($(ls "$SAMPLE_DIR"/*_2.fq.gz 2>/dev/null | sort))

# Ensure sequence files were actually found
if [[ ${#fq1[@]} -eq 0 || ${#fq2[@]} -eq 0 ]]; then
    echo "ERROR: No fq.gz found for $SAMPLE"
    exit 1
fi

echo ">>> Merging lanes for $SAMPLE"
# Concatenate all R1 files together, and all R2 files together
cat "${fq1[@]}" > "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fq.gz"
cat "${fq2[@]}" > "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fq.gz"

# Reformat sequence IDs to ensure compatibility with downstream paired-end tools.
# sed '1~4 ...' targets the header line of each fastq record (every 4th line starting from line 1)
# It replaces ' 1:' with '.1:' and appends '/1' or '/2' to denote read pair direction.
zcat "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fq.gz" | \
    sed '1~4 s/ 1:/.1:/;1~4 s/$/\/1/' | gzip > "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fixed.fq.gz"

zcat "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fq.gz" | \
    sed '1~4 s/ 2:/.1:/;1~4 s/$/\/2/' | gzip > "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fixed.fq.gz"

# Remove the original merged files and rename the fixed files to the standard names
rm -rf "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fq.gz" "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fq.gz"
mv "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fixed.fq.gz" "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fq.gz"
mv "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fixed.fq.gz" "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fq.gz"


# ==================== Step 2: Quality Control & Host Removal ====================
mkdir -p "$KNEADDATA_DIR/$SAMPLE"

# Activate the conda environment containing kneaddata
source "$CONDA_ACTIVATE" kneaddate

echo ">>> Running Kneaddata for $SAMPLE"
# Run kneaddata to perform quality trimming (via Trimmomatic) and 
# remove host (Human/Mouse) and ribosomal RNA contamination (via Bowtie2)
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
    --reorder \
    --run-fastqc-start \
    --run-fastqc-end

# Clean up intermediate fastq files to save disk space
if [[ $CLEAN_INTERMEDIATE == true ]]; then
    # Delete the raw merged fastq files as we now have the clean kneaddata output
    rm -f "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fq.gz" "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fq.gz"
    
    # Delete the contaminated/unmatched sequence outputs from Bowtie2 alignments.
    # Note: Asterisks (*) are kept OUTSIDE double quotes so Bash correctly evaluates them as wildcards.
    rm -f "$KNEADDATA_DIR/$SAMPLE/"*contam* "$KNEADDATA_DIR/$SAMPLE/"*unmatched*
fi


# ==================== Step 3: Taxonomic Profiling ====================
if [[ "$MODE" == "kraken" ]]; then
    echo ">>> Running Kraken2 + Bracken for $SAMPLE"
    # Activate the conda environment containing Kraken2 and Bracken
    source "$CONDA_ACTIVATE" kraken2.1.5
    mkdir -p "$KRAKEN_DIR" "$BRACKEN_DIR"

    # Run Kraken2 for k-mer based taxonomic classification
    kraken2 \
        --db "$DB_KRAKEN2" \
        --paired "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_1.fastq" "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_2.fastq" \
        --threads "$THREADS" \
        --use-names --report-zero-counts \
        --report "$KRAKEN_DIR/${SAMPLE}.report" \
        --output "$KRAKEN_DIR/${SAMPLE}.output"

    # Convert Kraken2 report to MetaPhlAn-style format for easier downstream plotting
    kreport2mpa.py -r "$KRAKEN_DIR/${SAMPLE}.report" --display-header -o "$KRAKEN_DIR/${SAMPLE}.mpa"

    rm -rf "$KRAKEN_DIR/${SAMPLE}.output" # The read-by-read classification output (.output) is huge and rarely needed downstream

    # Clean up massive intermediate files generated by Kraken and Kneaddata
    if [[ $CLEAN_FQ == true ]]; then
        rm -rf "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_1.fastq" "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_2.fastq"
    fi

    # Run Bracken to re-estimate relative abundances at Phylum (P), Genus (G), and Species (S) levels
    for tax in P G S; do
        mkdir -p "$BRACKEN_DIR/$tax"
        bracken \
            -d "$DB_KRAKEN2" \
            -i "$KRAKEN_DIR/${SAMPLE}.report" \
            -r "$READ_LEN" -l "$tax" -t 0 \
            -o "$BRACKEN_DIR/$tax/${SAMPLE}.brk"

        # Extract the fraction of total reads (abundance) and format it into a clean count table
        tail -n+2 "$BRACKEN_DIR/$tax/${SAMPLE}.brk" | LC_ALL=C sort | cut -f6 | sed "1 s/^/$SAMPLE\n/" > "$BRACKEN_DIR/$tax/${SAMPLE}.count"
    done

elif [[ "$MODE" == "metaphlan" ]]; then
    echo ">>> Running MetaPhlAn4 for $SAMPLE"
    # Activate the conda environment containing MetaPhlAn 4.2.2
    source "$CONDA_ACTIVATE" mpa4.2.2
    mkdir -p "$META_DIR"

    # Run MetaPhlAn4 for marker-gene based taxonomic profiling
    metaphlan \
        "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_1.fastq","$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_2.fastq" \
        --mapout "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}.bowtie2.bz2" \
        --input_type fastq \
        --db_dir "$DB_METAPHLAN" \
        --index "$METAPHLAN_INDEX" \
        --nproc "$THREADS" \
        --offline \
        -o "$META_DIR/${SAMPLE}_metaphlan_bugs_list.tsv"

    # Convert MetaPhlAn SGB (Species-level Genome Bins) profile to standard GTDB taxonomy
    sgb_to_gtdb_profile.py -i "$META_DIR/${SAMPLE}_metaphlan_bugs_list.tsv" -o "$META_DIR/${SAMPLE}_gtdb_bugs_list.tsv"

    # Clean up large intermediate fastq files from Kneaddata to save disk space
    if [[ $CLEAN_FQ == true ]]; then
        rm -rf "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_1.fastq" "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_2.fastq"
    fi

else
    # Catch-all for unsupported modes (e.g. if user inputs 'humann3' instead of kraken/metaphlan)
    echo "ERROR: Unknown mode $MODE"
    exit 1
fi

# Log successful completion. 
# Because 'set -e' is active at the top, this line will ONLY print if every preceding step was 100% successful.
echo ">>> Sample $SAMPLE finished successfully."
