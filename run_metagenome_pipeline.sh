#!/usr/bin/env bash

# ==================== INPUT ====================
SAMPLE_DIR=$1     # sample directory
MODE=$2           # kraken or metaphlan
MAIN_DIR=$3       # output directory
THREADS=$4        # number of cpu

CLEAN_FQ=true     # whether to clean intermediate fq files

if [[ -z "$SAMPLE_DIR" || -z "$MODE" || -z "$MAIN_DIR" || -z "$THREADS" ]]; then
    echo "Usage: $0 <sample_dir> <mode: kraken|metaphlan> <output_dir> <Threads>"
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

DB_DIR="/home/han_xiao/xiaohan/Database"
DB_KNEADDATA_HUMAN="${DB_DIR}/kneaddata/hg_39"
DB_KNEADDATA_MOUSE="${DB_DIR}/kneaddata/mouse_C57BL_6NJ"
DB_KNEADDATA_SILVA="${DB_DIR}/kneaddata/SILVA_128_LSUParc_SSUParc_ribosomal_RNA"
DB_KRAKEN2="${DB_DIR}/kraken2/k2_pluspf_20250402"
DB_METAPHLAN="${DB_DIR}/mpa422"
METAPHLAN_INDEX="mpa_vJan25_CHOCOPhlAnSGB_202503"

READ_LEN=150
TRIMMOMATIC="/home/han_xiao/xiaohan/Software/Trimmomatic-0.39"
CONDA_ACTIVATE="/home/han_xiao/xiaohan/Software/miniconda3/bin/activate"

# ==================== Step 1: Merge fq.gz ====================
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

zcat "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fq.gz" | \
    sed '1~4 s/ 1:/.1:/;1~4 s/$/\/1/' | gzip > "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fixed.fq.gz"

zcat "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fq.gz" | \
    sed '1~4 s/ 2:/.1:/;1~4 s/$/\/2/' | gzip > "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fixed.fq.gz"

rm -rf "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fq.gz" "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fq.gz"
mv "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fixed.fq.gz" "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fq.gz"
mv "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fixed.fq.gz" "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fq.gz"


# ==================== Step 2: Kneaddata host removal ====================
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
    --reorder \
    --run-fastqc-start \
    --run-fastqc-end

if [[ $CLEAN_FQ == true ]]; then
    rm -f "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_1.fq.gz" "$RAW_MERGED_DIR/$SAMPLE/${SAMPLE}_2.fq.gz"
    # rm -f "$RAW_MERGED_DIR/$sample/${sample}_1.fastq" "$RAW_MERGED_DIR/$sample/${sample}_2.fastq"
    rm -f "$KNEADDATA_DIR/$SAMPLE/*contam*" "$KNEADDATA_DIR/$SAMPLE//*unmatched*"
fi


# ==================== Step 3: Taxonomic Profiling ====================
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

    if [[ $CLEAN_FQ == true ]]; then
            rm -rf "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_1.fastq" "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_2.fastq"
			rm -rf "$KRAKEN_DIR/${SAMPLE}.output"
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
    mkdir -p "$META_DIR"

    metaphlan \
        "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_1.fastq","$KNEADDATA_DIR/$SAMPLE/${SAMPLE}_paired_2.fastq" \
        --mapout "$KNEADDATA_DIR/$SAMPLE/${SAMPLE}.bowtie2.bz2" \
        --input_type fastq \
        --db_dir "$DB_METAPHLAN" \
        --index "$METAPHLAN_INDEX" \
        --nproc "$THREADS" \
        --offline \
        -o "$META_DIR/${SAMPLE}_metaphlan_bugs_list.tsv"

    sgb_to_gtdb_profile.py -i "$META_DIR/${SAMPLE}_metaphlan_bugs_list.tsv" -o "$META_DIR/${SAMPLE}_gtdb_bugs_list.tsv"

else
    echo "ERROR: Unknown mode $MODE"
    exit 1
fi

echo ">>> Sample $SAMPLE finished successfully."

