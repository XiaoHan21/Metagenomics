#!/bin/bash
# ===============================================
#  Submit all metagenome samples to HPC via qsub
# ===============================================

core_pipeline="/home/han_xiao/xiaohan/Work/02-SpatialGF/02.script/Stool_Metagenomics/run_metagenome_pipeline.sh"

RAW_DIR="/home/han_xiao/xiaohan/Work/02-SpatialGF/01.data/Stool_metagenomics/01.RawData"
OUT_DIR="/home/han_xiao/xiaohan/Work/02-SpatialGF/03.result/21.Stool_metagenomics/"
PIPELINE=${1:-kraken}   # options：kraken 或 metaphlan
THREADS=${2:-8}

cd $OUT_DIR

for f in "$RAW_DIR"/*; do
  SAMPLE=$(basename "$f")
  JOB="job_${SAMPLE}_${PIPELINE}.sh"

<<'COM'

  cat > "$JOB" <<EOF
#!/bin/bash
#PBS -N ${SAMPLE}
#PBS -P as_lkc_sunny.wong
#PBS -l select=1:ncpus=$THREADS:mem=50g
#PBS -l walltime=240:00:00
#PBS -e ${SAMPLE}.e
#PBS -o ${SAMPLE}.o

cd ${OUT_DIR}
bash ${core_pipeline} "$RAW_DIR/$SAMPLE" "$PIPELINE" "$OUT_DIR" "$THREADS"
EOF

  echo "Submitting $JOB..."
  # qsub "$JOB"
COM

cd ${OUT_DIR}
bash ${core_pipeline} "$RAW_DIR/$SAMPLE" "$PIPELINE" "$OUT_DIR" "$THREADS"

done

