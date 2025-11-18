#!/bin/bash
# ===============================================
#  Submit all metagenome samples via PBS job array
#  - Automatically generate a job array
#  - Max simultaneous jobs = 5
# ===============================================

core_pipeline="/home/kaiyee001/scripts/run_metagenome_pipeline.sh"

RAW_DIR="/home/kaiyee001/AIPS_Raw_Read/Pending_run"
OUT_DIR="/home/kaiyee001/Metagenomics_Results"
PIPELINE=${1:-metaphlan}   # options: kraken 或 metaphlan
THREADS=${2:-16}

mkdir -p $OUT_DIR

# 收集所有样本名
mapfile -t SAMPLES < <(ls "$RAW_DIR")
NUM=${#SAMPLES[@]}

# 生成 job array 脚本
ARRAY_JOB="run_metagenome_array.sh"

cat > "$ARRAY_JOB" <<EOF
#!/bin/bash
#PBS -N meta_array_${PIPELINE}
#PBS -P as_lkc_sunny.wong
#PBS -l select=1:ncpus=${THREADS}:mem=100g
#PBS -l walltime=240:00:00
#PBS -q qintel_wfly
#PBS -J 1-${NUM}%5
#PBS -e ${OUT_DIR}/array_job.e
#PBS -o ${OUT_DIR}/array_job.o

# 加载样本数组
SAMPLES=(${SAMPLES[@]})

# 当前 job 对应的样本（PBS_ARRAY_INDEX 从 1 开始）
SAMPLE=\${SAMPLES[\$((PBS_ARRAY_INDEX - 1))]}

cd ${OUT_DIR}

echo "Running sample: \$SAMPLE"
bash ${core_pipeline} "${RAW_DIR}/\${SAMPLE}" "${PIPELINE}" "${OUT_DIR}" "${THREADS}"
EOF

# 提交 array job
echo "Submitting job array with ${NUM} samples (max 5 concurrent)..."
qsub "$ARRAY_JOB"
