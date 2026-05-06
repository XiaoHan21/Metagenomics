#!/bin/bash
# ===============================================
#  Submit all metagenome samples (Smart Array/Single job)
#  - Automatically switches between array and single job
#  - Prevents log file overwriting with timestamps
# ===============================================

# Define the core analysis script and directories
core_pipeline="/home/han_xiao/xiaohan/Software/Pipeline/Metagenomics/02.script/run_metagenome_pipeline.sh"
RAW_DIR="/home/han_xiao/xiaohan/Software/Pipeline/Metagenomics/01.data/02.demo/fq"
OUT_DIR="/home/han_xiao/xiaohan/Software/Pipeline/Metagenomics/03.result/02.demo/"

# Parse command-line arguments (defaults: pipeline=humann, threads=16)
PIPELINE=${1:-humann}
THREADS=${2:-16}

# Set the maximum number of concurrent jobs for the PBS array
MAX_JOB_NUM=5

# Ensure the main output directory exists
mkdir -p "$OUT_DIR"

# Collect all sample names from the raw data directory into an array
mapfile -t SAMPLES < <(ls "$RAW_DIR")
NUM=${#SAMPLES[@]}

# Exit immediately if no samples are found
if [ "$NUM" -eq 0 ]; then
    echo "Error: No samples found in $RAW_DIR"
    exit 1
fi

# ================= CORE FIX: Smart Array Detection =================
# If only 1 sample is found, the -J parameter is omitted to prevent a PBS syntax error.
if [ "$NUM" -eq 1 ]; then
    PBS_ARRAY_LINE=""
    echo "Only 1 sample detected. Submitting as a standard single job..."
else
    PBS_ARRAY_LINE="#PBS -J 1-${NUM}%${MAX_JOB_NUM}"
    echo "Submitting job array with ${NUM} samples (max ${MAX_JOB_NUM} concurrent)..."
fi
# ===================================================================

# Define the name of the generated submission script
ARRAY_JOB="run_metagenome_array.sh"

# Generate the PBS job script
# Note: Inside <<EOF, variables with an escaped \$ will be preserved 
# in the child script and evaluated only when running on the compute node.
cat > "$ARRAY_JOB" <<EOF
#!/bin/bash
#PBS -N meta_array_${PIPELINE}
#PBS -P as_lkc_sunny.wong
#PBS -l select=1:ncpus=${THREADS}:mem=100g
#PBS -l walltime=240:00:00
#PBS -q qintel_wfly
${PBS_ARRAY_LINE}

# Load the sample array (hardcodes all sample names into the child script at submission)
SAMPLES=(${SAMPLES[@]})

# ================= Compatibility Handling =================
# For a single job submission, PBS_ARRAY_INDEX will be empty.
# We use the \${VAR:-1} syntax to default the index to 1 if it is empty.
CURRENT_INDEX=\${PBS_ARRAY_INDEX:-1}

# Extract the sample name corresponding to the current job (Bash arrays are 0-indexed)
SAMPLE=\${SAMPLES[\$((CURRENT_INDEX - 1))]}
# ==========================================================

# ================= Log Isolation and Redirection =================
# 1. Create a dedicated logs directory to keep the main output folder clean
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "\$LOG_DIR"

# 2. Get the actual execution timestamp on the compute node
TIMESTAMP=\$(date +"%Y%m%d_%H%M%S")

# 3. Redirect standard output (1) and standard error (2) in real-time to: sample_timestamp.o/e
exec 1> "\$LOG_DIR/\${SAMPLE}_\${TIMESTAMP}.o"
exec 2> "\$LOG_DIR/\${SAMPLE}_\${TIMESTAMP}.e"
# =================================================================

cd "${OUT_DIR}"

# Write job execution details to the log file
echo "=========================================================="
echo "Starting analysis for sample: \$SAMPLE"
echo "Start time: \$(date)"
echo "PBS Job ID: \$PBS_JOBID"
echo "Compute Node: \$(hostname)"
echo "Pipeline Mode: ${PIPELINE}"
echo "=========================================================="

# Run the core single-sample analysis script
bash "${core_pipeline}" "${RAW_DIR}/\${SAMPLE}" "${PIPELINE}" "${OUT_DIR}" "${THREADS}"

echo "=========================================================="
echo "Finished sample: \$SAMPLE"
echo "End time: \$(date)"
echo "=========================================================="
EOF

# Provide feedback and submit the script to the queue
echo "Successfully generated: $ARRAY_JOB"
qsub "$ARRAY_JOB"
