#!/bin/bash
#SBATCH -J pbrun_metrics_a100
#SBATCH --gres=gpu:a100:1        # request 1 A100
#SBATCH --cpus-per-task=32       # match V100 node cores
#SBATCH --mem=180G               # total RAM for the job
#SBATCH -t 01:59:00              # walltime
#SBATCH -o pbrun_met_%j.out      # stdout
#SBATCH -e pbrun_met_%j.err      # stderr
#SBATCH --mail-type=END,FAIL     # email at job end/fail
#SBATCH --mail-user=qaedi65@gmail.com
# -- optionally pin to a particular node:
# #SBATCH --nodelist=cac107

set -euo pipefail
export LANG=C

# ---------- User-editable variables ----------
# Container and runtime
SIF_IMAGE="/global/project/hpcg6049/somatic_pipeline/containers/NXF_SINGULARITY_CACHEDIR/parabrick_4.5.1.sif"
# Where you normally run (bind this)
WORK_DIR="$(pwd)"

# Sample & files
SAMPLE_ID="KQ0202-DNA"

# Input BAM (Output from the previous fq2bam run)
INPUT_BAM="/global/scratch/hpc6049/251031_A01939_0053_BHGK7JDSXF/fq2bam_gpu_outs/KQ0202-DNA.bam"
REF_FILE="/global/project/hpcg6049/somatic_pipeline/data/references/GRCh38.d1.vd1.fa"

# Output Directory
OUT_DIR="/global/scratch/hpc6049/251031_A01939_0053_BHGK7JDSXF/fq2bam_gpu_outs"

# Parabricks options 
# Left empty because the docs confirm no specific extra flags are needed for standard WGS metrics
PBRUN_EXTRA_OPTS="" 

# ---------- End user-editable variables ----------

echo "Job starting on $(hostname) at $(date)"
echo "SAMPLE_ID=${SAMPLE_ID}"
echo "SIF_IMAGE=${SIF_IMAGE}"

module purge
module load apptainer

# create output dir if missing
mkdir -p "${OUT_DIR}"
chmod u+w "${OUT_DIR}" || true

# Build bind list (dedupe parent directories)
declare -A _dirs
for f in "$WORK_DIR" "$REF_FILE" "$INPUT_BAM" "$OUT_DIR"; do
  [[ -z "${f:-}" ]] && continue
  d=$(dirname "$f")
  _dirs["$d"]=1
done

bindlist=""
for d in "${!_dirs[@]}"; do
  # ensure host:container same path mapping
  bindlist+="${d}:${d},"
done
bindlist="${bindlist%,}"    # trim trailing comma

echo "Will bind these directories:"
for d in "${!_dirs[@]}"; do echo "  $d"; done
echo "Bind string: $bindlist"
echo

# Pre-flight checks INSIDE container
check_cmd=""
for f in "$REF_FILE" "$INPUT_BAM" "$OUT_DIR"; do
  esc_f="${f//\'/\'\\\'\'}"
  check_cmd+="ls -ld '${esc_f}' || echo 'MISSING: ${esc_f}'; "
done
check_cmd+="echo '--- CONTAINER nvidia-smi ---'; nvidia-smi || true; echo '--- CONTAINER ENV ---'; env | grep -E 'CUDA|NVIDIA|LD_LIBRARY_PATH' || true; id; groups;"

echo "Running pre-flight checks inside container..."
apptainer exec --nv --bind "${bindlist}" "${SIF_IMAGE}" bash -lc "${check_cmd}"

echo "Capturing host and container GPU state..."
nvidia-smi -q > nvidia_host_q_${SLURM_JOB_ID:-local}.txt || true
apptainer exec --nv "${SIF_IMAGE}" nvidia-smi -q > nvidia_container_q_${SLURM_JOB_ID:-local}.txt || true

# Launch Parabricks, capture full logs
PBRUN_LOG="${OUT_DIR}/${SAMPLE_ID}_bammetrics.${SLURM_JOB_ID:-local}.log"
echo "Starting Parabricks bammetrics run; logs -> ${PBRUN_LOG}"
set -x

# NOTE: set -o pipefail ensures that if pbrun fails, the script sees the failure code, not the tee success code
apptainer exec --nv --bind "${bindlist}" "${SIF_IMAGE}" \
  bash -lc "set -o pipefail; \
    echo 'Parabricks bammetrics starting on $(hostname) at $(date)'; \
    nvidia-smi; \
    pbrun bammetrics \
      --ref '${REF_FILE}' \
      --bam '${INPUT_BAM}' \
      --out-metrics-file '${OUT_DIR}/${SAMPLE_ID}_wgs_metrics.txt' \
      ${PBRUN_EXTRA_OPTS} \
      2>&1 | tee '${PBRUN_LOG}'"

rc=${PIPESTATUS[0]:-0}
set +x

if [[ $rc -ne 0 ]]; then
  echo "Parabricks bammetrics failed with rc=$rc. See ${PBRUN_LOG} and Parabricks stderr for details."
  exit $rc
else
  echo "Parabricks bammetrics finished successfully at $(date)."
  echo "Metrics written to: ${OUT_DIR}/${SAMPLE_ID}_wgs_metrics.txt"
  exit 0
fi