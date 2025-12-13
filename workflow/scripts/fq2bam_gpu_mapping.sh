#!/bin/bash
#SBATCH -J pbrun_ma32_a100
#SBATCH --gres=gpu:a100:1       # request 1 A100
#SBATCH --cpus-per-task=32       # match V100 node cores (adjust if needed)
#SBATCH --mem=180G               # total RAM for the job
#SBATCH -t 05:59:00              # walltime
#SBATCH -o pbrun_%j.out          # stdout
#SBATCH -e pbrun_%j.err          # stderr
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

# Sample & files (edit these)
SAMPLE_ID="KQ0202-DNA"
FASTQ1="/global/scratch/hpc6049/251031_A01939_0053_BHGK7JDSXF/fastp_outs/trimmed_fqs/KQ0202-DNA.trimmed_1.fastq.gz"
FASTQ2="/global/scratch/hpc6049/251031_A01939_0053_BHGK7JDSXF/fastp_outs/trimmed_fqs/KQ0202-DNA.trimmed_2.fastq.gz"
REF_FILE="/global/project/hpcg6049/somatic_pipeline/data/references/GRCh38.d1.vd1.fa"
KNOWN_SITES_MILLS="/global/project/hpcg6049/somatic_pipeline/data/references/Mills.d1vd1.ready.vcf.gz"
KNOWN_SITES_1000G="/global/project/hpcg6049/somatic_pipeline/data/references/1000G.indels.d1vd1.ready.vcf.gz"
KNOWN_SITES_DBSNP="//global/project/hpcg6049/somatic_pipeline/data/references/dbSNP.GCF40.d1vd1.ready.vcf.gz"
OUT_DIR="/global/scratch/hpc6049/251031_A01939_0053_BHGK7JDSXF/fq2bam_gpu_outs"



# Parabricks options (tweak as needed)
BWA_OPTIONS="-Y"
PBRUN_EXTRA_OPTS=""   # e.g. "--num-streams 2" if your version supports lowering streams

# Export these if you want to force single GPU (uncomment to use)
# export CUDA_VISIBLE_DEVICES=0

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
for f in "$WORK_DIR" "$REF_FILE" "$KNOWN_SITES_MILLS" "$KNOWN_SITES_1000G" "$KNOWN_SITES_DBSNP" "$FASTQ1" "$FASTQ2" "$OUT_DIR"; do
  [[ -z "${f:-}" ]] && continue
  d=$(dirname "$f")
  _dirs["$d"]=1
done

bindlist=""
for d in "${!_dirs[@]}"; do
  # ensure host:container same path mapping (keeps paths identical inside container)
  bindlist+="${d}:${d},"
done
bindlist="${bindlist%,}"   # trim trailing comma

echo "Will bind these directories:"
for d in "${!_dirs[@]}"; do echo "  $d"; done
echo "Bind string: $bindlist"
echo

# Pre-flight checks INSIDE container: existence & permissions + nvidia-smi
check_cmd=""
for f in "$REF_FILE" "$KNOWN_SITES_MILLS" "$KNOWN_SITES_1000G" "$KNOWN_SITES_DBSNP" "$FASTQ1" "$FASTQ2" "$OUT_DIR"; do
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
PBRUN_LOG="${OUT_DIR}/${SAMPLE_ID}_pbrun.${SLURM_JOB_ID:-local}.log"
echo "Starting Parabricks run; logs -> ${PBRUN_LOG}"
set -x

apptainer exec --nv --bind "${bindlist}" "${SIF_IMAGE}" \
  bash -lc "echo 'Parabricks starting on $(hostname) at $(date)'; \
    nvidia-smi; \
    pbrun fq2bam \
      --ref '${REF_FILE}' \
      --in-fq '${FASTQ1}' '${FASTQ2}' \"@RG\\tID:${SAMPLE_ID}\\tLB:${SAMPLE_ID}\\tPL:ILLUMINA\\tSM:${SAMPLE_ID}\\tPU:${SAMPLE_ID}\" \
      --knownSites '${KNOWN_SITES_MILLS}' \
      --knownSites '${KNOWN_SITES_1000G}' \
      --knownSites '${KNOWN_SITES_DBSNP}' \
      --out-recal-file '${OUT_DIR}/${SAMPLE_ID}_BQSR_REPORT.txt' \
      --bwa-options='${BWA_OPTIONS}' ${PBRUN_EXTRA_OPTS} \
      --out-bam '${OUT_DIR}/${SAMPLE_ID}.bam' 2>&1 | tee '${PBRUN_LOG}'"

rc=${PIPESTATUS[0]:-0}
set +x

if [[ $rc -ne 0 ]]; then
  echo "Parabricks failed with rc=$rc. See ${PBRUN_LOG} and Parabricks stderr for details."
  echo "Saved GPU info to nvidia_host_q_${SLURM_JOB_ID:-local}.txt and nvidia_container_q_${SLURM_JOB_ID:-local}.txt"
  exit $rc
else
  echo "Parabricks finished successfully at $(date). Output written to ${OUT_DIR}/${SAMPLE_ID}.bam"
  exit 0
fi
