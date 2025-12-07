#!/bin/bash
#
# submit_test_run1.sh
#
# Submit with:
#   sbatch submit_test_run1.sh
#

#SBATCH --job-name=test_run1_smk
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --account=def-hpcg6049_cpu
#SBATCH --output=test_run1_smk_%j.out
#SBATCH --error=test_run1_smk_%j.err

set -euo pipefail

echo "Snakemake driver job started on $(hostname) at $(date)"

# ---------------------------------------------------------------------------
# 1. Make sure the driver log directory exists
# ---------------------------------------------------------------------------
mkdir -p slurm_logs

# ---------------------------------------------------------------------------
# 2. Go to the Snakemake workflow directory
# ---------------------------------------------------------------------------
cd /global/project/hpcg6049/somatic_pipeline/somatic_pipeline_sm

# ---------------------------------------------------------------------------
# 3. Activate conda environment that has snakemake installed
# ---------------------------------------------------------------------------
source /global/software/python/anaconda3-2024.06-1/etc/profile.d/conda.sh
conda activate snakemake

# Optional: print snakemake version for sanity
snakemake --version

# ---------------------------------------------------------------------------
# 4. Run Snakemake for the hardcoded test batch in Snakefile
# ---------------------------------------------------------------------------
snakemake \
  --snakefile Snakefile \
  --profile profiles/slurm \
  --cores 50 \
  --rerun-incomplete \
  --printshellcmds \
  --latency-wait 60

echo "Snakemake driver job finished at $(date)"
