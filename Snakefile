import os
import pandas as pd

configfile: "config/config.yaml"

# --- Validation ---
if not config.get("manifest"):
    raise ValueError("ERROR: --config manifest=... is required.")
if not config.get("fastq_dir"):
    raise ValueError("ERROR: --config fastq_dir=... is required.")
if not config.get("batch_out"):
    config["batch_out"] = os.path.dirname(config["fastq_dir"].rstrip("/"))

config["r_script_merge"] = "workflow/scripts/merge_fastq_lanes.R"

# --- Target Aggregation ---
def get_downstream_inputs(wildcards):
    checkpoint_output = checkpoints.merge_fastq_lanes.get(**wildcards).output.groups_tsv
    try:
        df = pd.read_csv(checkpoint_output, sep="\t", comment="#", header=None, names=["prefix", "read", "n_files", "out_file"])
    except pd.errors.EmptyDataError:
        return []

    if df.empty: return []
    samples = df["prefix"].unique()
    
    targets = []
    
    # 1. Fastp/FastQC
    targets.extend([os.path.join(config["batch_out"], "fastp_outs", "reports", f"{s}.fastp.json") for s in samples])
    if config.get("run_fastqc", True):
        targets.extend([os.path.join(config["batch_out"], "fastqc_outs", f"{s}.trimmed_1_fastqc.html") for s in samples])
        targets.extend([os.path.join(config["batch_out"], "fastqc_outs", f"{s}.trimmed_2_fastqc.html") for s in samples])

    # 2. Final BAMs (triggers Align -> MarkDup -> BQSR)
    targets.extend([os.path.join(config["batch_out"], "fq2bam_outs", s, f"{s}.analysis_ready.bam") for s in samples])
    targets.extend([os.path.join(config["batch_out"], "fq2bam_outs", s, f"{s}.analysis_ready.bai") for s in samples])

    # 3. Metrics (triggers Flagstat / WGS Metrics)
    targets.extend([os.path.join(config["batch_out"], "fq2bam_outs", s, f"{s}.dedup.bam.flagstat") for s in samples])
    targets.extend([os.path.join(config["batch_out"], "fq2bam_outs", s, f"{s}.wgs_metrics.txt") for s in samples])

    return targets

# --- Includes ---
include: "workflow/rules/merge.smk"
include: "workflow/rules/fastp.smk"
include: "workflow/rules/fastqc.smk"
include: "workflow/rules/align.smk"       # <--- NEW
include: "workflow/rules/process_bam.smk" # <--- NEW
include: "workflow/rules/metrics.smk"     # <--- NEW

# --- Rule All ---
rule all:
    input:
        get_downstream_inputs
