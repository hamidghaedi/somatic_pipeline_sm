import os
import glob
import csv
import datetime

# =============================================================================
# Manifest-driven configuration
# =============================================================================

# Path to the manifest CSV for this test run
MANIFEST_CSV = "/global/scratch/hpc6049/TEST_DNA-RUN1_2025-11-16_NNNNNNNNNN/data/test_manifest.csv"

# Batch ID to run (must match the manifest 'batch_id' column)
TARGET_BATCH_ID = "TEST_DNA-RUN1_2025-11-16_NNNNNNNNNN"


def load_manifest(manifest_csv, batch_id):
    """
    Load manifest rows for a given batch_id and build a per-sample dictionary.

    Key = sample_submitter_id
    Value = dict with metadata that we may use later.
    """
    samples = {}

    # test_manifest is a CSV (comma-separated)
    with open(manifest_csv, newline="") as f:
        reader = csv.DictReader(f, delimiter=",")
        for row in reader:
            if row.get("batch_id", "") != batch_id:
                continue

            sample_id = row.get("sample_submitter_id", "").strip()
            if not sample_id:
                continue

            lanes_str = row.get("lanes", "")
            lanes = [x.strip() for x in lanes_str.split(";") if x.strip()]

            samples[sample_id] = {
                "project_id":          row.get("project_id", ""),
                "batch_id":            row.get("batch_id", ""),
                "case_id":             row.get("case_submitter_id", ""),
                "sample_id":           sample_id,
                "tumor_descriptor":    row.get("tumor_descriptor", ""),
                "analyte_type":        row.get("analyte_type", ""),
                "fastq_dir":           row.get("fastq_local_dir", ""),
                "lanes":               lanes,
                "i7_index":            row.get("i7_index", ""),
                "i5_index":            row.get("i5_index", ""),
                "library_layout":      row.get("library_layout", ""),
                "specimen_type":       row.get("specimen_type", ""),
                "tissue_preservation": row.get("tissue_preservation_method", ""),
                # optional fields for MultiQC header
                "flowcell_id":         row.get("flowcell_id", row.get("flowcell", "")),
                "run_date":            row.get("run_date", row.get("sequencing_run_date", "")),
            }

    if not samples:
        raise ValueError(
            f"No samples found in manifest {manifest_csv} for batch_id={batch_id}"
        )

    return samples


MANIFEST_SAMPLES = load_manifest(MANIFEST_CSV, TARGET_BATCH_ID)
SAMPLES = sorted(MANIFEST_SAMPLES.keys())

# Derive ROOT from fastq_local_dir (parent of ".../fastq")
_any_sample = next(iter(MANIFEST_SAMPLES.values()))
FASTQ_DIR = _any_sample["fastq_dir"]
if not FASTQ_DIR:
    raise ValueError("fastq_local_dir is empty in manifest for selected batch.")

ROOT = os.path.dirname(FASTQ_DIR)

# =============================================================================
# Paths and constants
# =============================================================================

# Where this workflow (Snakefile + config dir) lives
WORKFLOW_DIR = workflow.basedir

# I/O layout relative to ROOT (per-batch scratch area)
MERGED_DIR       = f"{ROOT}/fastq_merged"
FASTQC_DIR       = f"{ROOT}/fastqc_raw"
FASTQC_LOGS_DIR  = f"{FASTQC_DIR}/logs"
FASTP_BASE       = f"{ROOT}/fastp_outs"
FASTP_TRIM_DIR   = f"{FASTP_BASE}/trimmed_fqs"
FASTP_REPORT_DIR = f"{FASTP_BASE}/reports"
FASTP_LOG_DIR    = f"{FASTP_BASE}/logs"
FQ2BAM_BASE      = f"{ROOT}/fq2bam_cpu_outs"
MULTIQC_OUTDIR   = f"{ROOT}/multiqc_out"

# Container images (SIFs)
CONTAINER_BASE = "/global/project/hpcg6049/somatic_pipeline/containers/NXF_SINGULARITY_CACHEDIR"

FASTQC_SIF   = f"{CONTAINER_BASE}/fastqc-0.12.1--hdfd78af_0.sif"
FASTP_SIF    = f"{CONTAINER_BASE}/fastp-1.0.1--heae3180_0.sif"
BWA_SIF      = f"{CONTAINER_BASE}/bwa-0.7.18--h577a1d6_2.sif"
SAMTOOLS_SIF = f"{CONTAINER_BASE}/samtools-1.22.1--h96c455f_0.sif"
GATK_SIF     = f"{CONTAINER_BASE}/gatk4-4.6.2.0--py310hdfd78af_1.sif"
MULTIQC_SIF  = f"{CONTAINER_BASE}/multiqc-1.32--pyhdfd78af_0.sif"

# References
REF_BASE = "/global/project/hpcg6049/somatic_pipeline/data/references"

REF_FASTA         = f"{REF_BASE}/GRCh38.d1.vd1.fa"
KNOWN_SITES_MILLS = f"{REF_BASE}/Mills.d1vd1.ready.vcf.gz"
KNOWN_SITES_1KG   = f"{REF_BASE}/1000G.indels.d1vd1.ready.vcf.gz"
KNOWN_SITES_DBSNP = f"{REF_BASE}/dbSNP.GCF40.d1vd1.ready.vcf.gz"

# MultiQC assets in the workflow repo
BASE_MULTIQC_CONFIG_TEMPLATE = f"{WORKFLOW_DIR}/config/multiqc_config.base.yaml"
LOGO_PATH = f"{WORKFLOW_DIR}/config/mohcc_logo_tiny.jpg"
CSS_PATH  = f"{WORKFLOW_DIR}/config/mohcc_multiqc.css"

# Batch-specific MultiQC config and report paths
MULTIQC_BATCH_CONFIG = f"{MULTIQC_OUTDIR}/multiqc_{TARGET_BATCH_ID}.yaml"
MULTIQC_REPORT       = f"{MULTIQC_OUTDIR}/{TARGET_BATCH_ID}_multiqc_report.html"

READS = ["R1", "R2"]

# =============================================================================
# Helper functions
# =============================================================================

def r1_lanes_for_sample(wc):
    """
    Collect all R1 FASTQs for a given sample using fastq_dir from the manifest.
    Pattern:
      <fastq_dir>/{sample}_S*_L*_R1_001.fastq.gz
    """
    info = MANIFEST_SAMPLES[wc.sample]
    fq_dir = info["fastq_dir"]
    pattern = os.path.join(fq_dir, f"{wc.sample}_S*_L*_R1_001.fastq.gz")
    files = sorted(glob.glob(pattern))
    if not files:
        raise ValueError(f"No R1 FASTQs found for sample {wc.sample} with pattern {pattern}")
    return files


def r2_lanes_for_sample(wc):
    """
    Collect all R2 FASTQs for a given sample using fastq_dir from the manifest.
    Pattern:
      <fastq_dir>/{sample}_S*_L*_R2_001.fastq.gz
    """
    info = MANIFEST_SAMPLES[wc.sample]
    fq_dir = info["fastq_dir"]
    pattern = os.path.join(fq_dir, f"{wc.sample}_S*_L*_R2_001.fastq.gz")
    files = sorted(glob.glob(pattern))
    if not files:
        raise ValueError(f"No R2 FASTQs found for sample {wc.sample} with pattern {pattern}")
    return files


# -------------------------------------------------------------------------
# Dynamic runtime for alignment: each 20 GiB of R1 ~= 12 hours
# -------------------------------------------------------------------------

BLOCK_GIB       = 20
HOURS_PER_BLOCK = 12
MIN_HOURS       = 12  # minimum walltime per sample


def compute_align_hours_for_sample(sample):
    """
    Compute alignment walltime in hours for a given sample based on
    total R1 lane size.
    """
    info = MANIFEST_SAMPLES[sample]
    fq_dir = info["fastq_dir"]
    pattern = os.path.join(fq_dir, f"{sample}_S*_L*_R1_001.fastq.gz")
    files = glob.glob(pattern)

    if not files:
        return MIN_HOURS

    total_bytes = 0
    for f in files:
        try:
            total_bytes += os.path.getsize(f)
        except OSError:
            continue

    if total_bytes <= 0:
        return MIN_HOURS

    size_gib = total_bytes / (1024.0 ** 3)
    blocks = int(-(-size_gib // BLOCK_GIB))  # ceil(size_gib / BLOCK_GIB)
    if blocks < 1:
        blocks = 1
    hours = blocks * HOURS_PER_BLOCK
    if hours < MIN_HOURS:
        hours = MIN_HOURS
    return int(hours)


# runtime is in MINUTES for the slurm profile (runtime: 240 = 4h)
ALIGN_RUNTIME_MIN = {
    sample: compute_align_hours_for_sample(sample) * 60
    for sample in SAMPLES
}

# -------------------------------------------------------------------------
# Helper: write a batch-specific MultiQC config YAML
# -------------------------------------------------------------------------

import yaml  # requires PyYAML (comes with snakemake env)

def write_multiqc_config(out_path):
    """
    Create a batch-specific MultiQC config based on a base template
    and manifest metadata.
    """
    # Use the first sample as representative for batch-level fields
    info0 = next(iter(MANIFEST_SAMPLES.values()))

    project_id = info0.get("project_id") or "UNKNOWN"
    batch_id   = info0.get("batch_id") or TARGET_BATCH_ID
    flowcell   = info0.get("flowcell_id") or "UNKNOWN"
    run_date   = info0.get("run_date") or "UNKNOWN"
    report_date = datetime.date.today().isoformat()

    # Load the base config (if missing, raise a clear error)
    if not os.path.exists(BASE_MULTIQC_CONFIG_TEMPLATE):
        raise FileNotFoundError(
            f"Base MultiQC config not found at {BASE_MULTIQC_CONFIG_TEMPLATE}"
        )

    with open(BASE_MULTIQC_CONFIG_TEMPLATE) as f:
        cfg = yaml.safe_load(f) or {}

    # Override / inject the dynamic bits
    cfg.setdefault("title", "MOHCC-WGTS Quality Control Report")
    cfg["subtitle"] = f"Batch ID: {batch_id}"

    cfg["report_header_info"] = [
        {"Project": project_id},
        {"Batch ID": batch_id},
        {"Flowcell ID": flowcell},
        {"Run Date": run_date},
        {"Report Date": report_date},
    ]

    cfg["custom_logo"] = LOGO_PATH
    cfg.setdefault("custom_logo_title", "Marathon of Hope Cancer Centres Network")
    cfg["custom_css_files"] = [CSS_PATH]

    # Ensure output dir exists
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    with open(out_path, "w") as out_f:
        yaml.safe_dump(cfg, out_f, sort_keys=False)


# =============================================================================
# rule all
# =============================================================================

rule all:
    input:
        # merged raw FASTQs (temp, but still part of DAG)
        expand(f"{MERGED_DIR}/{{sample}}_R1.fastq.gz", sample=SAMPLES),
        expand(f"{MERGED_DIR}/{{sample}}_R2.fastq.gz", sample=SAMPLES),

        # FastQC on merged raw reads
        expand(f"{FASTQC_DIR}/{{sample}}_{{read}}_fastqc.html",
               sample=SAMPLES, read=READS),
        expand(f"{FASTQC_DIR}/{{sample}}_{{read}}_fastqc.zip",
               sample=SAMPLES, read=READS),

        # fastp-trimmed FASTQs
        expand(f"{FASTP_TRIM_DIR}/{{sample}}.trimmed_1.fastq.gz", sample=SAMPLES),
        expand(f"{FASTP_TRIM_DIR}/{{sample}}.trimmed_2.fastq.gz", sample=SAMPLES),

        # fq2bam: final analysis-ready BAM + index
        expand(f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.analysis_ready.bam", sample=SAMPLES),
        expand(f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.analysis_ready.bam.bai", sample=SAMPLES),

        # QC metrics
        expand(f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.dedup.bam.flagstat", sample=SAMPLES),
        expand(f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.wgs_metrics.txt", sample=SAMPLES),

        # MultiQC report
        MULTIQC_REPORT

# =============================================================================
# Pre-processing (merge lanes, FastQC, fastp)
# =============================================================================

rule merge_lanes:
    input:
        R1 = r1_lanes_for_sample,
        R2 = r2_lanes_for_sample
    output:
        R1 = temp(f"{MERGED_DIR}/{{sample}}_R1.fastq.gz"),
        R2 = temp(f"{MERGED_DIR}/{{sample}}_R2.fastq.gz")
    threads: 2
    resources:
        mem_mb = 4000
    message:
        "Merging lanes for sample {wildcards.sample}"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{MERGED_DIR}"

        # R1
        cat {input.R1} > {output.R1}

        # R2
        cat {input.R2} > {output.R2}
        """


rule fastqc_raw:
    input:
        fq = f"{MERGED_DIR}/{{sample}}_{{read}}.fastq.gz"
    output:
        html = f"{FASTQC_DIR}/{{sample}}_{{read}}_fastqc.html",
        zip  = f"{FASTQC_DIR}/{{sample}}_{{read}}_fastqc.zip"
    log:
        f"{FASTQC_LOGS_DIR}/{{sample}}_{{read}}.fastqc.log"
    threads: 8
    resources:
        mem_mb = 8000
    params:
        sif    = FASTQC_SIF,
        outdir = FASTQC_DIR
    message:
        "Running FastQC on {input.fq}"
    shell:
        r"""
        set -euo pipefail
        module --force purge
        module load StdEnv/2023
        module load apptainer

        mkdir -p "{params.outdir}"
        mkdir -p "{FASTQC_LOGS_DIR}"

        apptainer exec --bind /global/project,/global/scratch {params.sif} \
          fastqc \
            --threads {threads} \
            --outdir {params.outdir} \
            --nogroup \
            {input.fq} \
          > {log} 2>&1
        """


rule fastp_sample:
    input:
        r1 = lambda w: f"{MERGED_DIR}/{w.sample}_R1.fastq.gz",
        r2 = lambda w: f"{MERGED_DIR}/{w.sample}_R2.fastq.gz"
    output:
        trimmed_r1 = f"{FASTP_TRIM_DIR}/{{sample}}.trimmed_1.fastq.gz",
        trimmed_r2 = f"{FASTP_TRIM_DIR}/{{sample}}.trimmed_2.fastq.gz",
        json       = f"{FASTP_REPORT_DIR}/{{sample}}.fastp.json",
        html       = f"{FASTP_REPORT_DIR}/{{sample}}.fastp.html"
    log:
        f"{FASTP_LOG_DIR}/{{sample}}.fastp.log"
    params:
        sif = FASTP_SIF
    threads: 16
    resources:
        mem_mb = 16000
    message:
        "Running fastp on merged FASTQs for sample {wildcards.sample}"
    shell:
        r"""
        set -euo pipefail
        module --force purge
        module load StdEnv/2023
        module load apptainer

        mkdir -p "{FASTP_TRIM_DIR}"
        mkdir -p "{FASTP_REPORT_DIR}"
        mkdir -p "{FASTP_LOG_DIR}"

        apptainer exec --bind /global/project,/global/scratch {params.sif} \
          fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.trimmed_r1} \
            --out2 {output.trimmed_r2} \
            --json {output.json} \
            --html {output.html} \
            --thread {threads} \
          > {log} 2>&1
        """

# =============================================================================
# Alignment + fq2bam (BWA -> sorted BAM -> MarkDuplicates -> BQSR)
# =============================================================================

rule align_bwa_mem:
    input:
        r1 = lambda w: f"{FASTP_TRIM_DIR}/{w.sample}.trimmed_1.fastq.gz",
        r2 = lambda w: f"{FASTP_TRIM_DIR}/{w.sample}.trimmed_2.fastq.gz"
    output:
        sorted_bam = temp(f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.sorted.bam")
    log:
        f"{FQ2BAM_BASE}/{{sample}}/logs/align.log"
    params:
        sif_bwa = BWA_SIF,
        sif_sam = SAMTOOLS_SIF,
        ref     = REF_FASTA,
        rg      = lambda w: f"@RG\\tID:{w.sample}\\tLB:{w.sample}\\tPL:ILLUMINA\\tSM:{w.sample}\\tPU:{w.sample}"
    threads: 32
    resources:
        mem_mb  = 180000,
        runtime = lambda w, attempt: ALIGN_RUNTIME_MIN[w.sample]
    message:
        "Aligning {wildcards.sample} with BWA-MEM and sorting with samtools"
    shell:
        r"""
        set -euo pipefail
        module --force purge
        module load StdEnv/2023
        module load apptainer

        mkdir -p "$(dirname {output.sorted_bam})"
        mkdir -p "$(dirname {log})"

        echo "--- Starting BWA MEM Alignment ---" > {log}
        echo "BWA Container: {params.sif_bwa}" >> {log}
        echo "Reference: {params.ref}" >> {log}
        echo "Requested runtime (minutes): {resources.runtime}" >> {log}

        # Fail fast if BWA index is missing
        if ! ls {params.ref}.bwt >/dev/null 2>&1; then
             echo "ERROR: BWA index not found: {params.ref}.bwt" >> {log}
             exit 1
        fi

        apptainer exec --bind /global/project,/global/scratch {params.sif_bwa} \
            bwa mem \
              -t {threads} \
              -Y \
              -R '{params.rg}' \
              {params.ref} \
              {input.r1} {input.r2} \
              2>> {log} \
        | apptainer exec --bind /global/project,/global/scratch {params.sif_sam} \
            samtools sort \
              -@ {threads} \
              -o {output.sorted_bam} \
              - \
              2>> {log}

        if [[ ! -s {output.sorted_bam} ]]; then
            echo "ERROR: Output BAM is missing or empty." >> {log}
            exit 1
        fi

        echo "Alignment completed successfully." >> {log}
        """


rule mark_duplicates:
    input:
        sorted_bam = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.sorted.bam"
    output:
        dedup_bam = temp(f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.dedup.bam"),
        metrics   = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.dedup_metrics.txt"
    log:
        f"{FQ2BAM_BASE}/{{sample}}/logs/markdup.log"
    params:
        sif       = GATK_SIF,
        java_opts = "-Xmx56g"
    threads: 8
    resources:
        mem_mb  = 64000,
        runtime = 1440
    message:
        "Marking duplicates for {wildcards.sample} with GATK MarkDuplicates"
    shell:
        r"""
        set -euo pipefail
        module --force purge
        module load StdEnv/2023
        module load apptainer

        mkdir -p "$(dirname {output.dedup_bam})"
        mkdir -p "$(dirname {log})"

        apptainer exec --bind /global/project,/global/scratch {params.sif} \
          gatk MarkDuplicates \
            --java-options "{params.java_opts}" \
            -I {input.sorted_bam} \
            -O {output.dedup_bam} \
            -M {output.metrics} \
          > {log} 2>&1
        """


rule index_dedup:
    input:
        dedup_bam = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.dedup.bam"
    output:
        bai = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.dedup.bam.bai"
    log:
        f"{FQ2BAM_BASE}/{{sample}}/logs/index_dedup.log"
    params:
        sif = SAMTOOLS_SIF
    threads: 1
    resources:
        mem_mb  = 4096,
        runtime = 60
    message:
        "Indexing deduplicated BAM for {wildcards.sample}"
    shell:
        r"""
        set -euo pipefail
        module --force purge
        module load StdEnv/2023
        module load apptainer

        apptainer exec --bind /global/project,/global/scratch {params.sif} \
          samtools index {input.dedup_bam} \
          > {log} 2>&1
        """


rule bqsr_recal_table:
    input:
        bam = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.dedup.bam",
        bai = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.dedup.bam.bai"
    output:
        recal_table = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.recal_table.txt"
    log:
        f"{FQ2BAM_BASE}/{{sample}}/logs/bqsr_recal.log"
    params:
        sif       = GATK_SIF,
        ref       = REF_FASTA,
        mills     = KNOWN_SITES_MILLS,
        dbsnp     = KNOWN_SITES_DBSNP,
        kg_indels = KNOWN_SITES_1KG,
        java_opts = "-Xmx56g"
    threads: 8
    resources:
        mem_mb  = 64000,
        runtime = 1380
    message:
        "Running GATK BaseRecalibrator for {wildcards.sample}"
    shell:
        r"""
        set -euo pipefail
        module --force purge
        module load StdEnv/2023
        module load apptainer

        apptainer exec --bind /global/project,/global/scratch {params.sif} \
          gatk BaseRecalibrator \
            --java-options "{params.java_opts}" \
            -R {params.ref} \
            -I {input.bam} \
            --known-sites {params.mills} \
            --known-sites {params.kg_indels} \
            --known-sites {params.dbsnp} \
            -O {output.recal_table} \
          > {log} 2>&1
        """


rule apply_bqsr:
    input:
        bam         = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.dedup.bam",
        recal_table = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.recal_table.txt",
        bai         = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.dedup.bam.bai"
    output:
        final_bam = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.analysis_ready.bam",
        final_bai = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.analysis_ready.bam.bai"
    log:
        f"{FQ2BAM_BASE}/{{sample}}/logs/apply_bqsr.log"
    params:
        sif       = GATK_SIF,
        ref       = REF_FASTA,
        java_opts = "-Xmx56g"
    threads: 8
    resources:
        mem_mb  = 64000,
        runtime = 1380
    message:
        "Applying BQSR to {wildcards.sample} and producing analysis-ready BAM"
    shell:
        r"""
        set -euo pipefail
        module --force purge
        module load StdEnv/2023
        module load apptainer

        apptainer exec --bind /global/project,/global/scratch {params.sif} \
          gatk ApplyBQSR \
            --java-options "{params.java_opts}" \
            -R {params.ref} \
            -I {input.bam} \
            --bqsr-recal-file {input.recal_table} \
            -O {output.final_bam} \
            --create-output-bam-index true \
          > {log} 2>&1

        # Ensure index exists at expected path
        if [[ ! -s {output.final_bai} ]]; then
          apptainer exec --bind /global/project,/global/scratch {SAMTOOLS_SIF} \
            samtools index {output.final_bam}
        fi
        """

# =============================================================================
# QC metrics: samtools flagstat + GATK CollectWgsMetrics
# =============================================================================

rule samtools_flagstat:
    input:
        bam = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.dedup.bam",
        bai = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.dedup.bam.bai"
    output:
        txt = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.dedup.bam.flagstat"
    log:
        f"{FQ2BAM_BASE}/{{sample}}/logs/flagstat.log"
    params:
        sif = SAMTOOLS_SIF
    threads: 2
    resources:
        mem_mb  = 8000,
        runtime = 480
    message:
        "Running samtools flagstat on {wildcards.sample}.dedup.bam"
    shell:
        r"""
        set -euo pipefail
        module --force purge
        module load StdEnv/2023
        module load apptainer

        apptainer exec --bind /global/project,/global/scratch {params.sif} \
          samtools flagstat {input.bam} \
          > {output.txt} 2> {log}
        """


rule collect_wgs_metrics:
    input:
        bam = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.dedup.bam",
        bai = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.dedup.bam.bai"
    output:
        metrics = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.wgs_metrics.txt"
    log:
        f"{FQ2BAM_BASE}/{{sample}}/logs/wgs_metrics.log"
    params:
        sif       = GATK_SIF,
        ref       = REF_FASTA,
        java_opts = "-Xmx56g"
    threads: 16
    resources:
        mem_mb  = 64000,
        runtime = 1439
    message:
        "Collecting WGS metrics for {wildcards.sample}.dedup.bam"
    shell:
        r"""
        set -euo pipefail
        module --force purge
        module load StdEnv/2023
        module load apptainer

        apptainer exec --bind /global/project,/global/scratch {params.sif} \
          gatk CollectWgsMetrics \
            --java-options "{params.java_opts}" \
            -I {input.bam} \
            -O {output.metrics} \
            -R {params.ref} \
          > {log} 2>&1
        """

# =============================================================================
# MultiQC: config + report
# =============================================================================

rule multiqc_config:
    output:
        MULTIQC_BATCH_CONFIG
    run:
        write_multiqc_config(output[0])


rule multiqc:
    input:
        expand(f"{FASTP_REPORT_DIR}/{{sample}}.fastp.json", sample=SAMPLES),
        expand(f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.dedup.bam.flagstat", sample=SAMPLES),
        expand(f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.wgs_metrics.txt", sample=SAMPLES),
        MULTIQC_BATCH_CONFIG
    output:
        html = MULTIQC_REPORT
    log:
        f"{MULTIQC_OUTDIR}/multiqc.log"
    params:
        sif         = MULTIQC_SIF,
        fastp_dir   = FASTP_BASE,
        align_dir   = FQ2BAM_BASE,
        outdir      = MULTIQC_OUTDIR,
        batch_name  = TARGET_BATCH_ID,
        config      = MULTIQC_BATCH_CONFIG,
        # MultiQC: -n <name> -> <name>.html and <name>_data/
        report_name = f"{TARGET_BATCH_ID}_multiqc_report",
        data_dir    = f"{MULTIQC_OUTDIR}/{TARGET_BATCH_ID}_multiqc_report_data"
    threads: 2
    resources:
        mem_mb  = 16000,
        runtime = 60
    message:
        "Running MultiQC for batch {params.batch_name}"
    shell:
        r"""
        set -euo pipefail
        module --force purge
        module load StdEnv/2023
        module load apptainer

        mkdir -p "{params.outdir}"

        echo "--- MultiQC configuration ---" > {log}
        echo "Container: {params.sif}"       >> {log}
        echo "FASTP Dir: {params.fastp_dir}" >> {log}
        echo "Align Dir: {params.align_dir}" >> {log}
        echo "Config:    {params.config}"    >> {log}
        echo "Report:    {output.html}"      >> {log}
        echo ""                              >> {log}

        # Clean any previous report so MultiQC doesn't create *_1.html
        rm -f  "{output.html}"
        rm -rf "{params.data_dir}"

        apptainer exec \
          --bind /global/project,/global/scratch \
          {params.sif} \
          multiqc \
            -f \
            -n "{params.report_name}" \
            -o "{params.outdir}" \
            -c "{params.config}" \
            "{params.fastp_dir}" \
            "{params.align_dir}" \
          >> {log} 2>&1

        # Sanity check: MultiQC must produce the expected HTML
        if [[ ! -f "{output.html}" ]]; then
          echo "ERROR: MultiQC did not produce expected HTML: {output.html}" >> {log}
          exit 1
        fi

        echo "MultiQC completed successfully for batch {params.batch_name}" >> {log}
        """