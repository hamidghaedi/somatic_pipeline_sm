import os
import glob
import csv

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
# Paths and constants (derived from ROOT and your reference scripts)
# =============================================================================

MERGED_DIR       = f"{ROOT}/fastq_merged"
FASTQC_DIR       = f"{ROOT}/fastqc_raw"
FASTQC_LOGS_DIR  = f"{FASTQC_DIR}/logs"
FASTP_BASE       = f"{ROOT}/fastp_outs"
FASTP_TRIM_DIR   = f"{FASTP_BASE}/trimmed_fqs"
FASTP_REPORT_DIR = f"{FASTP_BASE}/reports"
FASTP_LOG_DIR    = f"{FASTP_BASE}/logs"
FQ2BAM_BASE      = f"{ROOT}/fq2bam_cpu_outs"

# Container images (SIFs) - exactly as on disk
FASTQC_SIF    = "/global/project/hpcg6049/somatic_pipeline/containers/NXF_SINGULARITY_CACHEDIR/fastqc-0.12.1--hdfd78af_0.sif"
FASTP_SIF     = "/global/project/hpcg6049/somatic_pipeline/containers/NXF_SINGULARITY_CACHEDIR/fastp-1.0.1--heae3180_0.sif"
BWA_SIF       = "/global/project/hpcg6049/somatic_pipeline/containers/NXF_SINGULARITY_CACHEDIR/bwa-0.7.18--h577a1d6_2.sif"
SAMTOOLS_SIF  = "/global/project/hpcg6049/somatic_pipeline/containers/NXF_SINGULARITY_CACHEDIR/samtools-1.22.1--h96c455f_0.sif"
GATK_SIF      = "/global/project/hpcg6049/somatic_pipeline/containers/NXF_SINGULARITY_CACHEDIR/gatk4-4.6.2.0--py310hdfd78af_1.sif"

# References
REF_FASTA         = "/global/project/hpcg6049/somatic_pipeline/data/references/GRCh38.d1.vd1.fa"
KNOWN_SITES_MILLS = "/global/project/hpcg6049/somatic_pipeline/data/references/Mills.d1vd1.ready.vcf.gz"
KNOWN_SITES_1KG   = "/global/project/hpcg6049/somatic_pipeline/data/references/1000G.indels.d1vd1.ready.vcf.gz"
KNOWN_SITES_DBSNP = "/global/project/hpcg6049/somatic_pipeline/data/references/dbSNP.GCF40.d1vd1.ready.vcf.gz"

# For bwa mem we follow your fq2bam CPU script: use the FASTA path as the BWA reference
BWA_REF = REF_FASTA

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

# =============================================================================
# rule all
# =============================================================================

rule all:
    input:
        # merged raw FASTQs
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
        expand(f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.wgs_metrics.txt", sample=SAMPLES)

# =============================================================================
# Pre-processing (merge lanes, FastQC, fastp)
# =============================================================================

rule merge_lanes:
    input:
        R1 = r1_lanes_for_sample,
        R2 = r2_lanes_for_sample
    output:
        R1 = f"{MERGED_DIR}/{{sample}}_R1.fastq.gz",
        R2 = f"{MERGED_DIR}/{{sample}}_R2.fastq.gz"
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
        sorted_bam = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.sorted.bam"
    log:
        f"{FQ2BAM_BASE}/{{sample}}/logs/align.log"
    params:
        sif_bwa = BWA_SIF,
        sif_sam = SAMTOOLS_SIF,
        ref     = BWA_REF,
        rg      = r"@RG\tID:{sample}\tLB:{sample}\tPL:ILLUMINA\tSM:{sample}\tPU:{sample}"
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
        dedup_bam = f"{FQ2BAM_BASE}/{{sample}}/{{sample}}.dedup.bam",
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
