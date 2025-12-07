import os

# -----------------------------------------------------------------------------
# fastp: trim merged FASTQs per sample
# -----------------------------------------------------------------------------
rule fastp_sample:
    input:
        # Merged FASTQs produced by merge_fastq_lanes.R
        r1 = lambda w: os.path.join(
            config["batch_out"], "fastq_merged", f"{w.sample}_R1.fastq.gz"
        ),
        r2 = lambda w: os.path.join(
            config["batch_out"], "fastq_merged", f"{w.sample}_R2.fastq.gz"
        )
    output:
        # Trimmed FASTQs go under fastp_outs/trimmed_fqs
        trimmed_r1 = os.path.join(
            config["batch_out"],
            "fastp_outs",
            "trimmed_fqs",
            "{sample}.trimmed_1.fastq.gz",
        ),
        trimmed_r2 = os.path.join(
            config["batch_out"],
            "fastp_outs",
            "trimmed_fqs",
            "{sample}.trimmed_2.fastq.gz",
        ),
        # JSON + HTML reports go under fastp_outs/reports
        json = os.path.join(
            config["batch_out"],
            "fastp_outs",
            "reports",
            "{sample}.fastp.json",
        ),
        html = os.path.join(
            config["batch_out"],
            "fastp_outs",
            "reports",
            "{sample}.fastp.html",
        )
    log:
        # Keep logs in fastp_outs/logs
        os.path.join(
            config["batch_out"],
            "fastp_outs",
            "logs",
            "{sample}.fastp.log",
        )
    params:
        sif = config["containers"]["fastp"],
        # Optional extra fastp args from config, if present
        extra = (
            config["fastp"]["extra"]
            if "fastp" in config and "extra" in config["fastp"]
            else ""
        )
    threads:
        config["resources"]["fastp"]["cpus"]
    resources:
        mem_mb = config["resources"]["fastp"]["mem_mb"],
        runtime = config["resources"]["fastp"]["runtime"],
        cpus = config["resources"]["fastp"]["cpus"]
    shell:
        r"""
        set -euo pipefail
        module --force purge
        module load StdEnv/2023
        module load apptainer

        # Ensure directory structure exists
        mkdir -p "$(dirname {output.trimmed_r1})"
        mkdir -p "$(dirname {output.json})"
        mkdir -p "$(dirname {log})"

        apptainer exec --bind /global/project,/global/scratch {params.sif} \
          fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.trimmed_r1} \
            --out2 {output.trimmed_r2} \
            --json {output.json} \
            --html {output.html} \
            --thread {threads} \
            {params.extra} \
          > {log} 2>&1
        """

