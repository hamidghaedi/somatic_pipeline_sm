import os

# -----------------------------------------------------------------------------
# fastqc: QC on trimmed FASTQs
# -----------------------------------------------------------------------------
rule fastqc_sample:
    input:
        # Read trimmed FASTQs from fastp_outs/trimmed_fqs
        r1 = os.path.join(
            config["batch_out"],
            "fastp_outs",
            "trimmed_fqs",
            "{sample}.trimmed_1.fastq.gz",
        ),
        r2 = os.path.join(
            config["batch_out"],
            "fastp_outs",
            "trimmed_fqs",
            "{sample}.trimmed_2.fastq.gz",
        )
    output:
        # Standard FastQC naming: <basename>_fastqc.{html,zip}
        html_r1 = os.path.join(
            config["batch_out"],
            "fastqc_outs",
            "{sample}.trimmed_1_fastqc.html",
        ),
        zip_r1 = os.path.join(
            config["batch_out"],
            "fastqc_outs",
            "{sample}.trimmed_1_fastqc.zip",
        ),
        html_r2 = os.path.join(
            config["batch_out"],
            "fastqc_outs",
            "{sample}.trimmed_2_fastqc.html",
        ),
        zip_r2 = os.path.join(
            config["batch_out"],
            "fastqc_outs",
            "{sample}.trimmed_2_fastqc.zip",
        )
    log:
        os.path.join(
            config["batch_out"],
            "fastqc_outs",
            "logs",
            "{sample}.fastqc.log",
        )
    params:
        sif = config["containers"]["fastqc"],
        outdir = os.path.join(config["batch_out"], "fastqc_outs"),
        extra = (
            config["fastqc"]["extra"]
            if "fastqc" in config and "extra" in config["fastqc"]
            else ""
        )
    threads:
        config["resources"]["fastqc"]["cpus"]
    resources:
        mem_mb = config["resources"]["fastqc"]["mem_mb"],
        runtime = config["resources"]["fastqc"]["runtime"],
        cpus = config["resources"]["fastqc"]["cpus"]
    shell:
        r"""
        set -euo pipefail
        module --force purge
        module load StdEnv/2023
        module load apptainer

        mkdir -p "{params.outdir}"
        mkdir -p "$(dirname {log})"

        apptainer exec --bind /global/project,/global/scratch {params.sif} \
          fastqc \
            --threads {threads} \
            --outdir {params.outdir} \
            {input.r1} {input.r2} \
            {params.extra} \
          > {log} 2>&1
        """

