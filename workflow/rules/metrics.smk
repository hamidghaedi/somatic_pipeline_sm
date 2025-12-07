rule samtools_flagstat:
    input:
        bam = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.dedup.bam"),
        bai = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.dedup.bam.bai")
    output:
        txt = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.dedup.bam.flagstat")
    log:
        os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "logs", "flagstat.log")
    params:
        sif = config["containers"]["samtools"]
    resources:
        mem_mb = 8000,
        runtime = 240,
        cpus = 2
    shell:
        """
        set -euo pipefail
        module --force purge
        module load StdEnv/2023
        module load apptainer
        apptainer exec --bind /global/project,/global/scratch {params.sif} \
            samtools flagstat {input.bam} > {output.txt} 2> {log}
        """

rule collect_wgs_metrics:
    input:
        bam = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.dedup.bam"),
        bai = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.dedup.bam.bai")
    output:
        metrics = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.wgs_metrics.txt")
    log:
        os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "logs", "wgs_metrics.log")
    params:
        sif = config["containers"]["gatk"],
        ref = config["references"]["fasta"],
        java_opts = "-Xmx56g"
    resources:
        mem_mb = 64000,
        runtime = 480,
        cpus = 8
    threads: 8
    shell:
        """
        set -euo pipefail
        module --force purge
        module load apptainer
        
        apptainer exec --bind /global/project,/global/scratch {params.sif} \
            gatk CollectWgsMetrics \
            --java-options "{params.java_opts}" \
            -I {input.bam} \
            -O {output.metrics} \
            -R {params.ref} \
            > {log} 2>&1
        """
