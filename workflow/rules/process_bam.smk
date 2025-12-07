# Mark Duplicates
rule mark_duplicates:
    input:
        sorted_bam = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.sorted.bam")
    output:
        dedup_bam = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.dedup.bam"),
        metrics = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.dedup_metrics.txt")
    log:
        os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "logs", "markdup.log")
    params:
        sif = config["containers"]["gatk"],
        java_opts = "-Xmx170g"
    threads: 8
    resources:
        mem_mb = config["resources"]["markdup"]["mem_mb"],
        runtime = config["resources"]["markdup"]["runtime"],
        cpus = config["resources"]["markdup"]["cpus"]
    shell:
        """
        set -euo pipefail
        module --force purge
        module load StdEnv/2023
        module load apptainer
        
        apptainer exec --bind /global/project,/global/scratch {params.sif} \
            gatk MarkDuplicates \
            --java-options "{params.java_opts}" \
            -I {input.sorted_bam} \
            -O {output.dedup_bam} \
            -M {output.metrics} \
            > {log} 2>&1
        """

# Index Dedup BAM (needed for BQSR)
rule index_dedup:
    input:
        dedup_bam = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.dedup.bam")
    output:
        bai = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.dedup.bam.bai")
    log:
        os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "logs", "index_dedup.log")
    params:
        sif = config["containers"]["samtools"]
    resources:
        mem_mb = 4096,
        runtime = 60,
        cpus = 1
    shell:
        """
        set -euo pipefail
        module --force purge
        module load apptainer
        apptainer exec --bind /global/project,/global/scratch {params.sif} \
            samtools index {input.dedup_bam} > {log} 2>&1
        """

# BaseRecalibrator
rule bqsr_recal_table:
    input:
        bam = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.dedup.bam"),
        bai = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.dedup.bam.bai"),
    output:
        recal_table = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.recal_table.txt")
    log:
        os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "logs", "bqsr_recal.log")
    params:
        sif = config["containers"]["gatk"],
        ref = config["references"]["fasta"],
        mills = config["references"]["known_sites"]["mills"],
        dbsnp = config["references"]["known_sites"]["dbsnp"],
        kg_indels = config["references"]["known_sites"]["kg_indels"],
        java_opts = "-Xmx170g"
    threads: 8
    resources:
        mem_mb = config["resources"]["bqsr"]["mem_mb"],
        runtime = config["resources"]["bqsr"]["runtime"],
        cpus = config["resources"]["bqsr"]["cpus"]
    shell:
        """
        set -euo pipefail
        module --force purge
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

# ApplyBQSR
rule apply_bqsr:
    input:
        bam = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.dedup.bam"),
        recal_table = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.recal_table.txt"),
        bai = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.dedup.bam.bai")
    output:
        final_bam = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.analysis_ready.bam"),
        final_bai = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.analysis_ready.bai")
    log:
        os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "logs", "apply_bqsr.log")
    params:
        sif = config["containers"]["gatk"],
        ref = config["references"]["fasta"],
        java_opts = "-Xmx170g"
    threads: 8
    resources:
        mem_mb = config["resources"]["bqsr"]["mem_mb"],
        runtime = config["resources"]["bqsr"]["runtime"],
        cpus = config["resources"]["bqsr"]["cpus"]
    shell:
        """
        set -euo pipefail
        module --force purge
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
        """
