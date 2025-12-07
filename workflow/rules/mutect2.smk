# -----------------------------------------------------------------------------
# Checkpoint: Plan Variant Calling
# -----------------------------------------------------------------------------
checkpoint plan_variant_calling:
    input:
        manifest = os.path.splitext(config["manifest"])[0] + ".with_sample_id.tsv",
        r_script = config["r_script_pairing"]
    output:
        plan = os.path.join(config["batch_out"], "mutect2_outs", "plan.tsv")
    params:
        sif = config["containers"]["tidyverse"],
        bam_base = os.path.join(config["batch_out"], "fq2bam_outs")
    log:
        os.path.join(config["batch_out"], "mutect2_outs", "logs", "pairing_plan.log")
    resources:
        mem_mb = 4096,
        runtime = 240,
        cpus = 1
    shell:
        """
        set -euo pipefail
        module --force purge
        module load StdEnv/2023
        module load apptainer

        mkdir -p $(dirname {output.plan})
        
        echo "--- Running Manifest Pairing R Script ---" > {log}
        apptainer exec --bind /global/project,/global/scratch {params.sif} \
            Rscript {input.r_script} \
            {input.manifest} \
            {params.bam_base} \
            {output.plan} \
            >> {log} 2>&1
        """

# -----------------------------------------------------------------------------
# Rule: Mutect2 Paired
# -----------------------------------------------------------------------------
rule mutect2_paired:
    input:
        tumor_bam = os.path.join(config["batch_out"], "fq2bam_outs", "{tumor}", "{tumor}.analysis_ready.bam"),
        tumor_bai = os.path.join(config["batch_out"], "fq2bam_outs", "{tumor}", "{tumor}.analysis_ready.bai"),
        normal_bam = os.path.join(config["batch_out"], "fq2bam_outs", "{normal}", "{normal}.analysis_ready.bam"),
        normal_bai = os.path.join(config["batch_out"], "fq2bam_outs", "{normal}", "{normal}.analysis_ready.bai")
    output:
        vcf_unfiltered = os.path.join(config["batch_out"], "mutect2_outs", "{tumor}__vs__{normal}", "{tumor}__vs__{normal}.unfiltered.vcf.gz"),
        vcf_filtered = os.path.join(config["batch_out"], "mutect2_outs", "{tumor}__vs__{normal}", "{tumor}__vs__{normal}.filtered.vcf.gz"),
        pileup_tumor = os.path.join(config["batch_out"], "mutect2_outs", "{tumor}__vs__{normal}", "{tumor}.pileups.table"),
        pileup_normal = os.path.join(config["batch_out"], "mutect2_outs", "{tumor}__vs__{normal}", "{normal}.pileups.table"),
        contamination = os.path.join(config["batch_out"], "mutect2_outs", "{tumor}__vs__{normal}", "{tumor}__vs__{normal}.contamination.table")
    log:
        os.path.join(config["batch_out"], "mutect2_outs", "{tumor}__vs__{normal}", "logs", "mutect2.log")
    params:
        sif = config["containers"]["gatk"],
        ref = config["references"]["fasta"],
        # FIXED: Use germline resource from config
        germline = config["references"]["germline"]
    resources:
        mem_mb = config["resources"]["mutect2"]["mem_mb"],
        runtime = config["resources"]["mutect2"]["runtime"],
        cpus = config["resources"]["mutect2"]["cpus"]
    threads:
        config["resources"]["mutect2"]["cpus"]
    shell:
        """
        set -euo pipefail
        module --force purge
        module load apptainer

        echo "--- Starting Mutect2 Paired ---" > {log}
        
        echo "[1/4] Mutect2" >> {log}
        apptainer exec --bind /global/project,/global/scratch {params.sif} \
            gatk Mutect2 \
            -R {params.ref} \
            -I {input.tumor_bam} -tumor {wildcards.tumor} \
            -I {input.normal_bam} -normal {wildcards.normal} \
            -O {output.vcf_unfiltered} \
            >> {log} 2>&1

        # 2. Pileup Summaries
        # FIXED: Use germline resource (-V, -L) instead of output VCF
        echo "[2/4] GetPileupSummaries (Tumor)" >> {log}
        apptainer exec --bind /global/project,/global/scratch {params.sif} \
            gatk GetPileupSummaries \
            -I {input.tumor_bam} \
            -V {params.germline} \
            -L {params.germline} \
            -O {output.pileup_tumor} \
            >> {log} 2>&1
        
        echo "[2/4] GetPileupSummaries (Normal)" >> {log}
        apptainer exec --bind /global/project,/global/scratch {params.sif} \
            gatk GetPileupSummaries \
            -I {input.normal_bam} \
            -V {params.germline} \
            -L {params.germline} \
            -O {output.pileup_normal} \
            >> {log} 2>&1

        # 3. Calculate Contamination
        echo "[3/4] CalculateContamination" >> {log}
        apptainer exec --bind /global/project,/global/scratch {params.sif} \
            gatk CalculateContamination \
            -I {output.pileup_tumor} \
            -matched {output.pileup_normal} \
            -O {output.contamination} \
            >> {log} 2>&1

        # 4. Filter
        echo "[4/4] FilterMutectCalls" >> {log}
        apptainer exec --bind /global/project,/global/scratch {params.sif} \
            gatk FilterMutectCalls \
            -R {params.ref} \
            -V {output.vcf_unfiltered} \
            --contamination-table {output.contamination} \
            -O {output.vcf_filtered} \
            >> {log} 2>&1
        
        echo "Done." >> {log}
        """

# -----------------------------------------------------------------------------
# Rule: Mutect2 Tumor Only
# -----------------------------------------------------------------------------
rule mutect2_tumor_only:
    input:
        tumor_bam = os.path.join(config["batch_out"], "fq2bam_outs", "{tumor}", "{tumor}.analysis_ready.bam"),
        tumor_bai = os.path.join(config["batch_out"], "fq2bam_outs", "{tumor}", "{tumor}.analysis_ready.bai")
    output:
        vcf_unfiltered = os.path.join(config["batch_out"], "mutect2_outs", "{tumor}__tumor_only", "{tumor}.unfiltered.vcf.gz"),
        vcf_filtered = os.path.join(config["batch_out"], "mutect2_outs", "{tumor}__tumor_only", "{tumor}.filtered.vcf.gz")
    log:
        os.path.join(config["batch_out"], "mutect2_outs", "{tumor}__tumor_only", "logs", "mutect2.log")
    params:
        sif = config["containers"]["gatk"],
        ref = config["references"]["fasta"]
    resources:
        mem_mb = config["resources"]["mutect2"]["mem_mb"],
        runtime = config["resources"]["mutect2"]["runtime"],
        cpus = config["resources"]["mutect2"]["cpus"]
    threads:
        config["resources"]["mutect2"]["cpus"]
    shell:
        """
        set -euo pipefail
        module --force purge
        module load apptainer

        echo "--- Starting Mutect2 Tumor-Only ---" > {log}

        # 1. Mutect2 Calling
        echo "[1/2] Mutect2" >> {log}
        apptainer exec --bind /global/project,/global/scratch {params.sif} \
            gatk Mutect2 \
            -R {params.ref} \
            -I {input.tumor_bam} -tumor {wildcards.tumor} \
            -O {output.vcf_unfiltered} \
            >> {log} 2>&1

        # 2. Filter (No contamination table)
        echo "[2/2] FilterMutectCalls" >> {log}
        apptainer exec --bind /global/project,/global/scratch {params.sif} \
            gatk FilterMutectCalls \
            -R {params.ref} \
            -V {output.vcf_unfiltered} \
            -O {output.vcf_filtered} \
            >> {log} 2>&1
            
        echo "Done." >> {log}
        """
