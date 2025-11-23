rule align_bwa_mem2:
    input:
        r1 = os.path.join(config["batch_out"], "fastp_outs", "trimmed_fqs", "{sample}.trimmed_1.fastq.gz"),
        r2 = os.path.join(config["batch_out"], "fastp_outs", "trimmed_fqs", "{sample}.trimmed_2.fastq.gz")
    output:
        sorted_bam = os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "{sample}.sorted.bam")
    log:
        os.path.join(config["batch_out"], "fq2bam_outs", "{sample}", "logs", "align.log")
    params:
        sif_bwa = config["containers"]["bwamem2"],  # Points to bwa-0.7.18.sif
        sif_sam = config["containers"]["samtools"],
        idx_prefix = config["references"]["bwamem2_idx"],
        rg = r"@RG\tID:{sample}\tLB:{sample}\tPL:ILLUMINA\tSM:{sample}\tPU:{sample}"
    threads: 
        config["resources"]["align"]["cpus"]
    resources:
        mem_mb = config["resources"]["align"]["mem_mb"],
        runtime = config["resources"]["align"]["runtime"],
        cpus = config["resources"]["align"]["cpus"]
    shell:
        """
        set -euo pipefail
        module --force purge
        module load apptainer

        mkdir -p $(dirname {output.sorted_bam})
        mkdir -p $(dirname {log})

        echo "--- Starting BWA MEM Alignment ---" > {log}
        echo "BWA Container: {params.sif_bwa}" >> {log}
        
        # Verify Index exists to fail fast (BWA Classic expects .bwt, .sa, etc)
        if ! ls {params.idx_prefix}.bwt >/dev/null 2>&1; then
             echo "ERROR: BWA index not found at prefix: {params.idx_prefix}" >> {log}
             echo "Expected {params.idx_prefix}.bwt, .sa, etc." >> {log}
             exit 1
        fi

        # 1. Run BWA MEM (Classic) -> stdout
        # 2. Pipe to Samtools Sort -> Output File
        # Note: We specifically use 'bwa mem' here, not 'bwa-mem2', matching your container.
        
        apptainer exec --bind /global/project,/global/scratch {params.sif_bwa} \
            bwa mem \
            -t {threads} \
            -Y \
            -R '{params.rg}' \
            {params.idx_prefix} \
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
