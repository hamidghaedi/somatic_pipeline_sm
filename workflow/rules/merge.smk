import os

# -------------------------------------------------------------------------
# Derive the path of the updated manifest that the R script will write.
#   <manifest_root>.with_sample_id.tsv
# -------------------------------------------------------------------------
_manifest_dir = os.path.dirname(config["manifest"])
_manifest_base = os.path.basename(config["manifest"])
_manifest_root, _manifest_ext = os.path.splitext(_manifest_base)
manifest_with_sample_id = os.path.join(
    _manifest_dir,
    f"{_manifest_root}.with_sample_id.tsv"
)

# -------------------------------------------------------------------------
# Checkpoint: merge_fastq_lanes
#   1) Run R script to add sample_id to manifest
#   2) Merge lane-level FASTQs -> per-sample R1/R2 in fastq_merged/
#   3) Write merge_groups.tsv
# -------------------------------------------------------------------------
checkpoint merge_fastq_lanes:
    input:
        fq_dir   = lambda w: config["fastq_dir"],
        manifest = lambda w: config["manifest"],
        r_script = config["r_script_merge"]
    output:
        groups_tsv = os.path.join(
            config["batch_out"],
            "fastq_merged",
            "merge_groups.tsv"
        ),
        merged_dir = directory(
            os.path.join(config["batch_out"], "fastq_merged")
        ),
        new_manifest = manifest_with_sample_id
    log:
        os.path.join(
            config["batch_out"],
            "fastq_merged",
            "merge_fastq.log"
        )
    params:
        sif = config["containers"]["tidyverse"]
    resources:
        mem_mb = config["resources"]["merge"]["mem_mb"],
        runtime = config["resources"]["merge"]["runtime"],
        cpus = config["resources"]["merge"]["cpus"]
    threads:
        config["resources"]["merge"]["cpus"]
    shell:
        r"""
        set -euo pipefail

        module --force purge
        module load StdEnv/2023
        module load apptainer

        # Ensure output locations exist
        mkdir -p "$(dirname {output.groups_tsv})"
        mkdir -p "$(dirname {output.new_manifest})"
        mkdir -p "{output.merged_dir}"

        echo "--- Starting merge_fastq_lanes on $(hostname) ---" > {log}
        echo "FASTQ dir   : {input.fq_dir}" >> {log}
        echo "Merged dir  : {output.merged_dir}" >> {log}
        echo "Manifest in : {input.manifest}" >> {log}
        echo "Threads     : {threads}" >> {log}

        # --------------------------------------------------------------
        # 1) Run R script to create the updated manifest with sample_id
        # --------------------------------------------------------------
        apptainer exec --bind /global/project,/global/scratch {params.sif} \
          Rscript {input.r_script} \
            {input.fq_dir} \
            {output.merged_dir} \
            {input.manifest} \
          >> {log} 2>&1

        echo "--- R manifest update completed ---" >> {log}
        echo "Expected updated manifest: {output.new_manifest}" >> {log}

        # --------------------------------------------------------------
        # 2) Merge lanes per sample (serial, safe)
        #    <prefix>_L00X_R1_001.fastq.gz -> <prefix>_R1.fastq.gz
        #    <prefix>_L00X_R2_001.fastq.gz -> <prefix>_R2.fastq.gz
        # --------------------------------------------------------------
        fq_dir="{input.fq_dir}"
        out_dir="{output.merged_dir}"

        cd "$fq_dir"

        ls *_R1_001.fastq.gz \
          | grep -v '^Undetermined' \
          | sed 's/_L00[0-9]_R1_001.fastq.gz$//' \
          | sort -u \
          > prefixes.txt

        echo "Found $(wc -l < prefixes.txt) sample prefixes for merging." >> {log}

        while read -r sample; do
          [ -n "$sample" ] || continue

          r1_files=$(ls "$fq_dir"/"$sample"_L00*"_R1_001.fastq.gz")
          r2_files=$(ls "$fq_dir"/"$sample"_L00*"_R2_001.fastq.gz")

          echo "$r1_files -> $out_dir/$sample"_R1.fastq.gz >> {log}
          cat $r1_files > "$out_dir"/"$sample"_R1.fastq.gz

          echo "$r2_files -> $out_dir/$sample"_R2.fastq.gz >> {log}
          cat $r2_files > "$out_dir"/"$sample"_R2.fastq.gz
        done < prefixes.txt

        echo "--- Merging completed ---" >> {log}

        # --------------------------------------------------------------
        # 3) Build merge_groups.tsv for bookkeeping
        # --------------------------------------------------------------
        {
          echo -e "prefix\tread\tn_files\tout_file"

          # R1
          for f in "{output.merged_dir}"/*_R1.fastq.gz; do
            [ -f "$f" ] || continue
            base=$(basename "$f" "_R1.fastq.gz")
            n=$(ls "{input.fq_dir}"/"$base"_L00*"_R1_001.fastq.gz" | wc -l)
            echo -e "$base\tR1\t$n\t$f"
          done

          # R2
          for f in "{output.merged_dir}"/*_R2.fastq.gz; do
            [ -f "$f" ] || continue
            base=$(basename "$f" "_R2.fastq.gz")
            n=$(ls "{input.fq_dir}"/"$base"_L00*"_R2_001.fastq.gz" | wc -l)
            echo -e "$base\tR2\t$n\t$f"
          done
        } > "{output.groups_tsv}"

        rm -f prefixes.txt

        echo "--- merge_fastq_lanes checkpoint finished successfully ---" >> {log}
        """
