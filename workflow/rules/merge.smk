import os

# -------------------------------------------------------------------------
# Derive the path of the updated manifest that the R script will write.
# The R script writes:
#   dirname(manifest_in) / <basename_without_ext>.with_sample_id.tsv
# where manifest_in == config["manifest"]
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
#   - merges multi-lane FASTQs into per-sample R1/R2
#   - writes an audit table "merge_groups.tsv" into fastq_merged/
#   - writes an updated manifest next to the original manifest
# -------------------------------------------------------------------------
checkpoint merge_fastq_lanes:
    input:
        fq_dir   = lambda w: config["fastq_dir"],
        manifest = lambda w: config["manifest"],
        r_script = config["r_script_merge"]
    output:
        # Merge audit table (written by R into the merged_dir)
        groups_tsv = os.path.join(
            config["batch_out"],
            "fastq_merged",
            "merge_groups.tsv"
        ),
        # Directory of merged FASTQs (R's <out_dir> argument)
        merged_dir = directory(
            os.path.join(config["batch_out"], "fastq_merged")
        ),
        # Updated manifest written by the R script next to the original
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
        2
    shell:
        r"""
        set -euo pipefail
        module --force purge
        module load apptainer

        # Ensure output locations exist
        mkdir -p "$(dirname {output.groups_tsv})"
        mkdir -p "$(dirname {output.new_manifest})"

        echo "--- Starting Merge Step on $(hostname) ---" > {log}

        apptainer exec --bind /global/project,/global/scratch {params.sif} \
            Rscript {input.r_script} \
              {input.fq_dir} \
              {output.merged_dir} \
              {input.manifest} \
            >> {log} 2>&1
        """

