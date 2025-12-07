#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tools)
})

# --------------------------------------------------------------------
# Helper: robust manifest reader (CSV / TSV / TXT)
# --------------------------------------------------------------------
read_manifest_safe <- function(file_path) {
  ext <- tolower(tools::file_ext(file_path))

  if (ext == "csv") {
    message("Detected CSV manifest. Reading with read_csv() ...")
    return(readr::read_csv(file_path, show_col_types = FALSE, lazy = FALSE))
  } else if (ext %in% c("tsv", "txt")) {
    message("Detected TSV/TXT manifest. Reading with read_tsv() ...")
    return(readr::read_tsv(file_path, show_col_types = FALSE, lazy = FALSE))
  } else {
    stop("Unsupported manifest extension: .", ext,
         "\nPlease provide a manifest ending in .csv, .tsv, or .txt")
  }
}

# --------------------------------------------------------------------
# Arguments
# We keep the same call signature as before for compatibility:
#   Rscript merge_fastq_lanes.R <fastq_dir> <merged_dir> <manifest>
#
# This script is now *only* responsible for writing the updated manifest.
# <merged_dir> is accepted but ignored.
# --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage:\n",
      "  Rscript merge_fastq_lanes.R <fastq_dir> <merged_dir> <manifest_file>\n\n",
      "NOTE: In this version the script only updates the manifest.\n",
      "      The <merged_dir> argument is accepted for compatibility\n",
      "      but not used.\n",
      sep = "")
  quit(status = 1)
}

fastq_dir   <- normalizePath(args[1], mustWork = TRUE)
# merged_dir <- args[2]  # accepted but not used
manifest_in <- normalizePath(args[length(args)], mustWork = TRUE)

message("FASTQ dir (for context): ", fastq_dir)
message("Manifest input: ", manifest_in)

# --------------------------------------------------------------------
# Read manifest and add 'sample_id' if needed
# --------------------------------------------------------------------
man <- read_manifest_safe(manifest_in)

if ("sample_id" %in% names(man)) {
  message("Column 'sample_id' already present in manifest; leaving as is.")
  man_out <- man
} else if ("sample_submitter_id" %in% names(man)) {
  message("Deriving 'sample_id' from 'sample_submitter_id'.")
  man_out <- man %>% mutate(sample_id = sample_submitter_id)
} else {
  stop("Manifest does not contain 'sample_submitter_id' or 'sample_id'.\n",
       "Cannot derive 'sample_id'.")
}

# --------------------------------------------------------------------
# Write manifest as TSV next to the original file
# --------------------------------------------------------------------
base <- sub("(\\.[^.]+)$", "", basename(manifest_in))
manifest_out_path <- file.path(
  dirname(manifest_in),
  paste0(base, ".with_sample_id.tsv")
)

readr::write_tsv(man_out, manifest_out_path)
cat("Saved updated manifest to: ", manifest_out_path, "\n")
cat("Done.\n")
