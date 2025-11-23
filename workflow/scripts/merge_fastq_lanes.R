#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tools) # For file_ext
})

# --- Helper Function: Robust File Reading ---
read_manifest_safe <- function(file_path) {
  ext <- tolower(tools::file_ext(file_path))
  
  if (ext == "csv") {
    message("Detected CSV format. Reading with read_csv...")
    return(readr::read_csv(file_path, show_col_types = FALSE, lazy = FALSE))
  } else if (ext %in% c("tsv", "txt")) {
    message("Detected TSV format. Reading with read_tsv...")
    return(readr::read_tsv(file_path, show_col_types = FALSE, lazy = FALSE))
  } else {
    stop("Unsupported file extension: .", ext, "\nPlease provide a manifest ending in .csv, .tsv, or .txt")
  }
}

# --- 1. Parse Arguments (Strict) ---
args <- commandArgs(trailingOnly = TRUE)

# Helper for flags
dry_run <- any(args %in% "--dry-run")
include_undetermined <- any(args %in% "--include-undetermined")

# Filter out flags to find positional args
pos_args <- args[!stats::setNames(args %in% c("--dry-run", "--include-undetermined"), NULL)]

if (length(pos_args) < 3) {
  cat("Usage:\n",
      "  Rscript merge_fastq_lanes.R <input_fastq_dir> <output_merged_dir> <manifest_file> [--dry-run]\n\n",
      "Example:\n",
      "  Rscript merge_fastq_lanes.R ./raw_fastq ./merged_fastq ./data/manifest.csv\n",
      sep = "")
  quit(status = 1)
}

inp_dir     <- normalizePath(pos_args[1], mustWork = TRUE)
out_dir     <- pos_args[2] # Don't normalize yet, might not exist
manifest_in <- normalizePath(pos_args[3], mustWork = TRUE)

# Ensure output directory exists
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 2. Scan FASTQ Files ---
# Pattern: <prefix>_L00X_R[12]_001.fastq.gz
pat <- "^(?<prefix>.+)_L(?<lane>00[1-9])_(?<read>R[12])_001\\.fastq\\.gz$"

all_files <- list.files(inp_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE)

if (length(all_files) == 0) stop("No FASTQ files found in: ", inp_dir)

df <- tibble(path = all_files, file = basename(all_files)) %>%
  mutate(m = str_match(file, pat)) %>%
  filter(!is.na(m[,1])) %>%
  transmute(
    path,
    file,
    prefix = m[, "prefix"],
    lane   = m[, "lane"],
    read   = m[, "read"]
  )


df <- df %>% filter(!startsWith(prefix, "Undetermined"))


if (nrow(df) == 0) stop("No matching FASTQ files found after filtering.")

# --- 3. Merge Logic ---
groups <- df %>%
  arrange(prefix, read, lane) %>%
  group_by(prefix, read) %>%
  summarise(files = list(path), n_files = n(), .groups = "drop") %>%
  mutate(out_file = file.path(out_dir, paste0(prefix, "_", read, ".fastq.gz")))

# Write merge audit log
groups_tsv <- file.path(out_dir, "merge_groups.tsv")
write_lines(paste0("# input: ", inp_dir, "\n# date: ", Sys.time()), groups_tsv)
write_tsv(groups %>% select(prefix, read, n_files, out_file), groups_tsv, append = TRUE)

cat("Found ", nrow(groups), " merge groups.\n")

# Execute Merges
merge_one <- function(in_files, out_file) {
  cmd <- sprintf("cat %s > %s", paste(shQuote(in_files), collapse = " "), shQuote(out_file))
  
  if (dry_run) {
    message("[DRY-RUN] ", cmd)
  } else {
    if (file.exists(out_file) && file.info(out_file)$size > 0) {
      message("SKIP (exists): ", basename(out_file))
    } else {
      message("MERGING: ", basename(out_file))
      if (system(cmd) != 0) stop("Merge failed: ", out_file)
    }
  }
}

pwalk(groups %>% select(files, out_file), ~ merge_one(..1, ..2))

# --- 4. Update Manifest ---
cat("\n--- Updating Manifest ---\n")

lookup <- groups %>%
    distinct(prefix) %>%
    mutate(prefix_base = sub("_S[0-9]+$", "", prefix)) %>%
    # If multiple prefixes share same base (rare), keep the first (stable order above)
    group_by(prefix_base) %>%
    slice(1) %>%
    ungroup()
#
colnames(lookup) <- c("sample_id", "sample_submitter_id")
    

# 1. Read using helper
man <- read_manifest_safe(manifest_in)

# 2. Map sample_id based on exact match
if (!"sample_submitter_id" %in% names(man)) {
  warning("Manifest missing 'sample_submitter_id' column. Cannot map sample_id.")
  man_out_df <- man
} else {
  man_out_df <- man %>%
    left_join(., lookup, by="sample_submitter_id")
  
  n_matched <- sum(!is.na(man_out_df$sample_id))
  cat("Mapped", n_matched, "samples to FASTQ files based on exact 'sample_submitter_id' match.\n")

  if (n_matched == 0) {
    warning("The manifest file 'sample_submitter_id' column is not matching any file name.\nDouble check manifest and input fastq directory!", call. = FALSE)
  }
}

# 3. Write Output (Always TSV)
out_name <- sub("(\\.[^.]+)$", "", basename(manifest_in))
manifest_out_path <- file.path(dirname(manifest_in), paste0(out_name, ".with_sample_id.tsv"))

write_tsv(man_out_df, manifest_out_path)

cat("Saved updated manifest to: ", manifest_out_path, "\n")
cat("Done.\n")