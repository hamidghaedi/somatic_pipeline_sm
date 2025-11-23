#!/usr/bin/env Rscript
# /global/project/hpcg6049/somatic_pipeline/somatic_pipeline_nf/bin/manifest_pairing_to_pairs.R

# Parse a WGS manifest and produce a submission plan for variant calling.
# - Infers Tumor vs Normal using tumor_descriptor: "Normal" => Normal, else Tumor.
# - Pairs per case_submitter_id:
#     * if has >=1 Normal and >=1 Tumor  -> pair EACH tumor with the FIRST normal
#     * if has Tumor only                -> tumor-only
# - Verifies BAM files exist (analysis_ready.bam, else dedup.bam, else sorted.bam, else parabricks .bam).
# - Writes a TSV plan with one row per job.
#
# -----------------------------
# OUTPUT COLUMNS (in plan TSV)
# -----------------------------
# case_submitter_id : Case key used for pairing (patients/cohorts).
# design            : "paired" if a Normal exists for the case; otherwise "tumor_only".
# tumor_id          : Tumor sample_id for this job row.
# normal_id         : Normal sample_id (first normal in the case) or NA if tumor_only.
#
# bams_tumor_ar     : Expected path to tumor analysis-ready BAM   (<sample>/<sample>.analysis_ready.bam)
# bams_tumor_dd     : Expected path to tumor de-duplicated BAM    (<sample>/<sample>.dedup.bam)
# bams_tumor_so     : Expected path to tumor sorted BAM           (<sample>/<sample>.sorted.bam)
# bams_tumor_pb     : Expected path to tumor Parabricks BAM       (<sample>/<sample>.bam)
#
# bams_normal_ar    : Expected path to normal analysis-ready BAM  (<sample>/<sample>.analysis_ready.bam)  (paired only; NA otherwise)
# bams_normal_dd    : Expected path to normal de-duplicated BAM   (<sample>/<sample>.dedup.bam)          (paired only; NA otherwise)
# bams_normal_so    : Expected path to normal sorted BAM          (<sample>/<sample>.sorted.bam)         (paired only; NA otherwise)
# bams_normal_pb    : Expected path to normal Parabricks BAM      (<sample>/<sample>.bam)                (paired only; NA otherwise)
#
# tumor_ar_exists   : TRUE if tumor analysis-ready BAM exists on disk.
# tumor_dd_exists   : TRUE if tumor dedup BAM exists on disk.
# normal_ar_exists  : TRUE if normal analysis-ready BAM exists (paired); TRUE by convention for tumor_only.
# normal_dd_exists  : TRUE if normal dedup BAM exists (paired); TRUE by convention for tumor_only.
#
# status            : Pipeline readiness status for the job:
#   - "OK"          : All required analysis-ready BAM(s) are present
#                    (paired: tumor AR & normal AR; tumor_only: tumor AR).
#   - "NEEDS_BQSR" : Analysis-ready BAM missing, but a fallback BAM exists (dedup present on either side).
#                    You can run BQSR+ApplyBQSR checkpoint to produce analysis-ready BAM(s) before calling.
#   - "MISSING_BAM": Neither AR nor DD exists where required. Upstream alignment+dedup must run first.
#
# tumor_bam_ar      : Alias of bams_tumor_ar (convenience).
# tumor_bam_dd      : Alias of bams_tumor_dd.
# normal_bam_ar     : Alias of bams_normal_ar.
# normal_bam_dd     : Alias of bams_normal_dd.
#
# NOTE:
# - Callers should prefer *analysis_ready* BAM(s).
# - If status == "NEEDS_BQSR", generate analysis-ready BAM(s) from available dedup BAM(s), then proceed.
# - If status == "MISSING_BAM", alignment stage must finish first.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage:\n",
      "  manifest_pairing_to_pairs.R <manifest.tsv> <bam_base_dir> <out_plan.tsv>\n",
      "Example:\n",
      "  manifest_pairing_to_pairs.R ",
      "/global/project/hpcg6049/somatic_pipeline/data/input/MOHCC-WGS_DNA-RUN14_2025-10-31_manifest.tsv ",
      "/global/scratch/hpc6049/fq2bam_outs ",
      "/global/scratch/hpc6049/mutect2_outs/plan.tsv\n", sep = "")
  quit(status = 2)
}

manifest_path <- args[[1]] %>% str_trim()
bam_base      <- args[[2]] %>% str_trim()
out_plan      <- args[[3]] %>% str_trim()

stopifnot(file.exists(manifest_path))
if(!dir.exists(bam_base)) dir.create(bam_base, recursive = TRUE)

# ----- helpers -----
trim <- function(x) str_replace_all(x, "^\\s+|\\s+$", "")

pick_bams <- function(sample_id, bam_base) {
  d  <- file.path(bam_base, sample_id)
  tibble::tibble(
    ar   = file.path(d, paste0(sample_id, ".analysis_ready.bam")),
    dd   = file.path(d, paste0(sample_id, ".dedup.bam")),
    so   = file.path(d, paste0(sample_id, ".sorted.bam")),
    pb   = file.path(d, paste0(sample_id, ".bam")) # optional Parabricks default
  )
}

# NA placeholder so unnest_wider always materializes columns for tumor-only rows
empty_bams <- tibble::tibble(
  ar = NA_character_,
  dd = NA_character_,
  so = NA_character_,
  pb = NA_character_
)

clean_manifest <- function(raw_df) {
  n_start <- nrow(raw_df)
  cleaned_df <- raw_df %>%
    mutate(
      case_submitter_id = str_trim(case_submitter_id),
      sample_id         = str_trim(sample_id),
      tumor_descriptor  = str_trim(tumor_descriptor)
    )
  n_before_filter <- nrow(cleaned_df)
  cleaned_df <- cleaned_df %>%
    filter(
      !is.na(case_submitter_id) & case_submitter_id != "",
      !is.na(sample_id) & sample_id != ""
    )
  n_after_filter <- nrow(cleaned_df)
  n_dropped <- n_before_filter - n_after_filter
  if (n_dropped > 0) cat(paste("MANIFEST WARN: Dropped", n_dropped, "rows due to missing/empty 'case_submitter_id' or 'sample_id'.\n"))
  n_to_std <- cleaned_df %>%
    filter(str_detect(case_submitter_id, "_")) %>%
    nrow()
  if (n_to_std > 0) {
    cat(paste(" MANIFEST INFO: Standardized", n_to_std, "case IDs that contained underscores.\n"))
  }
  cleaned_df <- cleaned_df %>%
    mutate(
      case_submitter_id = str_replace_all(case_submitter_id, "_", "-")
    )
  return(cleaned_df)
}

# ----- read & normalize manifest -----
mf <- suppressMessages(readr::read_tsv(
  manifest_path,
  col_types = readr::cols(.default = readr::col_character())
))

required_cols <- c("case_submitter_id", "sample_id", "tumor_descriptor")
missing <- setdiff(required_cols, names(mf))
if (length(missing) > 0) {
  stop("Manifest missing required columns: ", paste(missing, collapse = ", "))
}

mf <- mf %>%
  clean_manifest() %>%
  mutate(
    case_submitter_id = trim(case_submitter_id),
    sample_id         = trim(sample_id),
    tumor_descriptor  = trim(tumor_descriptor)
  ) %>%
  filter(!is.na(case_submitter_id), case_submitter_id != "",
         !is.na(sample_id), sample_id != "")

# classify role
mf <- mf %>%
  mutate(tumor_descriptor = str_trim(tumor_descriptor)) %>%
  mutate(role = case_when(
    tumor_descriptor == "Normal" ~ "NORMAL",
    tumor_descriptor == ""       ~ "UNKNOWN",
    is.na(tumor_descriptor)      ~ "UNKNOWN",
    TRUE                         ~ "TUMOR"
  ))

# ----- build pairs per case -----
cases <- mf %>%
  filter(role != "UNKNOWN") %>%
  group_by(case_submitter_id) %>%
  summarise(
    tumors  = list(sample_id[role == "TUMOR"]),
    normals = list(sample_id[role == "NORMAL"]),
    .groups = "drop"
  )
print(cases %>% as.data.frame())
cat("\n")

plan_rows <- cases %>%
  mutate(
    normal_id = purrr::map_chr(normals, ~ .x[1], .default = NA_character_)
  ) %>%
  mutate(
    design = if_else(is.na(normal_id), "tumor_only", "paired")
  ) %>%
  tidyr::unnest(tumors) %>%
  rename(tumor_id = tumors) %>%
  select(
    case_submitter_id,
    design,
    tumor_id,
    normal_id
  )

# ----- attach BAM paths & existence checks -----
plan <- plan_rows %>%
  mutate(bams        = map(tumor_id,  ~ pick_bams(.x, bam_base))) %>%
  tidyr::unnest_wider(bams, names_sep = "_tumor_") %>%  # bams_tumor_ar, bams_tumor_dd, ...
  mutate(bams_normal = if_else(design == "paired",
                               map(normal_id, ~ pick_bams(.x, bam_base)),
                               list(empty_bams))) %>%   # ensure columns exist for tumor_only rows
  tidyr::unnest_wider(bams_normal, names_sep = "_") %>%  # bams_normal_ar, bams_normal_dd, ...
  mutate(
    tumor_ar_exists   = !is.na(bams_tumor_ar)   & file.exists(bams_tumor_ar),
    tumor_dd_exists   = !is.na(bams_tumor_dd)   & file.exists(bams_tumor_dd),
    normal_ar_exists  = if_else(design=="paired",
                                !is.na(bams_normal_ar) & file.exists(bams_normal_ar),
                                TRUE),
    normal_dd_exists  = if_else(design=="paired",
                                !is.na(bams_normal_dd) & file.exists(bams_normal_dd),
                                TRUE),
    status = case_when(
      design=="paired"     & tumor_ar_exists & normal_ar_exists                  ~ "OK",
      design=="paired"     & (!tumor_ar_exists | !normal_ar_exists) &
                              (tumor_dd_exists | normal_dd_exists)              ~ "NEEDS_BQSR",
      design=="tumor_only" & tumor_ar_exists                                     ~ "OK",
      design=="tumor_only" & !tumor_ar_exists & tumor_dd_exists                  ~ "NEEDS_BQSR",
      TRUE                                                                        ~ "MISSING_BAM"
    ),

    # convenience aliases
    tumor_bam_ar  = bams_tumor_ar,
    tumor_bam_dd  = bams_tumor_dd,
    normal_bam_ar = bams_normal_ar,
    normal_bam_dd = bams_normal_dd
  ) %>%
  arrange(case_submitter_id, desc(design), tumor_id)

readr::write_tsv(plan, out_plan)

# also print a brief summary to stdout
summary_tbl <- plan %>% count(design, status, name = "n")
cat("# Plan written to:", out_plan, "\n")
print(summary_tbl)

# show problematic entries with best-available BAMs (AR -> DD -> SO -> PB)
missing_tbl <- plan %>% filter(status != "OK") %>%
  mutate(
    tumor_bam_best  = dplyr::coalesce(tumor_bam_ar, tumor_bam_dd, bams_tumor_so, bams_tumor_pb),
    normal_bam_best = dplyr::coalesce(normal_bam_ar, normal_bam_dd, bams_normal_so, bams_normal_pb)
  )

if (nrow(missing_tbl) > 0) {
  cat("\n# Entries with issues (status != OK):\n")
  print(missing_tbl %>%
          select(case_submitter_id, design, tumor_id, normal_id,
                 tumor_bam_best, normal_bam_best, status))
}
