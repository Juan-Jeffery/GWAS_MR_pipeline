#!/usr/bin/env Rscript
library(data.table)

# ─────────────── Function: Read and Format Exposure ───────────────
read_and_format_exposure_tsv <- function(file_path) {
  # Load data using fast read
  data_raw <- fread(file_path, sep = "\t", header = TRUE)
  
  required_cols <- c(
    "variant_id", "chromosome", "base_pair_location", "p_value",
    "beta", "standard_error", "effect_allele", "other_allele", "effect_allele_frequency"
  )
  
  # Validate required columns
  missing_cols <- setdiff(required_cols, colnames(data_raw))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  } 
  
  # Filter rows with critical missing values
  data <- data_raw[
    !is.na(variant_id) & variant_id != "" &
    !is.na(chromosome) &
    !is.na(base_pair_location) &
    !is.na(p_value)
  ]
  
  # Rename columns for TwoSampleMR compatibility
  setnames(data, 
           old = required_cols,
           new = c("SNP", "chr.exposure", "pos.exposure", "pval.exposure", 
                   "beta.exposure", "se.exposure", "effect_allele.exposure", 
                   "other_allele.exposure", "eaf.exposure"))
  
  return(data[, .(SNP, chr.exposure, pos.exposure, pval.exposure, 
                  beta.exposure, se.exposure, effect_allele.exposure, 
                  other_allele.exposure, eaf.exposure)])
}

# ─────────────── Execution ───────────────
exposure_path <- "/home/rd01/data/anemias_GCST90435802_gwascatalog/GCST90435802.tsv.gz"

# Load and filter by suggestive significance
exposure_df <- read_and_format_exposure_tsv(exposure_path)
filtered_df <- exposure_df[pval.exposure < 5e-4]

# Calculate summary statistics using data.table for speed
summary_stats <- filtered_df[, .(
  n_snp     = .N,
  mean_beta = mean(beta.exposure, na.rm = TRUE),
  mean_se   = mean(se.exposure, na.rm = TRUE),
  mean_eaf  = mean(eaf.exposure, na.rm = TRUE),
  mean_pval = mean(pval.exposure, na.rm = TRUE)
)]

print(summary_stats)