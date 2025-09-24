# ─────────────── read file ───────────────
read_and_format_exposure_tsv <- function(exposure_data_raw, exposure_variant_id) {
  data_raw <- fread(file.path(exposure_data_raw), sep = "\t", header = TRUE)
  required_cols <- c(
    exposure_variant_id, "chromosome", "base_pair_location", "p_value",
    "beta", "standard_error", "effect_allele", "other_allele", "effect_allele_frequency"
  )
  
  # check missing value
  missing_cols <- setdiff(required_cols, colnames(data_raw))
  if (length(missing_cols) > 0) {
    message("Missing required column(s):")
    print(missing_cols)
    message("Available columns in the input file:")
    print(colnames(data_raw))
    stop("Please ensure the input file contains all required columns.")
  } 
   
  # filter NA rows
  data <- data_raw[
    !is.na(data_raw$variant_id) & data_raw$variant_id != "" &
    !is.na(data_raw$chromosome) &
    !is.na(data_raw$base_pair_location) &
    !is.na(data_raw$p_value)
  ]
  data <- data[, .(
  SNP = variant_id,
  chr.exposure = chromosome,
  pos.exposure = base_pair_location,
  pval.exposure = p_value,
  beta.exposure = beta,
  se.exposure = standard_error,
  effect_allele.exposure = effect_allele,
  other_allele.exposure = other_allele,
  eaf.exposure = effect_allele_frequency
  )]
  return(data)
}

read_and_format_exposure_vcf <- function(exposure_data_raw) {
  exposure_vcf <- readVcf(file.path(exposure_data_raw))
  data <- gwasvcf_to_TwoSampleMR(vcf = exposure_vcf, type = "exposure")
  return(data)
}

read_and_format_outcome_tsv <- function(outcome_data_raw) {
    data <- fread(file.path(outcome_data_raw), sep = "\t")
    return(data)
}

read_and_format_outcome_vcf <- function(outcome_data_raw) {
    exposure_vcf <- readVcf(outcome_data_raw)
    data <- gwasvcf_to_TwoSampleMR(vcf = exposure_vcf, type = "outcome")
    return(data)
}

