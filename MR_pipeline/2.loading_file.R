# ─────────────── read file ───────────────
read_and_format_exposure_tsv <- function(exposure_data_raw, exposure_variant_id) {
  data_raw <- fread(file.path(exposure_data_raw), sep = "\t", header = TRUE)
  # check beta & odds_ratio 
  if ("beta" %in% colnames(data_raw)) {
    target_eff_col <- "beta"
    use_or <- FALSE
  } else if ("odds_ratio" %in% colnames(data_raw)) {
    target_eff_col <- "odds_ratio"
    use_or <- TRUE
  } else {
    stop("no 'beta' and 'odds_ratio' col please check file。")
  }
  
  # check required_cols
  required_cols <- c(
    exposure_variant_id, "chromosome", "base_pair_location", "p_value",
    target_eff_col, "standard_error", "effect_allele", "other_allele", "effect_allele_frequency"
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
    !is.na(get(exposure_variant_id)) & get(exposure_variant_id) != "" &
    !is.na(data_raw$chromosome) &
    !is.na(data_raw$base_pair_location) &
    !is.na(data_raw$p_value)
  ]

  # changge OR to beta
  if (use_or) {
    message("change odds_ratio into ln(OR)...")
    data[, beta := log(get(target_eff_col))]
  } else {
    data[, beta := get(target_eff_col)]
  }
  # 如果standard_error = null
  data[, se_numeric := as.numeric(as.character(standard_error))]
  data[, se.exposure := ifelse(
    is.na(se_numeric), 
    (log(as.numeric(ci_upper)) - log(as.numeric(ci_lower))) / (2 * 1.96), 
    se_numeric
  )]
  data <- data[!is.na(se.exposure) & se.exposure > 0]

  data <- data[, .(
    SNP = get(exposure_variant_id),
    chr.exposure = chromosome,
    pos.exposure = base_pair_location,
    pval.exposure = as.numeric(p_value), # Ensure numeric
    beta.exposure = as.numeric(beta),    # Ensure numeric
    se.exposure = as.numeric(se.exposure), # USE THE CALCULATED SE, NOT standard_error
    effect_allele.exposure = effect_allele,
    other_allele.exposure = other_allele,
    eaf.exposure = as.numeric(effect_allele_frequency)
  )]
  return(data)
}

read_and_format_exposure_csv <- function(exposure_data_raw, exposure_variant_id) {
  data_raw <- fread(file.path(exposure_data_raw), sep = ",", header = TRUE)
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
read_and_format_outcome_csv <- function(outcome_data_raw) {
    data <- fread(file.path(outcome_data_raw), sep = ",")
    return(data)
}
read_and_format_outcome_vcf <- function(outcome_data_raw) {
    exposure_vcf <- readVcf(outcome_data_raw)
    data <- gwasvcf_to_TwoSampleMR(vcf = exposure_vcf, type = "outcome")
    return(data)
}

