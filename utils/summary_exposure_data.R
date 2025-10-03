library(data.table)
library(dplyr)

# 你的函數（不用改）
read_and_format_exposure_tsv <- function(exposure_data_raw) {
  data_raw <- fread(file.path(exposure_data_raw), sep = "\t", header = TRUE)
  required_cols <- c(
    "variant_id", "chromosome", "base_pair_location", "p_value",
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

# ✅ 1. 載入檔案（改成你的檔案路徑）
exposure_data_path <- "/home/rd01/data/anemias_GCST90435802_gwascatalog/GCST90435802.tsv.gz"
exposure_df <- read_and_format_exposure_tsv(exposure_data_path)
high_pval_df <- exposure_df[exposure_df$pval.exposure < 5e-4, ]

# ✅ 2. 計算平均值
summary_stats <- high_pval_df %>%
  summarise(
    n_snp = n(),
    mean_beta = mean(beta.exposure, na.rm = TRUE),
    mean_se = mean(se.exposure, na.rm = TRUE),
    mean_eaf = mean(eaf.exposure, na.rm = TRUE),
    mean_pval = mean(pval.exposure, na.rm = TRUE)
  )

# ✅ 3. 顯示結果
print(summary_stats)
