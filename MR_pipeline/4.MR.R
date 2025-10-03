# ─────────────── Read & Harmonize ───────────────
process_outcome_data <- function(file_path, final_exposure_data, outcome_data, outcome_variant_id, exposure_name = "anemia", outcome_name = "T2D") {

  exposure_dat <- read_exposure_data(
    filename = file.path(file_path, "exposure.F.csv"),
    sep = ",",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    pval_col = "pval.exposure",
    eaf_col = "eaf.exposure",
    clump = FALSE
  )

  intersection_dat <- merge(exposure_dat, outcome_data, by.x = "SNP", by.y = outcome_variant_id)
  # variant_id, rs_id 
  write.csv(intersection_dat, file = file.path(file_path, "intersection_dat.csv"), row.names = FALSE)

  outcome_dat <- read_outcome_data(
    snps = intersection_dat$SNP,
    filename = file.path(file_path, "intersection_dat.csv"),
    sep = ",",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "p_value",
    eaf_col = "effect_allele_frequency"
  )
  write.csv(outcome_dat, file = file.path(file_path, "outcome_dat.csv"), row.names = FALSE)

  exposure_dat$exposure <- exposure_name
  outcome_dat$outcome <- outcome_name

  dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
  write.csv(dat, file = file.path(file_path, "dat.csv"), row.names = FALSE)

  dat_valid <- dat[dat$mr_keep == TRUE, ]
  write.csv(dat_valid, file = file.path(file_path, "table.SNP.csv"), row.names = FALSE)

  return(dat_valid)
}

# ─────────────── MR Analysis ───────────────
run_mr_analysis <- function(file_path, dat_valid) {
  mr_result <- mr(dat_valid, method_list = c(
    "mr_ivw",
    "mr_simple_median",
    "mr_penalised_weighted_median",
    "mr_two_sample_ml"
  ))
  mr_tab <- generate_odds_ratios(mr_result)
  write.csv(mr_tab, file = file.path(file_path, "table.MRresult.csv"), row.names = FALSE)

  heter_tab <- mr_heterogeneity(dat_valid)
  write.csv(heter_tab, file = file.path(file_path, "table.heterogeneity.csv"), row.names = FALSE)

  pleio_tab <- mr_pleiotropy_test(dat_valid)
  write.csv(pleio_tab, file = file.path(file_path, "table.pleiotropy.csv"), row.names = FALSE)

  return(mr_result)
}

# ─────────────── Plot Results ───────────────
plot_mr_results <- function(file_path, mr_result, dat_valid) {
  # 1. Scatter Plot
  pdf(file.path(file_path, "mr_scatter_plot.pdf"), width = 7, height = 7)
  print(mr_scatter_plot(mr_results = mr_result, dat = dat_valid))
  dev.off()
  # 2. Forest Plot
  res_single <- mr_singlesnp(dat_valid)
  pdf(file.path(file_path, "mr_forest_plot.pdf"), width = 7, height = 7)
  print(mr_forest_plot(res_single))
  dev.off()
  # 3. Funnel Plot
  pdf(file.path(file_path, "mr_funnel_plot.pdf"), width = 7, height = 7)
  print(mr_funnel_plot(res_single))
  dev.off()
  # 4. Leave-one-out Plot
  leave1 <- mr_leaveoneout(dat_valid)
  pdf(file.path(file_path, "mr_leaveoneout_plot.pdf"), width = 7, height = 7)
  print(mr_leaveoneout_plot(leave1))
  dev.off()
}
