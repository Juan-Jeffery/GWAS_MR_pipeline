# ─────────────── Function: Plot CMplot and filter by p-value ───────────────
plot_cmplot <- function(file_path, data, pval_threshold) {
  data_for_cmplot <- as.data.frame(data[, c("SNP", "chr.exposure", "pos.exposure", "pval.exposure")])
  colnames(data_for_cmplot) <- c("SNP", "CHR", "BP", "pvalue")
  data_for_cmplot <- na.omit(data_for_cmplot)

  CMplot(data_for_cmplot, plot.type = "m",
         LOG10 = TRUE,
         threshold = pval_threshold,
         threshold.lwd = 3,
         threshold.lty = 1,
         signal.cex = 0.2,
         chr.den.col = NULL,
         cex = 0.2,
         bin.size = 1e5,
         ylim = c(0, 10),
         file = "pdf",
         file.output = TRUE,
         width = 15,
         height = 9,
         verbose = TRUE)

  return(subset(data, pval.exposure < pval_threshold))
}

# ─────────────── Function: LD clumping ───────────────
ld_clump_exposure <- function(file_path, data, pop_n) {
  data$rsid <- data$SNP
  data$pval <- data$pval.exposure

  clumped <- ieugwasr::ld_clump(
    d = data,
    clump_kb = 1000,
    clump_r2 = 0.1,
    pop = pop_n
  )

  fwrite(clumped, file = file.path(file_path, "exposure.LD.csv"))
  return(clumped)
}

# ─────────────── Function: F-statistic filter ───────────────
filter_by_F_statistic <- function(file_path, data, exposure_n) {
  data$samplesize.exposure <- exposure_n
  data[, R2 := 2 * (beta.exposure)^2 * eaf.exposure * (1 - eaf.exposure)]
  data[, F := (samplesize.exposure - 2) * R2 / (1 - R2)]
  fwrite(data, file.path(file_path, "exposure.F.csv"))
  return(data)
}