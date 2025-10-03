plot_cmplot <- function(file_path, data, pval_threshold) {
  data_for_cmplot <- as.data.frame(data[, c("SNP", "chr.exposure", "pos.exposure", "pval.exposure")])
  colnames(data_for_cmplot) <- c("SNP", "CHR", "BP", "pvalue")
  
  data_for_cmplot$pvalue <- as.numeric(data_for_cmplot$pvalue)
  data_for_cmplot <- subset(data_for_cmplot, pvalue > 0 & pvalue < 1)
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

  return(
    subset(
      transform(data, pval.exposure = as.numeric(pval.exposure)),
      !is.na(pval.exposure) & pval.exposure > 0 & pval.exposure < 1 & pval.exposure < pval_threshold
    )
  )
}

# ─────────────── Function: LD clumping ───────────────
ld_clump_exposure <- function(file_path, data, pop_n, batch_size = 50000) {

  autosome <- data$chr.exposure %in% as.character(1:22)
  data_autosome <- data[autosome, ]

  # 轉換欄位型態
  data_autosome$rsid <- as.character(data_autosome$SNP)
  data_autosome$pval <- as.numeric(as.character(data_autosome$pval.exposure))
  data_autosome$chr.exposure <- as.integer(data_autosome$chr.exposure)
  data_autosome$pos.exposure <- as.integer(data_autosome$pos.exposure)
  data_autosome$effect_allele.exposure <- as.character(data_autosome$effect_allele.exposure)
  data_autosome$other_allele.exposure <- as.character(data_autosome$other_allele.exposure)

  # 移除 NA
  data_autosome <- na.omit(data_autosome)

  # 分批執行 LD clumping
  n <- nrow(data_autosome)
  clumped_list <- list()
  for (i in seq(1, n, by = batch_size)) {
    batch <- data_autosome[i:min(i + batch_size - 1, n), ]
    clumped_batch <- ieugwasr::ld_clump(
      d = batch,
      clump_kb = 1000,
      clump_r2 = 0.01,
      pop = pop_n
    )
    clumped_list[[length(clumped_list) + 1]] <- clumped_batch
  }

  # 合併所有 batch 結果
  merged_clumped <- data.table::rbindlist(clumped_list, fill = TRUE)
  merged_clumped <- unique(merged_clumped)  # 去重複 SNP
  # 移除id行 下一個LD才是正確的id
  if ("id" %in% colnames(merged_clumped)) {
    merged_clumped[, id := NULL]
  }

  # 再對合併結果做一次全局 LD clumping
  final_clumped <- ieugwasr::ld_clump(
    d = merged_clumped,
    clump_kb = 1000,
    clump_r2 = 0.01,
    pop = pop_n
  )

  # 寫出 LD 後檔案
  fwrite(final_clumped, file = file.path(file_path, "exposure.LD.csv"))

  return(final_clumped)
}


# ─────────────── Function: F-statistic filter ───────────────
filter_by_F_statistic <- function(file_path, data, exposure_n) {
  data$samplesize.exposure <- exposure_n
  data[, R2 := 2 * (beta.exposure)^2 * eaf.exposure * (1 - eaf.exposure)]
  data[, F := (samplesize.exposure - 2) * R2 / (1 - R2)]
  fwrite(data, file.path(file_path, "exposure.F.csv"))
  return(data)
}