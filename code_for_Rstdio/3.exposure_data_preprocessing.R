# ===== 設定輸出資料夾路徑 =====
row_data_path <- "/Users/jeffery/實驗/other/anemias_ukb_e_280_CSA_ieu"
file_path <- "/Users/jeffery/實驗/other/anemia_t2d_result"
setwd(file_path)

# ===== 套件載入 =====
library(TwoSampleMR)
library(gwasglue)
library(gwasvcf)
library(VariantAnnotation)
library(MendelianRandomization)
library(dplyr)
library(ggplot2)
library(CMplot)
library(ieugwasr)

# ===== Step 1: 讀取VCF並轉成TwoSampleMR格式 ===== (20)
exposure_vcf_path <- file.path(row_data_path, "ukb-e-280_CSA.vcf.gz")
exposure_vcf <- readVcf(exposure_vcf_path)
data <- gwasvcf_to_TwoSampleMR(vcf = exposure_vcf, type = "exposure")

# ===== Step 2: 用CMplot畫圖 (看p要切多少) =====
data_for_cmplot <- data[, c("SNP", "chr.exposure", "pos.exposure", "pval.exposure")]
colnames(data_for_cmplot) <- c("SNP", "CHR", "BP", "pvalue")

CMplot(data_for_cmplot, plot.type = "m",
       LOG10 = TRUE,
       threshold = 5e-03,
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

rm(data_for_cmplot)

# ===== Step 3: 篩選p-value =====
outTab <- subset(data, pval.exposure < 5e-03)
write.csv(outTab, file = file.path(file_path, "exposure.pvalue.csv"), row.names = FALSE)

rm(data)
rm(outTab)

# ===== Step 4: 讀入pvalue篩選後資料，準備LD clump =====
exposure_dat <- read_exposure_data(
  filename = file.path(file_path, "exposure.pvalue.csv"),
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  eaf_col = "eaf.exposure",
  samplesize_col = "samplesize.exposure",
  clump = FALSE
)

exposure_dat$rsid <- exposure_dat$SNP
exposure_dat$pval <- exposure_dat$pval.exposure

exposure_dat_clumped <- ieugwasr::ld_clump(
  d = exposure_dat,
  clump_kb = 10000,
  clump_r2 = 0.001,
  pop = "SAS"
)

write.csv(exposure_dat_clumped, file = file.path(file_path, "exposure.LD.csv"), row.names = FALSE)
rm(exposure_dat)
rm(exposure_dat_clumped)

# ===== Step 5: F值過濾 =====
dat <- read.csv(file.path(file_path, "exposure.LD.csv"), header = TRUE, sep = ",", check.names = FALSE)
N <- dat[1, "samplesize.exposure"]

dat <- transform(dat,R2 = 2 * (beta.exposure)^2 * eaf.exposure * (1 - eaf.exposure))
dat <- transform(dat,F = (N - 2) * R2 / (1 - R2))

dat_filtered <- dat[dat$F > 10, ]
write.csv(dat, file.path(file_path, "exposure.F.csv"), row.names = FALSE)

rm(dat)
rm(dat_filtered)
rm(N)

# ===== Step 6: 去除混雜因素 =====
dat <- read.csv(file.path(file_path, "exposure.F.csv"), header = TRUE, sep = ",", check.names = FALSE)

snpId <- dat$SNP
y <- seq_along(snpId)
chunks <- split(snpId, ceiling(y / 20))
outTab <- data.frame()

for (i in names(chunks)) {
  confounder <- phenoscanner(
    snpquery = chunks[[i]],
    catalogue = "GWAS",
    pvalue = 1e-04,
    proxies = "None",
    r2 = 0.8,
    build = 37
  )
  outTab <- rbind(outTab, confounder$results)
}

delSnp <- unique(outTab$SNP) # 需要刪除的SNP
write.csv(outTab, file.path(file_path, "confounder.result.csv"), row.names = FALSE)

dat_filtered <- dat[!dat$SNP %in% delSnp, ]
write.csv(dat_filtered, file.path(file_path, "exposure.confounder.csv"), row.names = FALSE)

# ===== 清理 =====
rm(outTab, dat, snpId, chunks, delSnp, confounder, y, i, dat_filtered)
