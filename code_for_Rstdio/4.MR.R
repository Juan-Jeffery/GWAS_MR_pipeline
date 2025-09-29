##### 設定路徑 #####
outcome_raw_path <- "/Users/jeffery/實驗/other/t2d_GCST90093001_Nasian_gwascatalog"
file_path <- "/Users/jeffery/實驗/other/anemia_t2d_result"
setwd(file_path)

##### 載入套件 #####
library(TwoSampleMR)
library(gwasglue)
library(VariantAnnotation)
library(gwasvcf)
library(dplyr)
library(ggplot2)
library(data.table)
library(R.utils)


##### 讀取 Outcome 檔案 (.tsv) #####
outcome_dat_raw <- fread(file.path(outcome_raw_path, "GCST90093109_buildGRCh37.tsv.gz"), sep = "\t")

##### 讀取 Exposure 檔案 (已過濾後的 CSV) #####
#exposure_raw_dat <- read.csv(file.path(file_path, "exposure.confounder.csv"), header = TRUE)
exposure_raw_dat <- read.csv(file.path(file_path, "exposure.F.csv"), header = TRUE, sep = ",", check.names = FALSE)

colnames(exposure_raw_dat)

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

##### 篩選共同 SNP 並儲存交集 #####
intersection_dat <- merge(exposure_dat, outcome_dat_raw, by.x = "SNP", by.y = "variant_id")
write.csv(intersection_dat, file = file.path(file_path, "intersection_dat.csv"), row.names = FALSE)

##### 建立 Outcome 物件 #####
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

##### 設定 Exposure / Outcome 名稱 #####
exposure_dat$exposure <- "anemia"
outcome_dat$outcome   <- "T2D"

##### 對齊 exposure 與 outcome 資料 #####
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
write.csv(dat, file = file.path(file_path, "dat.csv"), row.names = FALSE)

##### 保留有效 SNP 進行 MR #####
dat_valid <- dat[dat$mr_keep == TRUE, ]
write.csv(dat_valid, file = file.path(file_path, "table.SNP.csv"), row.names = FALSE)

##### MR 分析 #####
mrResult <- mr(dat_valid, method_list = c(
  "mr_ivw", "mr_simple_median",
  "mr_penalised_weighted_median", "mr_two_sample_ml"
))
mrTab <- generate_odds_ratios(mrResult)
write.csv(mrTab, file = file.path(file_path, "table.MRresult.csv"), row.names = FALSE)

##### 異質性檢驗 #####
heterTab <- mr_heterogeneity(dat_valid)
write.csv(heterTab, file = file.path(file_path, "table.heterogeneity.csv"), row.names = FALSE)

##### 多效性檢驗（Egger 截距）#####
pleioTab <- mr_pleiotropy_test(dat_valid)
write.csv(pleioTab, file = file.path(file_path, "table.pleiotropy.csv"), row.names = FALSE)

##### 作圖區塊 #####
# 1. 散點圖
pdf(file.path(file_path, "mr_scatter_plot.pdf"), width = 7, height = 7)
mr_scatter_plot(mr_results = mrResult, dat = dat_valid)
dev.off()

# 2. 單一 SNP 森林圖
res_single <- mr_singlesnp(dat_valid)
pdf(file.path(file_path, "mr_forest_plot.pdf"), width = 7, height = 7)
mr_forest_plot(res_single)
dev.off()

# 3. 漏斗圖
pdf(file.path(file_path, "mr_funnel_plot.pdf"), width = 7, height = 7)
mr_funnel_plot(res_single)
dev.off()

# 4. 留一法敏感性分析圖
leave1 <- mr_leaveoneout(dat_valid)
pdf(file.path(file_path, "mr_leaveoneout_plot.pdf"), width = 7, height = 7)
mr_leaveoneout_plot(leave1)
dev.off()

#（可選）MR-PRESSO：
# presso <- run_mr_presso(dat)
# write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file = file.path(file_path, "table.MR-PRESSO.csv"))
