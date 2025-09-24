#!/usr/bin/env Rscript
library(data.table)

start_time <- Sys.time()
cat("=== Script started at:", start_time, "===\n")

# ==== 參數設定 ====
tissue_file <- "/home/rd01/data/GTEx/GTEx_Analysis_v8_eQTL/Muscle_Skeletal.v8.signif_variant_gene_pairs.txt.gz"
vcf_file    <- "/home/rd01/biomart_ref/All_20180418.vcf.gz"
snp_out     <- "/home/rd01/test_result/test_catch_snp2/Muscle_Skeletal_TargetGenes_eQTL_result.tsv"
mr_out      <- "/home/rd01/test_result/test_catch_snp2/Muscle_Skeletal_TargetGenes_eQTL_MR_result.tsv"

# samplesize
samplesize  <- 803 #lung(515), Muscle(803)

# 定義目標基因 & 方向校正 
gene_directions <- data.table(
  gene_id = c("ENSG00000116016","ENSG00000135766","ENSG00000110628",
              "ENSG00000186951","ENSG00000162409","ENSG00000007171"),
  direction = c("down","down","up","down","up","up")
)
target_genes <- gene_directions$gene_id

# ==== Run ====
# 讀取 eQTL
cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Loading tissue eQTL data...\n")
eqtl <- fread(tissue_file)
cat("Total eQTL rows:", nrow(eqtl), "\n\n")
cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Filtering target genes...\n")
target_eqtl <- eqtl[gene_id %like% paste0("^(", paste(target_genes, collapse="|"), ")")]
cat("Target eQTL rows:", nrow(target_eqtl), "\n")
cat("Gene counts:\n")
print(table(sub("\\..*$","",target_eqtl$gene_id)))
cat("\n")

# 去掉版本號方便對應方向
target_eqtl[, gene_id_short := sub("\\..*$", "", gene_id)]
target_eqtl[, direction := fifelse(
  gene_id_short %in% gene_directions[direction == "down", gene_id], "down", "up"
)]
target_eqtl[, beta_corrected := ifelse(direction == "down", -1 * slope, slope)]
target_eqtl[, se_corrected := slope_se]

# ==== 解析 variant_id 拆 effect_allele / other_allele ====
#target_eqtl[, c("chr","pos","other_allele","effect_allele","b38") := tstrsplit(variant_id, "_")]
#target_eqtl[, pos := as.integer(pos)]
snp_split <- tstrsplit(target_eqtl$variant_id, "_", fixed=TRUE)
target_eqtl[, chr := gsub("^chr","", snp_split[[1]])]
target_eqtl[, pos := as.integer(snp_split[[2]])]
target_eqtl[, ref_from_id := snp_split[[3]]]
target_eqtl[, alt_from_id := snp_split[[4]]]
target_eqtl[, effect_allele := alt_from_id]
target_eqtl[, other_allele := ref_from_id]

fwrite(target_eqtl, snp_out, sep="\t", na="NA", quote=FALSE)


# ==== 計算 eaf（直接用 tissue maf） ====
target_eqtl[, eaf := maf]
cat("eaf summary:\n")
print(summary(target_eqtl$eaf))
cat("\n")

# ==== 載入 dbSNP VCF 並匹配 rsID ====
cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Loading dbSNP VCF...\n")
vcf <- fread(vcf_file, skip="#CHROM", select=1:5,
             col.names=c("chr","pos","rsid","ref","alt"))
vcf[, chr := gsub("^chr","", chr)]

cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Matching rsID (forward and reverse)...\n")

# 若沒有 effect_allele / other_allele，就用 SNP 名稱裡的 REF/ALT 先帶
if (!"effect_allele" %in% names(target_eqtl)) target_eqtl[, effect_allele := alt_from_id]
if (!"other_allele"  %in% names(target_eqtl)) target_eqtl[, other_allele  := ref_from_id]

m1 <- merge(target_eqtl, vcf, by.x=c("chr","pos","effect_allele","other_allele"),
            by.y=c("chr","pos","alt","ref"), all.x=TRUE, suffixes=c("", "_vcf"))
m2 <- merge(target_eqtl, vcf, by.x=c("chr","pos","other_allele","effect_allele"),
            by.y=c("chr","pos","alt","ref"), all.x=TRUE, suffixes=c("", "_vcf"))
#m1$rsid <- ifelse(is.na(m1$rsid), m2$rsid, m1$rsid)
m1$rsid <- ifelse(is.na(m1$rsid), m2$rsid[match(m1$variant_id, m2$variant_id)], m1$rsid)

cat("Matched rsID count:", sum(!is.na(m1$rsid)), "/", nrow(target_eqtl), "\n\n")

# ==== MR-ready 格式輸出 & deduplicate by rsID (min p-value) ====
cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Creating MR-ready dataset...\n")
mr_dt <- m1[!is.na(rsid), .(
  chromosome = chr,
  variant_id = rsid,
  base_pair_location = pos,
  effect_allele,
  other_allele,
  N = samplesize,
  effect_allele_frequency = eaf,
  beta = beta_corrected,
  standard_error = se_corrected,
  p_value = pval_nominal,
  gene_id
)]

# 去除重複 rsID，保留最小 p-value
mr_dt <- mr_dt[order(p_value)][, .SD[1], by=variant_id]
cat("After deduplicate (min p-value) rows:", nrow(mr_dt), "\n")
cat("Sample of MR-ready data:\n")
cat("Gene snp with rsID counts:\n")
print(table(sub("\\..*$","",mr_dt$gene_id)))
cat("\n")
print(head(mr_dt, 5))

# ==== 輸出檔案 ====
fwrite(mr_dt, mr_out, sep="\t", na="NA", quote=FALSE)

# ==== 統計 & 結束 ====
end_time <- Sys.time()
cat("\n=== Script finished at:", end_time, "===\n")
cat("Duration:", round(difftime(end_time, start_time, units="secs"), 2), "seconds\n")
