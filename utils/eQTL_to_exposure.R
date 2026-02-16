#!/usr/bin/env Rscript
library(data.table)
library(VariantAnnotation)
library(GenomicRanges)

start_time <- Sys.time()
cat("=== Script started at:", start_time, "===\n") 

# ==== 參數設定 ====
tissue_file <- "/home/rd01/data/GTEx/GTEx_Analysis_v8_eQTL/GTEx_Analysis_v8_eQTL/Muscle_Skeletal.v8.signif_variant_gene_pairs.txt.gz"
vcf_file    <- "/home/rd01/biomart_ref/All_20180418.vcf.gz"
snp_out     <- "/home/rd01/result/GTEx_snp/PRKAA1/PRKAA1_Muscle_Skeletal_TargetGenes_eQTL_result.tsv"
mr_out      <- "/home/rd01/result/GTEx_snp/PRKAA1/PRKAA1_Muscle_Skeletal_TargetGenes_eQTL_MR_result.tsv"

# sample size
samplesize  <- 803
# lung(515), Muscle(803), Whole_Blood(755), Brain_Hippocampus(150), Brain_Nucleus_accumbens_basal_ganglia(153), Nerve_Tibial(670)
# Brain_Cortex(184), Heart_Atrial_Appendage(322), Artery_Tibial(489), Artery_Aorta(432), Artery_Coronary(213), Heart_Left_Ventricle(334), Adipose_Subcutaneous(581), Brain_Caudate(172)

# 定義目標基因 & 方向校正 
gene_directions <- data.table(
  gene_id = c("ENSG00000162409"),
  direction = c("up")
)
target_genes <- gene_directions$gene_id

# ==== Run ====
# 讀取 eQTL
cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Loading tissue eQTL data...\n")
eqtl <- fread(tissue_file)
cat("Total eQTL rows:", nrow(eqtl), "\n\n")

# 過濾目標基因
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

# ==== 用 VariantAnnotation 只抓目標 SNP ==== 
cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Loading dbSNP VCF for target SNP only...\n")
gr <- GRanges(
  seqnames = target_eqtl$chr,
  ranges = IRanges(start = target_eqtl$pos, end = target_eqtl$pos)
)
param <- ScanVcfParam(which = gr, info = NA, geno = NA)
vcf <- readVcf(vcf_file, "hg38", param = param)

vcf_df <- data.table(
  chr  = as.character(seqnames(rowRanges(vcf))),
  pos  = start(rowRanges(vcf)),
  rsid = names(rowRanges(vcf)),
  ref  = as.character(ref(vcf)),
  alt  = as.character(unlist(alt(vcf)))
)

cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Matching rsID (forward and reverse)...\n")
if (!"effect_allele" %in% names(target_eqtl)) target_eqtl[, effect_allele := alt_from_id]
if (!"other_allele"  %in% names(target_eqtl)) target_eqtl[, other_allele  := ref_from_id]

m1 <- merge(target_eqtl, vcf_df,
            by.x=c("chr","pos","effect_allele","other_allele"),
            by.y=c("chr","pos","alt","ref"), all.x=TRUE)
m2 <- merge(target_eqtl, vcf_df,
            by.x=c("chr","pos","other_allele","effect_allele"),
            by.y=c("chr","pos","alt","ref"), all.x=TRUE)
m1$rsid <- ifelse(is.na(m1$rsid),
                  m2$rsid[match(m1$variant_id, m2$variant_id)],
                  m1$rsid)

cat("Matched rsID count:", sum(!is.na(m1$rsid)), "/", nrow(target_eqtl), "\n\n")

# ==== MR-ready 格式輸出 & deduplicate by rsID ====
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
mr_dt <- mr_dt[order(p_value)][, .SD[1], by=variant_id]

cat("After deduplicate (min p-value) rows:", nrow(mr_dt), "\n")
cat("Gene snp with rsID counts:\n")
print(table(sub("\\..*$","",mr_dt$gene_id)))
cat("\n")
print(head(mr_dt,5))

# ==== 輸出檔案 ====
fwrite(mr_dt, mr_out, sep="\t", na="NA", quote=FALSE)

# ==== 統計 & 結束 ====
end_time <- Sys.time()
cat("\n=== Script finished at:", end_time, "===\n")
cat("Duration:", round(difftime(end_time, start_time, units="secs"), 2), "seconds\n")
