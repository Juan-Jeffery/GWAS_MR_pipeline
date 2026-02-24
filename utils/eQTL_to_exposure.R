#!/usr/bin/env Rscript
library(data.table)
library(VariantAnnotation)
library(GenomicRanges)

# ==== 1. Configuration ====
# Input paths
tissue_file <- "/home/rd01/data/GTEx/GTEx_Analysis_v8_eQTL/GTEx_Analysis_v8_eQTL/Muscle_Skeletal.v8.signif_variant_gene_pairs.txt.gz"
vcf_file    <- "/home/rd01/biomart_ref/All_20180418.vcf.gz"

# Output paths
snp_out     <- "/home/rd01/result/GTEx_snp/PRKAA1/PRKAA1_Muscle_Skeletal_TargetGenes_eQTL_result.tsv"
mr_out      <- "/home/rd01/result/GTEx_snp/PRKAA1/PRKAA1_Muscle_Skeletal_TargetGenes_eQTL_MR_result.tsv"

# Metadata
sample_size <- 803
# lung(515), Muscle(803), Whole_Blood(755), Brain_Hippocampus(150), Brain_Nucleus_accumbens_basal_ganglia(153), Nerve_Tibial(670)
# Brain_Cortex(184), Heart_Atrial_Appendage(322), Artery_Tibial(489), Artery_Aorta(432), Artery_Coronary(213), Heart_Left_Ventricle(334), Adipose_Subcutaneous(581), Brain_Caudate(172)

gene_targets <- data.table(
  gene_id   = c("ENSG00000162409"),
  direction = c("up") # Options: 'up' or 'down' for beta correction
)

start_time <- Sys.time()
message("=== Script started at: ", start_time, " ===")

# ==== 2. Load & Filter eQTL Data ====
message("[", Sys.time(), "] Loading eQTL data...")
eqtl <- fread(tissue_file)

# Filter for target genes
pattern <- paste0("^(", paste(gene_targets$gene_id, collapse="|"), ")")
target_eqtl <- eqtl[gene_id %like% pattern]

# Directional beta correction
target_eqtl[, gene_id_short := sub("\\..*$", "", gene_id)]
target_eqtl[, beta_corrected := ifelse(
  gene_id_short %in% gene_targets[direction == "down", gene_id], -slope, slope
)]

# Parse variant_id (format: chr_pos_ref_alt_build)
snp_info <- tstrsplit(target_eqtl$variant_id, "_", fixed=TRUE)
target_eqtl[, `:=`(
  chr           = gsub("^chr", "", snp_info[[1]]),
  pos           = as.integer(snp_info[[2]]),
  effect_allele = snp_info[[4]],
  other_allele  = snp_info[[3]],
  eaf           = maf
)]

# ==== 3. Map rsIDs using dbSNP VCF ====
message("[", Sys.time(), "] Querying dbSNP VCF...")
gr_targets <- GRanges(
  seqnames = target_eqtl$chr,
  ranges   = IRanges(start = target_eqtl$pos, end = target_eqtl$pos)
)

# Subset VCF to specific coordinates
vcf_param <- ScanVcfParam(which = gr_targets)
vcf_data  <- readVcf(vcf_file, "hg38", param = vcf_param)

# Extract reference info from VCF
vcf_ref <- data.table(
  chr  = as.character(seqnames(rowRanges(vcf_data))),
  pos  = start(rowRanges(vcf_data)),
  rsid = names(rowRanges(vcf_data)),
  ref  = as.character(ref(vcf_data)),
  alt  = as.character(unlist(alt(vcf_data)))
)

# Bidirectional merging to handle strand/allele flips
m1 <- merge(target_eqtl, vcf_ref, 
            by.x = c("chr", "pos", "effect_allele", "other_allele"),
            by.y = c("chr", "pos", "alt", "ref"), all.x = TRUE)

m2 <- merge(target_eqtl, vcf_ref, 
            by.x = c("chr", "pos", "other_allele", "effect_allele"),
            by.y = c("chr", "pos", "alt", "ref"), all.x = TRUE)

# Prioritize m1 matches, fill gaps with m2
m1$rsid <- ifelse(is.na(m1$rsid), m2$rsid[match(m1$variant_id, m2$variant_id)], m1$rsid)
message("Matched rsIDs: ", sum(!is.na(m1$rsid)), " / ", nrow(target_eqtl))

# ==== 4. Format & Export MR Data ====
message("[", Sys.time(), "] Formatting final MR dataset...")

mr_final <- m1[!is.na(rsid), .(
  chromosome              = chr,
  variant_id              = rsid,
  base_pair_location      = pos,
  effect_allele,
  other_allele,
  N                       = sample_size,
  effect_allele_frequency = eaf,
  beta                    = beta_corrected,
  standard_error          = slope_se,
  p_value                 = pval_nominal,
  gene_id
)]

# Deduplicate by variant_id (keep lowest p-value)
mr_final <- mr_final[order(p_value)][, .SD[1], by = variant_id]

# Export results
fwrite(target_eqtl, snp_out, sep="\t", na="NA", quote=FALSE)
fwrite(mr_final, mr_out, sep="\t", na="NA", quote=FALSE)

message("=== Script finished in: ", round(difftime(Sys.time(), start_time, units="secs"), 2), "s ===")