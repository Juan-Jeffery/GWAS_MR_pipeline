#!/usr/bin/env Rscript
library(data.table)

# ==== 1. Configuration ====
extract_dir <- "/home/rd01/data/GTEx/GTEx_Analysis_v8_eQTL/"
target_gene <- "ENSG00000007171"
output_file <- "/home/rd01/test_result/test_all_tissue/ENSG00000007171_snp_counts.tsv"

# Find all eQTL significance files
files <- list.files(extract_dir, pattern = "signif_variant_gene_pairs.txt.gz$", 
                    full.names = TRUE, recursive = TRUE)

message("Found ", length(files), " tissue files. Starting processing...")

# ==== 2. Process Tissues ====
results_list <- lapply(files, function(f) {
  # Extract tissue name from filename
  tissue_name <- tstrsplit(basename(f), "\\.", fixed = TRUE)[[1]]
  message("Processing: ", tissue_name)
  
  # Read and filter using data.table for speed
  dt <- fread(f, sep = "\t", select = c("gene_id", "variant_id"))
  
  # Strip gene version and filter
  dt[, gene_base := sub("\\..*", "", gene_id)]
  subset_dt <- dt[gene_base == target_gene]
  
  # Count unique SNPs
  return(data.table(
    Tissue = tissue_name, 
    SNP_count = uniqueN(subset_dt$variant_id)
  ))
})

# ==== 3. Aggregate & Export ====
final_results <- rbindlist(results_list)

# Create output directory and save
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
fwrite(final_results, output_file, sep = "\t", quote = FALSE)

message("Finished. Results saved to: ", output_file)