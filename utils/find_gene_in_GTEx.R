library(data.table)
library(stringr)

# 路徑設定
tar_path <- "/home/rd01/data/GTEx/GTEx_Analysis_v8_eQTL.tar"
extract_dir <- "/home/rd01/data/GTEx/GTEx_Analysis_v8_eQTL/"
target_gene <- "ENSG00000007171"
output_file <- "/home/rd01/test_result/test_all_tissue/ENSG00000007171_snp_counts.tsv"

# Step 1: 解壓縮 tar 檔
#untar(tar_path, exdir = extract_dir)

# Step 2: 找出所有組織檔案
files <- list.files(extract_dir, pattern = "signif_variant_gene_pairs.txt.gz$", full.names = TRUE, recursive = TRUE)

results <- data.frame(Tissue = character(), SNP_count = integer(), stringsAsFactors = FALSE)

# Step 3: 迴圈處理每個檔案
for (f in files) {
  tissue <- strsplit(basename(f), "\\.")[[1]][1]  # tissue name
  message("處理中: ", tissue)
  
  # 讀取檔案
  dt <- fread(cmd = paste("zcat", f), sep = "\t", header = TRUE, data.table = FALSE)
  
  # 去掉基因版本號
  dt$gene_id_noversion <- sub("\\..*", "", dt$gene_id)
  
  # 篩選基因
  subset <- dt[dt$gene_id_noversion == target_gene, ]
  
  # 計算 SNP 數量
  snp_count <- length(unique(subset$variant_id))
  
  results <- rbind(results, data.frame(Tissue = tissue, SNP_count = snp_count, stringsAsFactors = FALSE))
}

# Step 4: 輸出結果
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.table(results, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

message("finish ", output_file)

# Step 5: 移除解壓縮資料夾
#unlink(extract_dir, recursive = TRUE)
message("remove dir ", extract_dir)
