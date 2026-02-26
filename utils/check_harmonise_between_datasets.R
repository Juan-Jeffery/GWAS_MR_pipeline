library(data.table)

# Standard T2D SNPs to check
check_snps <- c("rs7903146", "rs2237892", "rs13266634", "rs11549407")

datasets <- list(
  GCST90296697 = "data/t2d_GCST90296697_gwascatalog/GCST90296697.tsv.gz",
  GCST90435704 = "data/t2d_GCST90435704_gwascatalog/GCST90435704.tsv.gz",
  GCST90006934 = "data/t2d_GCST90006934_gwascatalog/GCST90006934_buildGRCh37.tsv.gz"
)

results_list <- list()

for (name in names(datasets)) {
  message("Checking dataset: ", name)
  
  # Read header once to find column names
  header <- colnames(fread(datasets[[name]], nrows = 0))
  beta_col <- grep("^beta$|^odds_ratio$", header, ignore.case = TRUE, value = TRUE)[1]
  ea_col <- grep("effect_allele", header, ignore.case = TRUE, value = TRUE)[1]
  oa_col <- grep("other_allele", header, ignore.case = TRUE, value = TRUE)[1]
  
  for (snp in check_snps) {
    # Search one by one using zgrep for safety
    cmd <- paste0("zgrep -w '", snp, "' ", datasets[[name]])
    row_text <- try(system(cmd, intern = TRUE), silent = TRUE)
    
    if (length(row_text) > 0 && !inherits(row_text, "try-error")) {
      # Use the first match if multiple rows return
      dt_row <- fread(text = row_text[1], header = FALSE)
      if (ncol(dt_row) == length(header)) {
        setnames(dt_row, header)
        
        val <- as.numeric(dt_row[[beta_col]])
        type_label <- ifelse(grepl("beta", beta_col, ignore.case = TRUE), "Beta", "OR")
        
        # Determine Direction
        direction <- "Unknown"
        if (type_label == "Beta") {
          direction <- ifelse(val > 0, "Increase", "Decrease")
        } else {
          direction <- ifelse(val > 1, "Increase", "Decrease")
        }
        
        results_list[[paste(name, snp)]] <- data.table(
          Study = name,
          SNP = snp,
          EA = dt_row[[ea_col]],
          OA = dt_row[[oa_col]],
          Value = val,
          Type = type_label,
          Direction = direction
        )
      }
    }
  }
}

final_table <- rbindlist(results_list)
setorder(final_table, SNP, Study)
print(final_table)