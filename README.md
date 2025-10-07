# 🧬 GWAS MR Pipeline

This pipeline performs **Mendelian Randomization (MR)** analysis using GWAS summary statistics.  
Please make sure the environment and parameters are correctly configured before running.

---

## 🚀 Environment Setup

Before execution, activate the Conda environment:

```bash
conda activate r440_env
```

---

## ⚙️ Parameter Settings

Make sure all parameters are properly filled in:

| Parameter | Description |
|------------|-------------|
| **GWAS API** | GWAS Catalog query API |
| **exposure_data_raw** | File path to the exposure dataset |
| **outcome_data_raw** | File path to the outcome dataset |
| **file_path** | Directory to store the output results |
| **pval_threshold** | P-value threshold for selecting significant SNPs |
| **exposure_n** | Sample size of the exposure dataset |
| **pop_n** | Population code (e.g., EUR, EAS, AFR) |
| **exposure_name** | Name of the exposure (used in output) |
| **outcome_name** | Name of the outcome (used in output) |
| **exposure_variant_id** | Variant ID column in the exposure data |
| **outcome_variant_id** | Variant ID column in the outcome data |

---

## 🧩 Execution

Run the main MR pipeline script:

```bash
nohup Rscript code/Run_MR_pipeline.R > /home/rd01/result/hb_t2d/697/Process_Detail.log 2>&1 &
```

Explanation:
- `nohup`: Keeps the process running in the background (even after logout)  
- `Rscript code/Run_MR_pipeline.R`: Path to the main R script  
- `>`: Redirects output to a log file  
- `2>&1`: Sends error messages to the same log file  
- `&`: Runs the process in the background  

---

## 📂 Output

After the analysis completes, results will be saved in the specified `file_path` directory, including:
- MR result tables  
- Process log file (`Process_Detail.log`)  
- Intermediate files (e.g., harmonized data, filtered SNP list)

---

## 📖 Notes

- Ensure all required R packages are installed before running `Run_MR_pipeline.R`  
- Exposure and outcome datasets should use the same genome build (recommended: GRCh37 or GRCh38)  
- When using custom GWAS files, confirm the column names match the pipeline requirements  
