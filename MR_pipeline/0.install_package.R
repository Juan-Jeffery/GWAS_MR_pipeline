# 1. 安裝 Bioconductor 套件管理器 
#install.packages("BiocManager", repos="https://cloud.r-project.org")

# 2. 裝 CRAN 上的套件
#install.packages(c(
#  "MendelianRandomization",
#  "CMplot",
#  "R.utils"
#), repos="https://cloud.r-project.org")

# 3. 裝 MendelianRandomization 依賴套件
#install.packages("gmp", repos = "https://cloud.r-project.org")
#install.packages("arrangements", repos = "https://cloud.r-project.org")
#install.packages("iterpc", repos = "https://cloud.r-project.org")
#install.packages("MendelianRandomization", repos="https://cloud.r-project.org")


# 4. 安裝 devtools (可用conda裝)
# install.packages("devtools")

# 5. 用 devtools 裝 TwoSampleMR 依賴套件
#install.packages("nloptr", repos = "https://cloud.r-project.org")
#install.packages("lme4", repos = "https://cloud.r-project.org")
#install.packages("meta", repos = "https://cloud.r-project.org")
#devtools::install_github("MRCIEU/TwoSampleMR")

# 6. 用 devtools 裝 mr.raps 套件
#devtools::install_github("qingyuanzhao/mr.raps")

# 7. 用 devtools 裝 主要執行套件
#Sys.setenv(GITHUB_PAT = "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")
#BiocManager::install("rtracklayer")
#BiocManager::install("VariantAnnotation")
#devtools::install_github("MRCIEU/gwasvcf")
#devtools::install_github("mrcieu/gwasglue")
#BiocManager::install("biomaRt")

