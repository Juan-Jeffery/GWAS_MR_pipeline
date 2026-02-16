# ─────────────── utils ───────────────
Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJqamo4NTM3MTRAZ21haWwuY29tIiwiaWF0IjoxNzcwMjY5ODY2LCJleHAiOjE3NzE0Nzk0NjZ9.WlmLSEY2hjqGXOI8kThLlcv3btbwphHXqDfi9Uh5jdYEzsGf9WtkZH4WhgpOjoCmB1NAv7Q084s8lnmPD2RKTDR6SQyFdFohvMeomWVDQGpGZq3gzUfPtud-sQbKPSuO5uS_vWGcYTJ_ZXA2w79tjaGj3auTrPVt-0YQ8FbDrTX2-wzT5-vu2xvN_ijRdezxOBb4muE6Gf7Ym57fjg-Dq8stXdGdaqI1216nVHg1mwrfdpw27E050eOg_7V7g60n9Yz979xgn00pY5lV1EhZbuLwO4wGmgwT9LSMNDUB8c_Kt3f90H4FlLDgwHaA44mv_K4HOdjH330-pDFs2GNMNA")
options(ieugwasr_api = "https://gwas-api.mrcieu.ac.uk/")


source("/home/rd01/code/MR_pipeline/1.loading_package.R")
source("/home/rd01/code/MR_pipeline/2.loading_file.R")
source("/home/rd01/code/MR_pipeline/3.exposure_data_preprocessing.R")
source("/home/rd01/code/MR_pipeline/4.MR.R")
source("/home/rd01/code/MR_pipeline/log.R")

# ─────────────── Parameters ───────────────
exposure_data_raw <- "/home/rd01/data/OSA/GCST90475825.tsv.gz"
outcome_data_raw <- "/home/rd01/data/t2d_GCST90296697_gwascatalog/GCST90296697.tsv.gz"
file_path <- "/home/rd01/result/osa_t2d/825_697"
setwd(file_path) 

pval_threshold <- 5e-6
exposure_n <- 430325
pop_n <- "EUR"
exposure_name <- "OSA"
outcome_name <- "T2D"
exposure_variant_id <- "rsid"
#exposure_variant_id <- "variant_id"
outcome_variant_id <- "rs_id" 
#704,934=variant_id, other=rs_id

# ─────────────── GWAS MR pipeline ───────────────
init_log()
# === Start process ===
start_time_all <- Sys.time()
log_msg("START", msg = "All process started")

tryCatch({
  load_packages()
  log_msg("SUCCESS", "load_packages", "All packages loaded successfully.")
}, error = function(e) {
  log_msg("ERROR", "load_packages", e$message)
  stop(e$message)
})

tryCatch({
  cat("=== Checking OpenGWAS API connection ===\n")
  cat("JWT:\n", ieugwasr::get_opengwas_jwt(), "\n")
  gwas_info <- ieugwasr::gwasinfo("ieu-a-2")
  user_info <- ieugwasr::user()
  print(gwas_info)
  print(user_info)
  log_msg("SUCCESS", "check_api_connection", "API connection successful.")
}, error = function(e) {
  log_msg("ERROR", "check_api_connection", e$message)
  stop("API connection failed. Please check if your JWT is valid and reset it using Sys.setenv(OPENGWAS_JWT = 'your_token')")
})


tryCatch({
  exposure_data <- read_and_format_exposure_tsv(exposure_data_raw, exposure_variant_id)
  log_msg("SUCCESS", "read_and_format_exposure_tsv", "Exposure data loaded")
}, error = function(e) {
  log_msg("ERROR", "read_and_format_exposure_tsv", e$message)
  stop(e$message)
})

tryCatch({
  outcome_data <- read_and_format_outcome_tsv(outcome_data_raw)
  log_msg("SUCCESS", "read_and_format_outcome_tsv", "Outcome data loaded")
}, error = function(e) {
  log_msg("ERROR", "read_and_format_outcome_tsv", e$message)
  stop(e$message)
})

tryCatch({
  filtered_exposure_data <- plot_cmplot(file_path, exposure_data, pval_threshold)
  log_msg("SUCCESS", "plot_cmplot", "Exposure filtered with p-value threshold")
}, error = function(e) {
  log_msg("ERROR", "plot_cmplot", e$message)
  stop(e$message)
})

tryCatch({
  clumped_exposure_data <- ld_clump_exposure(file_path, filtered_exposure_data, pop_n)
  log_msg("SUCCESS", "ld_clump_exposure", "LD clumping finished")
}, error = function(e) {
  log_msg("ERROR", "ld_clump_exposure", e$message)
  stop(e$message)
})

tryCatch({
  final_exposure_data <- filter_by_F_statistic(file_path, clumped_exposure_data, exposure_n)
  log_msg("SUCCESS", "filter_by_F_statistic", "Exposure SNPs filtered by F-statistic")
}, error = function(e) {
  log_msg("ERROR", "filter_by_F_statistic", e$message)
  stop(e$message)
})

tryCatch({
  dat_valid <- process_outcome_data(file_path, final_exposure_data, outcome_data, outcome_variant_id, exposure_name, outcome_name)
  log_msg("SUCCESS", "process_outcome_data", "Exposure and outcome harmonized")
}, error = function(e) {
  log_msg("ERROR", "process_outcome_data", e$message)
  stop(e$message)
})

tryCatch({
  mr_result <- run_mr_analysis(file_path, dat_valid)
  log_msg("SUCCESS", "run_mr_analysis", "MR analysis finished")
}, error = function(e) {
  log_msg("ERROR", "run_mr_analysis", "See process_detail.log for details")
  stop(e$message)
})

tryCatch({
  plot_mr_results(file_path, mr_result, dat_valid)
  log_msg("SUCCESS", "plot_mr_results", "MR plots saved")
  end_time_all <- Sys.time()
  log_msg("END", msg = paste("All process completed. Duration:", round(difftime(end_time_all, start_time_all, units = "mins"), 2), "mins"))
}, error = function(e) {
  log_msg("ERROR", "plot_mr_results", e$message)
  stop(e$message)
})

