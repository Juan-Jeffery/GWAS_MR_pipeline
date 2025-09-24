# ─────────────── log.R ───────────────
# 初始化並清空 log 檔案，並寫入開頭標語
init_log <- function(log_file = "Process.log") {
  if (file.exists(log_file)) {
    file.remove(log_file)
  }
  file.create(log_file)
  
  banner <- paste0(
    "=== Log started at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ===\n"
  )
  # 只寫檔案，不輸出到 console
  write(banner, file = log_file, append = TRUE)
}

# 記錄 log 訊息（附時間戳）
log_msg <- function(type = "INFO", func = NULL, msg, log_file = "Process.log") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s%s: %s\n", timestamp, type,
                  if (!is.null(func)) paste0(" in ", func) else "", msg)

  cat(line)  # 印到 console
  write(line, file = log_file, append = TRUE)  # 寫到 log 檔（保留舊紀錄）
}
