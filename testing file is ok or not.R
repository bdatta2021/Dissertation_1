

files <- list.files("Visium/processed_rds", pattern="\\.rds$", full.names=TRUE)
names <- basename(files)
for (i in seq_along(files)) {
  cat(sprintf("[%02d] Testing %s … ", i, names[i]))
  res <- tryCatch({
    readRDS(files[i])
    "OK\n"
  }, error = function(e) {
    paste0("FAILED: ", e$message, "\n")
  })
  cat(res)
}
