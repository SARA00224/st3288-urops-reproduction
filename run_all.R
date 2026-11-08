# Run both experiments in separate R sessions using the normal R library search path.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L || (length(args) == 1L && args != "--full")) {
  stop("Usage: Rscript --vanilla run_all.R [--full]")
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(script) != 1L) stop("Run this entry point with Rscript.")
root <- dirname(normalizePath(script, mustWork = TRUE))

packages <- c("data.table", "R.utils", "yaml", "jsonlite", "TwoSampleMR")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install the packages listed in README.md. Missing: ",
  paste(missing, collapse = ", "))

Sys.setenv(MR_INPUT_MODE = if (length(args)) "full" else "excerpt")
rscript <- file.path(R.home("bin"), "Rscript")
entries <- c("section5_1_simulation/run_all.R", "section5_2_mr/run_mr.R")
for (entry in entries) {
  message("Running ", entry)
  status <- system2(rscript, c("--vanilla", shQuote(file.path(root, entry))))
  if (status != 0L) stop(entry, " exited with status ", status)
}
message("Both experiments completed. See their results directories.")
