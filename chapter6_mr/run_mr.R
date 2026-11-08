# Run the numbered scripts in order. Each step can also be sourced separately.
m <- local({
  paths <- Filter(is.character, lapply(sys.frames(), function(frame) frame$ofile))
  script <- if (length(paths)) tail(paths, 1L)[[1L]] else
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(script) == 1L) dirname(normalizePath(script, mustWork = TRUE)) else getwd()
})
packages <- c("data.table", "R.utils", "yaml", "jsonlite", "TwoSampleMR")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
steps <- c("01_define_pairs.R", "02_read_gwas.R", "03_match_instruments.R",
  "04_validate_proxies.R", "05_run_mr_compare.R")
for (step in steps) source(file.path(m, step), encoding = "UTF-8")
if (interactive()) View(published_comparison)
