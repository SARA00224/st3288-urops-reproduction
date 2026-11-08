# Run with Rscript or source() using the path to this file.
local({
  # source() supplies ofile, while Rscript supplies a --file argument.
  source_files <- lapply(rev(sys.frames()), function(frame) frame$ofile)
  source_files <- Filter(function(path) is.character(path) &&
    length(path) == 1L && basename(path) == "run_all.R", source_files)
  script <- if (length(source_files)) source_files[[1L]] else
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  module <- if (length(script) == 1L && file.exists(script))
    dirname(normalizePath(script, winslash = "/")) else getwd()
  if (!all(file.exists(file.path(module, c("config.R", "R/functions.R")))))
    stop("Cannot locate the Section 5.1 scripts. Run Rscript with the path to run_all.R.")
  source(file.path(module, "config.R"), local = TRUE)
  source(file.path(module, "R", "functions.R"), local = TRUE)
  stopifnot(config$n_snps >= 2, length(config$cohort_sizes) == 2,
    all(config$cohort_sizes > 10), config$noise_sd > 0,
    config$int_offset > 0, config$int_offset < 1)
  for (folder in c("data/simulated", "data/processed", "results/tables", "results/figures"))
    dir.create(file.path(module, folder), recursive = TRUE, showWarnings = FALSE)
  set.seed(config$seed, kind = "Mersenne-Twister",
    normal.kind = "Inversion", sample.kind = "Rejection")
  maf <- runif(config$n_snps, .1, .45)
  raw <- lapply(names(config$cohort_sizes), function(k)
    simulate_cohort(config$cohort_sizes[[k]], k, config, maf))
  processed <- lapply(raw, preprocess_cohort, offset = config$int_offset)
  raw_data <- do.call(rbind, raw)
  processed_data <- do.call(rbind, processed)
  covariate_coefficients <- do.call(rbind, lapply(processed, attr,
    which = "covariate_coefficients"))
  simulation_variants <- data.frame(
    SNP=sprintf("sim_snp_%03d", seq_len(config$n_snps)),
    effect_allele="A", other_allele="G", generating_frequency=maf,
    true_beta_log=c(config$effect_log_scale, rep(0,config$n_snps-1)))
  variants <- sprintf("sim_snp_%03d", seq_len(config$n_snps))
  per_cohort <- do.call(rbind, lapply(processed, associate_cohort,
    variants = variants, protein = config$protein))
  meta <- meta_analyse(per_cohort)
  # Each cohort's contribution to the same SNP's meta estimate.
  meta_contributions <- per_cohort[,c("protein","cohort","SNP",
    "effect_allele","other_allele","beta","se","p","n")]
  meta_contributions$weight <- 1/meta_contributions$se^2
  meta_contributions$sum_weights <- ave(meta_contributions$weight,
    meta_contributions$SNP, FUN=sum)
  meta_contributions$weight_fraction <- meta_contributions$weight /
    meta_contributions$sum_weights
  meta_contributions$weighted_beta <- meta_contributions$weight * meta_contributions$beta
  meta_contributions$beta_contribution <- meta_contributions$weight_fraction *
    meta_contributions$beta
  write.csv(raw_data, file.path(module, "data/simulated/individuals.csv"), row.names=FALSE)
  write.csv(processed_data, file.path(module, "data/processed/individuals_processed.csv"), row.names=FALSE)
  write.csv(per_cohort, file.path(module, "results/tables/cohort_associations.csv"), row.names=FALSE)
  write.csv(meta, file.path(module, "results/tables/meta_associations.csv"), row.names=FALSE)
  write.csv(simulation_variants, file.path(module,"results/tables/simulation_variants.csv"), row.names=FALSE)
  write.csv(covariate_coefficients, file.path(module,"results/tables/covariate_coefficients.csv"), row.names=FALSE)
  write.csv(meta_contributions, file.path(module,"results/tables/meta_contributions.csv"), row.names=FALSE)
  # Adjust raw and log outcomes in the regression. The INT outcome is preprocessed.
  comparison <- do.call(rbind, lapply(processed, function(d) {
    do.call(rbind, lapply(c("protein_raw", "protein_log", "protein_INT"), function(y) {
      f <- if (y == "protein_INT") paste(y, "~ sim_snp_001") else
        paste(y, "~ sim_snp_001 + age + sex + processing_time + PC1 + PC2 + PC3")
      s <- coef(summary(lm(as.formula(f), data=d)))["sim_snp_001", ]
      data.frame(cohort=d$cohort[1], outcome=y, beta=s[1], se=s[2], p=s[4])
    }))
  }))
  write.csv(comparison, file.path(module, "results/tables/scale_comparison.csv"), row.names=FALSE)
  png(file.path(module, "results/figures/preprocessing.png"), width=1500, height=800)
  tryCatch({
    par(mfrow=c(2,4))
    for (d in processed) for (y in c("protein_raw", "protein_log", "protein_residual", "protein_INT"))
      hist(d[[y]], breaks=35, main=paste(d$cohort[1], y), xlab=y, col="steelblue")
  }, finally=dev.off())
  capture.output(sessionInfo(), file=file.path(module, "results/sessionInfo.txt"))
  dput(config, file=file.path(module, "results/config_used.R"))
  stopifnot(nrow(meta)==config$n_snps, all(is.finite(meta$beta)),
    all(meta$p >= 0 & meta$p <= 1), all(meta$n == sum(config$cohort_sizes)))
  cat("Completed: simulated data -> cohort preprocessing -> SNP regressions -> fixed-effect meta-analysis.\n")
  output_index <- data.frame(
    section=c("3","3.1","4","4.2","5","6","6","7"),
    file=c("data/simulated/individuals.csv","results/tables/simulation_variants.csv",
      "data/processed/individuals_processed.csv","results/tables/covariate_coefficients.csv",
      "results/tables/cohort_associations.csv","results/tables/meta_contributions.csv",
      "results/tables/meta_associations.csv","results/tables/scale_comparison.csv"),
    rows=c(nrow(raw_data),nrow(simulation_variants),nrow(processed_data),
      nrow(covariate_coefficients),nrow(per_cohort),nrow(meta_contributions),nrow(meta),nrow(comparison)))
  print(output_index,row.names=FALSE)
  # Expose only the named learning tables; keep simulation internals local.
  learning_tables <- list(raw_data=raw_data, simulation_variants=simulation_variants,
    processed=processed_data, covariate_coefficients=covariate_coefficients,
    cohort_results=per_cohort, meta_contributions=meta_contributions,
    meta_results=meta, comparison=comparison)
  list2env(learning_tables, envir=.GlobalEnv)
  if (interactive()) {
    View(raw_data, title="Section 3 - simulated participants")
    View(processed_data, title="Section 4 - protein preprocessing")
    View(per_cohort, title="Section 5 - SNP regressions")
    View(meta, title="Section 6 - meta analysis")
  }
  message("Saved eight tables and loaded them into the R environment. See README.md.")
})
