# Run the numbered scripts in order. Objects remain visible in Environment.
m <- local({
  paths <- Filter(is.character, lapply(sys.frames(), function(frame) frame$ofile))
  script <- if (length(paths)) tail(paths, 1L)[[1L]] else
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(script) == 1L) dirname(normalizePath(script, mustWork = TRUE)) else getwd()
})
dir.create(file.path(m, "data"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(m, "results/steps"), recursive=TRUE, showWarnings=FALSE)

stopifnot(requireNamespace("TwoSampleMR",quietly=TRUE))
if(!file.exists(file.path(m,"data/04_validation.rds"))) stop("Run 04_validate_proxies.R first.")
mr_input <- readRDS(file.path(m,"data/04_validation.rds"))$mr_input
# 1. Manual Wald ratio and the same first-order SE used by TwoSampleMR.
mr_results <- mr_input[,c("protein","target_snp","outcome_proxy_snp","outcome_id")]
mr_results$beta_manual <- mr_input$beta_Y/mr_input$beta_X
mr_results$se_manual <- mr_input$se_Y/abs(mr_input$beta_X)
mr_results$z <- mr_results$beta_manual/mr_results$se_manual
mr_results$log_p_manual <- log(2)+pnorm(-abs(mr_results$z),log.p=TRUE)
mr_results$p_manual <- exp(mr_results$log_p_manual)
mr_results$OR <- exp(mr_results$beta_manual)
mr_results$OR_lower <- exp(mr_results$beta_manual-qnorm(.975)*mr_results$se_manual)
mr_results$OR_upper <- exp(mr_results$beta_manual+qnorm(.975)*mr_results$se_manual)

# 2. Actual package call. Keep target and proxy IDs distinct; no fake same-rs merge.
mr_results$beta_package <- NA_real_;mr_results$se_package <- NA_real_
mr_results$p_package <- NA_real_
for(i in seq_len(nrow(mr_input))) {
  fit <- TwoSampleMR::mr_wald_ratio(
    b_exp=mr_input$beta_X[i],b_out=mr_input$beta_Y[i],
    se_exp=mr_input$se_X[i],se_out=mr_input$se_Y[i])
  mr_results$beta_package[i] <- fit$b
  mr_results$se_package[i] <- fit$se
  mr_results$p_package[i] <- fit$pval
}
stopifnot(max(abs(mr_results$beta_manual-mr_results$beta_package))<1e-12,
          max(abs(mr_results$se_manual-mr_results$se_package))<1e-12,
          max(abs(mr_results$log_p_manual-log(mr_results$p_package)))<1e-8)

# 3. Read published MR answers only AFTER the computation, for comparison.
paper_results <- read.csv(file.path(m,"sources/paper_mr_results.csv"),stringsAsFactors=FALSE)
paper_results <- paper_results[match(mr_results$protein,paper_results$Exposure),]
comparison <- data.frame(protein=mr_results$protein,
  beta_reproduced=mr_results$beta_package,beta_paper=paper_results$Beta,
  beta_difference=mr_results$beta_package-paper_results$Beta,
  se_reproduced=mr_results$se_package,se_paper=paper_results$SE,
  se_difference=mr_results$se_package-paper_results$SE,
  p_reproduced=mr_results$p_package,p_paper=paper_results$P_value)
comparison$beta_within_1e_6 <- abs(comparison$beta_difference)<1e-6
comparison$se_within_1e_6 <- abs(comparison$se_difference)<1e-6

# Compare with the fresh extraction from the publisher's original workbook.
published <- read.csv(file.path(m, "sources/published_workbook_records.csv"))
key <- match(mr_results$target_snp, published$original_snp)
stopifnot(!anyNA(key), !anyDuplicated(key))
published <- published[key, ]
stopifnot(identical(published$protein, mr_results$protein),
          all(published$published_method == "WR"),
          all(published$original_outcome_id == 12),
          all(published$published_exposure_ea == mr_input$target_ea),
          all(published$published_exposure_beta == mr_input$beta_X),
          all(published$published_exposure_se == mr_input$se_X))
published_comparison <- data.frame(
  protein = mr_results$protein, phenotype = mr_input$outcome,
  original_snp = mr_results$target_snp, outcome_snp = mr_results$outcome_proxy_snp,
  mr_method = "Wald ratio", published_beta = published$published_beta,
  reproduced_beta = mr_results$beta_package,
  absolute_beta_difference = abs(published$published_beta - mr_results$beta_package),
  published_se = published$published_se, reproduced_se = mr_results$se_package,
  absolute_se_difference = abs(published$published_se - mr_results$se_package),
  published_p = published$published_p, reproduced_p = mr_results$p_package,
  absolute_p_difference = abs(published$published_p - mr_results$p_package),
  source_table = published$source_table, source_sheet = published$source_sheet,
  source_row = published$source_row, source_cells = published$source_cells,
  verification_status = "partially verified",
  verification_note = "Published values verified directly; historical outcome proxy unresolved")

# 4. Separate sensitivity: the SAME proxy's own protein and disease betas.
# This changes the instrument; it is not the published target-denominator analysis.
protein_state <- readRDS(file.path(m,"data/02_protein.rds"))
design_state <- readRDS(file.path(m,"data/01_design.rds"))
same_proxy_results <- list()
for(i in seq_len(nrow(mr_input))) {
  x <- mr_input[i,];pc <- design_state$proxy_candidates[i,]
  pr <- protein_state$protein_records
  pr <- pr[pr$protein==x$protein & pr$pos==pc$proxy_pos &
    ((pr$effect_allele==x$proxy_ea & pr$other_allele==x$proxy_oa) |
     (pr$effect_allele==x$proxy_oa & pr$other_allele==x$proxy_ea)),]
  stopifnot(nrow(pr)==1)
  bx_proxy <- pr$beta*ifelse(pr$effect_allele==x$proxy_ea,1,-1)
  fit_proxy <- TwoSampleMR::mr_wald_ratio(
    b_exp=bx_proxy,b_out=x$beta_Y,se_exp=pr$se,se_out=x$se_Y)
  same_proxy_results[[i]] <- data.frame(protein=x$protein,SNP=x$outcome_proxy_snp,
    effect_allele=x$proxy_ea,beta_X_proxy=bx_proxy,beta_Y_proxy=x$beta_Y,
    beta=fit_proxy$b,se=fit_proxy$se,p=fit_proxy$pval)
}
same_proxy_results <- do.call(rbind,same_proxy_results)
write.csv(mr_results,file.path(m,"results/steps/05_mr_results.csv"),row.names=FALSE)
write.csv(comparison,file.path(m,"results/steps/05_paper_comparison.csv"),row.names=FALSE)
write.csv(published_comparison, file.path(m, "results/published_vs_reproduced.csv"), row.names = FALSE)
write.csv(same_proxy_results,file.path(m,"results/steps/05_same_proxy_sensitivity.csv"),row.names=FALSE)
saveRDS(list(mr_results=mr_results,comparison=comparison,same_proxy_results=same_proxy_results),
        file.path(m,"data/05_results.rds"))
capture.output(sessionInfo(),file=file.path(m,"results/steps/05_sessionInfo.txt"))
real_pairs <- mr_input;real_mr_results <- mr_results
print(comparison,row.names=FALSE)
message("05 complete. Use View(mr_results), View(comparison), View(same_proxy_results).")

