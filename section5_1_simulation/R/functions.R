# Functions use base R only. Synthetic dosages are observed 0/1/2 counts,
# a special case of dosage; imputation uncertainty and LD are not simulated.
simulate_cohort <- function(n, cohort, cfg, maf) {
  g <- sapply(maf, function(f) rbinom(n, 2, f))
  colnames(g) <- sprintf("sim_snp_%03d", seq_along(maf))
  d <- data.frame(id = paste0(cohort, "_", seq_len(n)), cohort = cohort,
    age = runif(n, 18, 70), sex = rbinom(n, 1, .5),
    processing_time = rbinom(n, 1, .2),
    PC1 = rnorm(n), PC2 = rnorm(n), PC3 = rnorm(n))
  simulation_noise <- rnorm(n, sd = cfg$noise_sd)
  d$protein_raw <- exp(2 + cfg$effect_log_scale * g[, 1] +
    .01 * (d$age - 45) + .12 * d$sex + .15 * d$processing_time +
    .08 * d$PC1 - .05 * d$PC2 + .03 * d$PC3 +
    simulation_noise)
  # Known simulation truth, not fitted contributions in a real dataset.
  d$genetic_log_component <- cfg$effect_log_scale * g[,1]
  d$covariate_log_component <- .01*(d$age-45) + .12*d$sex +
    .15*d$processing_time + .08*d$PC1 - .05*d$PC2 + .03*d$PC3
  d$simulation_noise <- simulation_noise
  cbind(d, as.data.frame(g))
}
preprocess_cohort <- function(d, offset = .5) {
  stopifnot(all(is.finite(d$protein_raw)), all(d$protein_raw > 0),
    !anyNA(d))
  d$protein_log <- log(d$protein_raw)
  fit <- lm(protein_log ~ age + sex + processing_time + PC1 + PC2 + PC3,
    data = d)
  d$protein_fitted <- fitted(fit)
  d$protein_residual <- residuals(fit)
  d$residual_rank <- rank(d$protein_residual, ties.method = "average")
  d$rank_probability <- (d$residual_rank - offset) /
    (nrow(d) + 1 - 2 * offset)
  coefficients <- coef(summary(fit))
  attr(d, "covariate_coefficients") <- data.frame(
    cohort=d$cohort[1], term=rownames(coefficients),
    estimate=coefficients[,1], se=coefficients[,2],
    t_statistic=coefficients[,3], p=coefficients[,4],
    n=nrow(d), residual_df=fit$df.residual, row.names=NULL)
  # Explicit teaching convention: (average rank - 0.5)/n.
  # Exact original Sun rank-offset convention remains to be verified.
  d$protein_INT <- qnorm(d$rank_probability)
  d
}
associate_cohort <- function(d, variants, protein) {
  do.call(rbind, lapply(variants, function(v) {
    g <- d[[v]]; z <- d$protein_INT
    stopifnot(all(is.finite(g)), all(g >= 0 & g <= 2), var(g) > 0)
    fit <- lm(z ~ g)
    s <- coef(summary(fit))["g", ]
    # Independent arithmetic check of OLS beta and SE.
    ssx <- sum((g - mean(g))^2)
    sgz <- sum((g - mean(g)) * (z - mean(z)))
    beta <- sgz / ssx
    sse <- sum(residuals(fit)^2)
    residual_df <- length(g) - 2
    residual_variance <- sse / residual_df
    se <- sqrt(residual_variance / ssx)
    stopifnot(isTRUE(all.equal(unname(s[1]), beta, tolerance = 1e-10)),
      isTRUE(all.equal(unname(s[2]), se, tolerance = 1e-10)))
    data.frame(protein = protein, cohort = d$cohort[1], SNP = v,
      effect_allele = "A", other_allele = "G", eaf = mean(g)/2,
      beta = unname(s[1]), se = unname(s[2]), p = unname(s[4]), n = nrow(d),
      mean_g=mean(g), mean_z=mean(z), s_gg=ssx, s_gz=sgz,
      intercept=unname(coef(fit)[1]), sse=sse, residual_df=residual_df,
      residual_variance=residual_variance, t_statistic=unname(s[3]))
  }))
}
meta_analyse <- function(tab) {
  do.call(rbind, lapply(split(tab, tab$SNP), function(x) {
    stopifnot(nrow(x) == 2, length(unique(x$cohort)) == 2,
      length(unique(x$effect_allele)) == 1,
      length(unique(x$other_allele)) == 1, all(x$se > 0))
    w <- 1/x$se^2
    b <- sum(w*x$beta)/sum(w); se <- sqrt(1/sum(w))
    data.frame(protein = x$protein[1], SNP = x$SNP[1],
      effect_allele = x$effect_allele[1], other_allele = x$other_allele[1],
      beta = b, se = se, p = 2*pnorm(-abs(b/se)), n = sum(x$n),
      same_direction = length(unique(sign(x$beta))) == 1,
      both_nominal = all(x$p < .05), sum_weights=sum(w), z_statistic=b/se)
  }))
}
