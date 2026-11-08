# Statistical methods for the single-protein SNP association simulation

The workflow generates participant-level data, preprocesses the protein phenotype separately in each cohort, estimates one association per SNP per cohort, and combines cohort estimates. Sections 3–6 give the derivations, Section 7 describes the auxiliary tables, and Section 8 gives instructions for reopening the outputs. The [README](../README.md) provides the experimental rationale and results overview.

## 1. Running the analysis and locating the tables

From the standalone package root, run the following command. From the parent project, use the relative path in the README.

~~~sh
Rscript --vanilla run_all.R
~~~

In an interactive R session at the standalone root, use source("run_all.R"). The script writes eight CSV files, loads their corresponding R objects and opens four principal tables in the interactive data viewer. The counts below apply to the default configuration.

| Stage | File | R object | Unit and row count |
|---|---|---|---|
| 3. Participant generation | [individuals.csv](../data/simulated/individuals.csv) | raw_data | Participant, 3,301 |
| 3. Variant settings | [simulation_variants.csv](../results/tables/simulation_variants.csv) | simulation_variants | SNP, 20 |
| 4. Protein preprocessing | [individuals_processed.csv](../data/processed/individuals_processed.csv) | processed | Participant, 3,301 |
| 4. Covariate model | [covariate_coefficients.csv](../results/tables/covariate_coefficients.csv) | covariate_coefficients | Coefficient within cohort, 14 |
| 5. Single-SNP regression | [cohort_associations.csv](../results/tables/cohort_associations.csv) | cohort_results | SNP within cohort, 40 |
| 6. Meta-analysis | [meta_associations.csv](../results/tables/meta_associations.csv) | meta_results | Pooled SNP estimate, 20 |
| 6. Weight contributions | [meta_contributions.csv](../results/tables/meta_contributions.csv) | meta_contributions | SNP contribution within cohort, 40 |
| 7. Scale comparison | [scale_comparison.csv](../results/tables/scale_comparison.csv) | comparison | First-SNP estimate per outcome scale and cohort, six |

Use View(raw_data), View(processed), View(cohort_results) and View(meta_results) to inspect the principal tables. The row unit changes from participant to association estimate after preprocessing. Section 8.1 shows how to reopen saved outputs.

## 2. Objective, notation and parameters

### 2.1 Objective and notation

For participant $i$ in cohort $k\in\{A,B\}$ and SNP $j$, the sequence is

$$
\begin{aligned}
\{P_{ik},G_{ijk},C_{ik}\}&\longrightarrow\{\hat\beta_{jk},SE_{jk},p_{jk}\}\\
&\longrightarrow\{\hat\beta_{j,\mathrm{meta}},SE_{j,\mathrm{meta}},p_{j,\mathrm{meta}}\}.
\end{aligned}
$$

$P_{ik}$ denotes positive relative protein abundance, $G_{ijk}$ the effect-allele dosage and $C_{ik}$ the covariate vector. $\hat\beta$ is the estimated association slope for the processed protein phenotype, $SE$ is its standard error and $p$ is the P-value for testing a zero slope. IL23R_simulated and sim_snp_001 are synthetic identifiers.

### 2.2 Simulation parameters

| Setting | Default | Role |
|---|---:|---|
| seed | 20261021 | Fixes the pseudorandom sequence |
| cohort_sizes | A = 2,481, B = 820 | Uses the post-QC analysis counts reported under “Study participants” |
| n_snps | 20 | Number of single-SNP tests per cohort |
| effect_log_scale | 0.30 | First SNP coefficient in the log-abundance generating model |
| noise_sd | 0.70 | Standard deviation of the generation error |
| int_offset | 0.5 | Plotting position offset for the inverse normal transformation |

The simulation distributions, effect sizes, seed and INT offset are project choices. The entry script specifies Mersenne–Twister, Inversion and Rejection and reads the parameters from [config.R](../config.R). The cohort sizes follow Sun et al. (2018), Methods, “Study participants”. [Primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC6697541/).

## 3. Generating participant-level data

The function simulate_cohort() generates the data written to data/simulated/individuals.csv.

### 3.1 Genotypes

Draw a generating A-allele frequency for each SNP:

$$
f_j\sim U(0.10,0.45).
$$

Both cohorts use the same generating frequencies. Conditional on $f_j$, draw independently across individuals and SNPs:

$$
G_{ijk}\sim\mathrm{Binomial}(2,f_j).
$$

Counts 0, 1 and 2 represent the number of A alleles. Observed cohort frequencies vary by sampling. A and G are fixed simulation labels. This model omits LD and imputation uncertainty. Imputed allele dosages may take fractional values between 0 and 2.

### 3.2 Covariates

Draw age from $U(18,70)$, sex from $\mathrm{Bernoulli}(0.5)$, processing time from $\mathrm{Bernoulli}(0.2)$, and PC1–PC3 independently from $N(0,1)$. Processing time is coded 0 for at most one day and 1 for more than one day. The PC columns are simulated adjustment variables independent of genotype, rather than components estimated from a genotype matrix. Chance sample correlations remain possible.

### 3.3 Protein abundance

Generate log abundance and then exponentiate it:

$$
L_{ik}=2+0.30G_{i1k}+h(C_{ik})+\varepsilon_{ik},
$$

$$
\varepsilon_{ik}\sim N(0,0.70^2),\qquad P_{ik}=\exp(L_{ik}).
$$

The covariate contribution is

$$
\begin{aligned}
h(C_{ik})={}&0.01(\mathrm{age}_{ik}-45)+0.12\,\mathrm{sex}_{ik}
+0.15\,\mathrm{processing\_time}_{ik}\\
&+0.08\,\mathrm{PC1}_{ik}-0.05\,\mathrm{PC2}_{ik}
+0.03\,\mathrm{PC3}_{ik}.
\end{aligned}
$$

Exponentiation ensures positive abundance. Holding the covariates fixed, one additional A allele at SNP 1 raises generated log abundance by 0.30 and multiplies abundance by $\exp(0.30)\approx1.35$. This prespecified signal serves as a positive control. The other 19 SNPs have zero coefficients in the generating model and can still yield nominal associations by chance.

### 3.4 Output fields

| Quantity | Field | Interpretation |
|---|---|---|
| Participant and cohort | id, cohort | Keys linking the raw and processed records |
| $G_{ijk}$ | sim_snp_001 through sim_snp_020 | Counts of effect allele A |
| $C_{ik}$ | age, sex, processing_time, PC1, PC2, PC3 | Covariates used for generation and adjustment |
| $P_{ik}$ | protein_raw | Positive simulated abundance |
| $0.30G_{i1k}$ | genetic_log_component | Known genetic contribution on the log scale |
| $h(C_{ik})$ | covariate_log_component | Known covariate contribution on the log scale |
| $\varepsilon_{ik}$ | simulation_noise | Realised generation error |

Up to floating-point precision, log(protein_raw) equals 2 plus the three saved generating contributions. The simulator knows these contributions for each participant. They are distinct from estimates obtained by fitting the simulated data.

## 4. Protein preprocessing within cohorts

The function preprocess_cohort() performs all three operations separately in A and B, preserving the sequence described in Sun et al., Methods, “Genome-wide association study”. Sun's ancestry covariates came from multidimensional scaling. The simulated PC columns and the specific INT plotting position below are implementation choices. [Primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC6697541/).

### 4.1 Natural log transformation

$$
L_{ik}=\ln(P_{ik}).
$$

The log transformation compresses large positive abundances and expresses multiplicative changes additively. In this simulation it recovers the generating log scale. Normality of log abundance in a measured dataset depends on its distribution.

### 4.2 Covariate regression and residualisation

Within each cohort fit an intercept and the six covariates:

$$
L_{ik}=a_k+\gamma_k^{\mathsf T}C_{ik}+e_{ik}.
$$

The fitted value is

$$
\hat L_{ik}=\hat a_k+\hat\gamma_k^{\mathsf T}C_{ik}.
$$

The residual is

$$
\hat e_{ik}=L_{ik}-\hat L_{ik}.
$$

Positive residuals lie above the fitted covariate mean. The residual carries the remaining abundance variation, including genetic signal and noise. The fit can absorb some variation shared with covariates. The coefficients $\gamma_k$ describe the covariates, and SNPs enter the subsequent association regressions.

### 4.3 Rank-based inverse normal transformation

Rank the residuals from smallest to largest within cohort, assigning average ranks to ties. For residual rank $r_{ik}$ and cohort size $n_k$, define

$$
u_{ik}=\frac{r_{ik}-a}{n_k+1-2a},\qquad a=0.5,\qquad Z_{ik}=\Phi^{-1}(u_{ik}).
$$

Here $\Phi$ is the standard normal CDF, and the default plotting position is $(r_{ik}-0.5)/n_k$. With four distinct residuals the probabilities are 0.125, 0.375, 0.625 and 0.875, yielding approximately −1.15, −0.32, 0.32 and 1.15. The offset avoids probabilities of zero and one. The scores preserve rank order while changing spacing and need not have sample variance exactly one. The source describes rank inverse normalisation without establishing the specific offset used here.

### 4.4 Output fields

The 3,301 processed rows retain their original order and raw columns. Cohort membership determines both the ranking group and probability denominator.

| Operation | Field | Quantity |
|---|---|---|
| Raw abundance | protein_raw | $P_{ik}$ |
| Log transformation | protein_log | $L_{ik}$ |
| Covariate prediction | protein_fitted | $\hat L_{ik}$ |
| Residualisation | protein_residual | $\hat e_{ik}$ |
| Within-cohort ranking | residual_rank | $r_{ik}$ |
| Probability mapping | rank_probability | $u_{ik}$ |
| Normal quantile mapping | protein_INT | $Z_{ik}$, response for Section 5 |

## 5. Single-SNP association testing

The function associate_cohort() fits each SNP separately within each cohort:

$$
Z_{ik}=\alpha_{jk}+\beta_{jk}G_{ijk}+\eta_{ijk}.
$$

The slope estimates a marginal association in INT units per additional A allele. Each fit uses one SNP predictor and an intercept. This is a base R implementation of the additive regression idea reported by Sun et al., Methods, “Genome-wide association study”. Sun used SNPTEST v2.5.2 with expected allele dosages to account for imputation uncertainty. [Primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC6697541/).

For a single fit, let $g_i=G_{ijk}$ and $z_i=Z_{ik}$ and suppress $j,k$. Define

$$
S_{gg}=\sum_i(g_i-\bar g)^2,
$$

$$
\hat\beta=\frac{\sum_i(g_i-\bar g)(z_i-\bar z)}{S_{gg}},\qquad
\hat\alpha=\bar z-\hat\beta\bar g.
$$

The numerator describes joint variation in genotype and response, and $S_{gg}$ measures genotype variation. A monomorphic SNP gives $S_{gg}=0$ and triggers the variance check.

### 5.1 Standard error

The association model residual is

$$
\hat\eta_i=z_i-\hat\alpha-\hat\beta g_i.
$$

This residual is on the INT scale, whereas the residual in Section 4.2 is on the log-abundance scale. The conventional homoskedastic OLS variance and standard error are

$$
\hat\sigma^2=\frac{\sum_i\hat\eta_i^2}{n-2},\qquad
SE(\hat\beta)=\sqrt{\frac{\hat\sigma^2}{S_{gg}}}.
$$

The residual degrees of freedom are $n-2$ because the model estimates an intercept and a slope. The script checks both slope and standard error against the values returned by lm().

### 5.2 Two-sided P-value

For $H_0:\beta=0$ use

$$
t=\frac{\hat\beta}{SE(\hat\beta)},\qquad
p=2\left[1-F_{t_{n-2}}(|t|)\right],
$$

where $F_{t_{n-2}}$ is the Student t CDF. This describes the mathematical probability computed by the standard lm() summary, which evaluates the tail numerically. Classical exact inference relies on the regression error assumptions. The implementation uses the conventional test after residualisation and INT without an additional finite-sample calibration step.

### 5.3 Relationship to joint covariate adjustment

Residualising the outcome and regressing it on the original genotype can give a different slope from a joint regression on genotype and covariates. The Frisch–Waugh–Lovell equivalence for linear regression requires residualising both outcome and genotype. Here INT further changes the outcome nonlinearly after the covariate regression.

### 5.4 Output fields

Each of the 40 rows in cohort_associations.csv represents one regression, identified by cohort and SNP.

| Quantity | Field |
|---|---|
| Genotype and response means | mean_g, mean_z |
| Genotype sum of squares and genotype–response cross-product | s_gg, s_gz |
| Slope and intercept | beta, intercept |
| Residual sum of squares and degrees of freedom | sse, residual_df |
| Residual variance estimate and slope standard error | residual_variance, se |
| Signed t statistic and two-sided P-value | t_statistic, p |
| Participants in the regression | n |

The table also records protein, cohort, SNP, effect_allele, other_allele and the sample effect-allele frequency:

$$
EAF=\frac{\sum_i g_i}{2n}.
$$

The simulated data are complete, so the sample count is constant across SNPs within each cohort. Participant-specific association residuals can be reconstructed using the saved INT response, genotype, intercept and slope. Their squared sum is already recorded in sse.

## 6. Fixed-effect inverse-variance meta-analysis

The function meta_analyse() groups results by SNP. It requires two distinct cohort records, identical allele labels across cohorts and positive standard errors. For each SNP the inverse-variance weight is

$$
w_{jk}=1/SE_{jk}^2,
$$

and the combined association and standard error are

$$
\hat\beta_{j,\mathrm{meta}}=
\frac{\sum_k w_{jk}\hat\beta_{jk}}{\sum_k w_{jk}},\qquad
SE_{j,\mathrm{meta}}=\sqrt{1/\sum_k w_{jk}}.
$$

The method gives greater weight to more precise estimates and assumes independent cohort estimates of a common effect on a comparable scale. It pools summary estimates after separate preprocessing and regression. Sun et al. used METAL for fixed-effect inverse-variance pooling, as described in Methods, “Meta-analysis and statistical significance”. [Primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC6697541/).

The pooled two-sided P-value uses a normal approximation:

$$
z_{\mathrm{meta}}=\hat\beta_{\mathrm{meta}}/SE_{\mathrm{meta}},\qquad
p_{\mathrm{meta}}=2\Phi(-|z_{\mathrm{meta}}|).
$$

This differs from the cohort t tests. The implementation checks allele labels and stops on inconsistency. Real-data harmonisation can require strand checks and effect-allele reorientation. When allele counting is reversed for a biallelic variant, the slope changes sign and EAF becomes $1-EAF$.

### 6.1 Output fields

| Field in meta_associations.csv | Calculation or interpretation |
|---|---|
| beta | Sum of the cohort beta_contribution values for the SNP |
| sum_weights | Sum of inverse-variance weights |
| se | Square root of the reciprocal of sum_weights |
| z_statistic | beta divided by se |
| p | Two-sided normal tail probability |
| n | Sum of cohort sample counts, 3,301 |
| same_direction | Equality of the two signs as returned by sign() |
| both_nominal | Both cohort P-values strictly below 0.05 |

The 20 rows retain protein, SNP and allele identifiers. The two flags document criteria related to Sun's discovery procedure. The simulation reports all rows and omits the original study's full multiple-testing decision rule. The pooled P-value derives from the pooled estimate and SE.

## 7. Auxiliary tables

### 7.1 Variant generation settings

simulation_variants.csv contains one record per SNP. SNP identifies the genotype column, effect_allele and other_allele give the synthetic labels, generating_frequency stores $f_j$, and true_beta_log stores 0.30 for SNP 1 and zero for the other SNPs. The generating frequency is a model parameter and eaf in the association table is a sample estimate.

### 7.2 Covariate coefficients

covariate_coefficients.csv contains 14 rows, one intercept and six coefficients in each cohort. The columns cohort and term identify a fitted coefficient. estimate holds $\hat a_k$ or a component of $\hat\gamma_k$. The remaining columns give se, t_statistic, p, n and residual_df. The covariate model has $n_k-7$ residual degrees of freedom at full rank. Use the corresponding cohort's intercept and coefficients with the participant's covariates to reconstruct protein_fitted.

### 7.3 Weight contributions

meta_contributions.csv contains 40 rows and retains each cohort's association and allele identifiers.

| Field | Calculation |
|---|---|
| weight | $1/\mathrm{se}^2$ |
| sum_weights | Sum of the two weights for the SNP |
| weight_fraction | weight divided by sum_weights |
| weighted_beta | weight multiplied by beta |
| beta_contribution | weight_fraction multiplied by beta |

The two weight fractions sum to one for each SNP. Sum the two beta_contribution values to recover the pooled beta, or sum weighted_beta and divide by sum_weights.

### 7.4 Scale comparison for the positive control

scale_comparison.csv contains six rows, covering protein_raw, protein_log and protein_INT for each cohort. It records cohort, outcome, beta, se and p. The raw and log models include the six covariates jointly. The INT model uses the preprocessed outcome and genotype only.

The generating coefficient 0.30 is conditional and applies to log abundance. The estimated meta-analysis beta of about 0.4111 is measured in residual INT units per A allele. The INT depends on ranks across participants within each cohort. No fixed conversion factor links these coefficients, and the three scale-comparison fits differ in model specification as well as response scale.

## 8. Reopening outputs and reporting

### 8.1 Read saved tables without rerunning

Run this code from the standalone package root. In the parent project, set m to the relative module directory.

~~~r
m <- "."
raw_data <- read.csv(file.path(m, "data/simulated/individuals.csv"))
processed <- read.csv(file.path(m, "data/processed/individuals_processed.csv"))
cohort_results <- read.csv(file.path(m, "results/tables/cohort_associations.csv"))
meta_results <- read.csv(file.path(m, "results/tables/meta_associations.csv"))
simulation_variants <- read.csv(file.path(m, "results/tables/simulation_variants.csv"))
covariate_coefficients <- read.csv(file.path(m, "results/tables/covariate_coefficients.csv"))
meta_contributions <- read.csv(file.path(m, "results/tables/meta_contributions.csv"))
comparison <- read.csv(file.path(m, "results/tables/scale_comparison.csv"))

View(raw_data)
View(processed)
View(cohort_results)
View(meta_results)
~~~

Use View() for the auxiliary tables or head() for a console preview.

### 8.2 Inspect a particular calculation

~~~r
View(processed[, c("id", "cohort", "protein_raw", "protein_log",
  "protein_fitted", "protein_residual", "residual_rank",
  "rank_probability", "protein_INT")])
View(subset(cohort_results, SNP == "sim_snp_001"))
View(subset(meta_contributions, SNP == "sim_snp_001"))
~~~

The [preprocessing figure](../results/figures/preprocessing.png) displays the raw, log, residual and INT distributions. The script saves the tables and figure automatically and overwrites identically named generated files on rerun. Save a separate copy when retaining outputs for a parameter comparison.

### 8.3 Code and sources

| File | Role |
|---|---|
| [config.R](../config.R) | Seed, cohort sizes, SNP count and simulation settings |
| [run_all.R](../run_all.R) | Runs the four stages and saves outputs |
| [R/functions.R](../R/functions.R) | Implements generation, preprocessing, regression and meta-analysis |
| [config_used.R](../results/config_used.R) | Parameters recorded at runtime |
| [sessionInfo.txt](../results/sessionInfo.txt) | R version and runtime information |
| [Reproducibility report](../REPRODUCIBILITY_REPORT.md) | Run checks and documentation audit |

The parent project also provides RUN_BETA.R and its compatibility aliases beta_results, protein_data and scale_comparison. The standalone PowerShell wrapper uses the same run_all.R and resolves Rscript through PATH.

The simulation covers selected statistical steps. Protein measurement, assay and genetic QC, ancestry estimation, phasing and imputation fall outside its scope. The verification checks the base R implementation against its formulas and saved baseline. The project has no numerical validation against SNPTEST or METAL. The downstream MR exercise belongs to a separate module.

Sun BB et al. *Genomic atlas of the human plasma proteome*. Nature 558, 73–79 (2018). Methods, “Study participants”, “Genome-wide association study” and “Meta-analysis and statistical significance”. [DOI 10.1038/s41586-018-0175-2](https://doi.org/10.1038/s41586-018-0175-2).

Zheng J et al. *Phenome-wide Mendelian randomization mapping the influence of the plasma proteome on complex diseases*. Nature Genetics 52, 1122–1131 (2020). This reference provides context for the separate downstream exercise. [DOI 10.1038/s41588-020-0682-6](https://doi.org/10.1038/s41588-020-0682-6).
