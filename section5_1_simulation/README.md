# Section 5.1 SNP–protein association simulation

## Research objective and scope

This NUS Undergraduate Research Opportunities Programme (UROP) exercise demonstrates how individual-level genotypes, protein abundance and covariates are converted into SNP–protein association summary statistics. It generates two synthetic subcohorts, preprocesses the protein phenotype within each, fits one additive model per SNP per subcohort, and combines the two estimates for each SNP. The label IL23R_simulated identifies a single illustrative phenotype. All participant records and variant identifiers are synthetic, so the results describe the simulated system and provide no evidence about real IL23R biology.

The full derivations, notation and links between equations and output fields are in [Statistical methods](docs/STATISTICAL_METHODS.md).

## Relationship to Sun et al. (2018)

Sun BB et al. *Genomic atlas of the human plasma proteome*. Nature 558, 73–79 (2018). The section names below refer to the Methods. [DOI 10.1038/s41586-018-0175-2](https://doi.org/10.1038/s41586-018-0175-2) · [Public full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC6697541/).

| Category | What this exercise uses |
|---|---|
| Reported design and methods | “Study participants” describes two non-overlapping INTERVAL subcohorts initially comprising 2,731 and 831 participants, with 2,481 and 820 remaining after genetic quality control. “Genome-wide association study” describes within-subcohort natural log transformation, regression on age, sex, blood draw to processing delay (≤1 day versus >1 day), and the first three ancestry components derived by multidimensional scaling, followed by rank inverse normalisation of residuals and additive tests on allelic dosages using SNPTEST v2.5.2. “Meta-analysis and statistical significance” describes fixed-effect inverse-variance pooling with METAL. |
| Simplified implementation | The simulation uses the reported post-QC sample sizes and preprocessing sequence. Base R linear regression on generated allele counts replaces SNPTEST on imputed dosages, and direct inverse-variance calculations replace METAL. The INT plotting position and tie rule below are explicit conventions of this implementation. |
| Simulation design choices | The SNP count and identifiers, frequency distribution, independence assumptions, covariate distributions, abundance model and coefficients, random seed and choice of a single positive control were specified for this exercise. |

The original significance strategy combined a meta-analysis threshold of $1.5\times10^{-11}$ with nominal evidence in both subcohorts and concordant effect directions. This exercise records the latter two indicators and reports every simulated SNP. The original multiple-testing threshold is cited as context for the discovery study.

## Running the experiment

Use R with its standard packages and PNG graphics support. R 4.4.1 was used for verification. All inputs are generated in the run, and no downloads, package installation or manual output directory setup are required.

From the root of the standalone Section 5.1 package:

~~~sh
Rscript --vanilla run_all.R
~~~

From the combined reproduction package root:

~~~sh
Rscript --vanilla section5_1_simulation/run_all.R
~~~

The entry script loads [config.R](config.R) and [R/functions.R](R/functions.R) from its own directory. It writes data and results beside the script, replacing previous generated files of the same name. In interactive R, run the equivalent path with source(), for example source("run_all.R") from the standalone root. The run loads the eight table objects listed below and opens the four main tables in the interactive data viewer.

## Simulation design

Generate 2,481 participants in subcohort A and 820 in subcohort B. Set seed 20261021 with the Mersenne–Twister, Inversion and Rejection random number settings. Draw one generating A-allele frequency $f_j\sim U(0.10,0.45)$ for each of 20 SNPs and share these frequencies across subcohorts. Conditional on $f_j$, draw $G_{ij}\sim\mathrm{Binomial}(2,f_j)$ independently across individuals and SNPs. Each $G_{ij}$ counts copies of the designated effect allele A. Observed cohort allele frequencies are estimated by $\overline G_j/2$ and can differ from $f_j$. Independent draws omit systematic linkage disequilibrium while allowing chance sample correlations.

| Covariate | Generation and role |
|---|---|
| Age | Uniform from 18 to 70 years. The abundance equation centres age at 45. |
| Sex | Bernoulli with probability 0.5, coded 0 and 1. |
| Processing time | Bernoulli with probability 0.2. Codes 0 and 1 represent at most one day and more than one day, respectively. |
| PC1–PC3 | Three independent standard normal variables representing adjustment terms. They are generated independently of genotype. |

Covariates and errors are drawn independently of genotypes. The same abundance model is used in both subcohorts. For participant $i$, generate relative abundance $P_i=\exp(L_i)$, where

$$
\begin{aligned}
L_i={}&2+0.30G_{i1}+0.01(\mathrm{Age}_i-45)+0.12\mathrm{Sex}_i
+0.15\mathrm{Processing}_i\\
&+0.08\mathrm{PC1}_i-0.05\mathrm{PC2}_i+0.03\mathrm{PC3}_i
+\varepsilon_i,\qquad \varepsilon_i\sim N(0,0.70^2).
\end{aligned}
$$

The known positive effect on sim_snp_001 provides a positive control for the direction and detectability of a deliberately introduced signal in this run. The remaining 19 SNPs have zero coefficients in the generating equation and illustrate finite-sample variation under a null design. The numerical choices specify an illustrative dataset. They were not estimated from Sun's participants or calibrated to achieve a target final coefficient. The generated PCs illustrate the mechanics of covariate adjustment without modelling ancestry or ancestry-related confounding.

## Protein preprocessing

Carry out each operation separately within A and B. This preserves the design of cohort-specific phenotype construction and association testing before pooling. It gives each cohort its own fitted covariate model and residual rank distribution.

1. **Natural log transformation.** Compute $L_i=\log P_i$. Taking logs compresses the upper tail of positive, right-skewed abundances and expresses multiplicative changes on an additive scale. Here the abundance was generated by exponentiation, so this step recovers the generated log scale.
2. **Covariate regression and residualisation.** Fit $L_i=\alpha_k+C_i^\mathsf T\gamma_k+e_i$ in cohort $k$, where $C_i$ contains the six covariates above. The residual $r_i=L_i-\widehat L_i$ is the variation left after subtracting the fitted covariate mean. Using it as the input to INT makes the ranks describe abundance variation after this adjustment. Residuals can still contain genetic signal, random error and variation unexplained by the fitted model.
3. **Rank-based inverse normal transformation (INT).** Rank the residuals within their cohort, using average ranks for ties. For rank $R_i$, cohort size $n_k$ and offset $a=0.5$, compute $u_i=(R_i-a)/(n_k+1-2a)=(R_i-0.5)/n_k$ and $Z_i=\Phi^{-1}(u_i)$, where $\Phi$ is the standard normal CDF. The offset keeps probabilities inside $(0,1)$. The transformation retains rank order and gives the marginal response an approximately normal scale while discarding original spacing. Ranking follows regression because regressing out covariates can change which individuals rank higher.

## SNP association testing

For each cohort $k$ and SNP $j$, fit the ordinary least squares model $Z_i=\alpha_{jk}+\beta_{jk}G_{ij}+e_i$. The predictor is the count of A alleles and the response is the transformed protein residual. The slope estimates the change in INT units per additional A allele. Separate fits yield one marginal association record per variant, the unit used in GWAS summary statistics.

The standard error quantifies sampling uncertainty under the fitted model. The cohort P-value uses a two-sided Student $t$ test for a zero slope. The slope, intercept, standard error and test statistic are derived in [Statistical methods, Section 5](docs/STATISTICAL_METHODS.md#5-single-snp-association-testing).

The covariate model operates on log abundance before the nonlinear transformation. This procedure yields a different estimand from a model that jointly regresses untransformed or log abundance on genotype and covariates.

## Fixed-effect meta-analysis

Pool cohort estimates for the same SNP using the fixed-effect model, which treats the two estimates as independent measurements of a common effect on the transformed scale. This preserves separate cohort preprocessing and demonstrates how to pool summary statistics when individual records are held separately. Stacking the individual data and fitting a pooled model can change the weighting and phenotype construction.

For cohort-specific estimates $b_k$ and standard errors $s_k$, let $w_k=1/s_k^2$. More precise estimates receive more weight.

$$
b_{\mathrm{meta}}=\frac{\sum_k w_k b_k}{\sum_k w_k}.
$$

The combined SE comes from the total weight, and the combined P-value uses a two-sided normal test. The code requires one record from each distinct cohort and identical effect and other allele labels across cohorts. The same_direction flag checks equality of the two estimated signs, and both_nominal indicates that both cohort P-values are below 0.05. They correspond to two criteria reported by Sun et al. and remain descriptive fields in this workflow. [Statistical methods, Section 6](docs/STATISTICAL_METHODS.md#6-fixed-effect-inverse-variance-meta-analysis).

## Results and interpretation

The verified run produced 3,301 individual records, 40 cohort-specific association records and 20 meta-analysis records. The positive control yielded the following values, rounded for display.

| Estimate for sim_snp_001 | N | Beta | SE | P-value |
|---|---:|---:|---:|---:|
| Cohort A | 2,481 | 0.408392 | 0.030111 | 1.783e-40 |
| Cohort B | 820 | 0.419134 | 0.052292 | 3.786e-15 |
| Fixed-effect meta-analysis | 3,301 | 0.411067 | 0.026094 | 6.535e-56 |

The generating coefficient 0.30 is the conditional change in log abundance per A allele. The meta-analysis coefficient 0.4111 is the association per A allele in units of the residual rank-normal score. Their numerical difference does not estimate bias because their scales and estimands differ. The scale comparison table illustrates this distinction. Its raw and log models include covariates directly, while its INT model uses the preprocessed response.

Three zero-effect SNPs (003, 009 and 020) have meta-analysis P-values below 0.05 in this run. Chance sample correlations and repeated testing can yield nominal results even when a SNP has no direct effect in the generating model. This single run illustrates that possibility and provides no estimate of the long-run false-positive rate.

| Generated file | R object | Content |
|---|---|---|
| data/simulated/individuals.csv | raw_data | 3,301 rows with identifiers, covariates, abundance, known generating components and 20 genotype columns |
| results/tables/simulation_variants.csv | simulation_variants | 20 rows with generating frequencies, allele labels and true log-scale coefficients |
| data/processed/individuals_processed.csv | processed | 3,301 rows with added log abundance, fitted value, residual, rank, probability and INT score |
| results/tables/covariate_coefficients.csv | covariate_coefficients | 14 fitted coefficient records, seven per cohort including intercept |
| results/tables/cohort_associations.csv | cohort_results | 40 slope, SE, P-value and diagnostic records |
| results/tables/meta_contributions.csv | meta_contributions | 40 cohort weight and weighted-effect contributions |
| results/tables/meta_associations.csv | meta_results | 20 pooled beta, SE, P-value, sample count and concordance records |
| results/tables/scale_comparison.csv | comparison | Six scale-specific regression records for the positive control |

Open [the distribution figure](results/figures/preprocessing.png) to inspect the four stages in each cohort. [config_used.R](results/config_used.R) and [sessionInfo.txt](results/sessionInfo.txt) record the applied settings and software environment. From the standalone root, load and verify the final table in R with:

~~~r
meta <- read.csv("results/tables/meta_associations.csv")
stopifnot(nrow(meta) == 20L, all(meta$n == 3301L))
meta[meta$SNP == "sim_snp_001", c("beta", "se", "p")]
~~~

See the [reproducibility report](REPRODUCIBILITY_REPORT.md) for checks against the earlier Section 5.1 results.

## Limitations

The exercise omits protein measurement and assay quality control, real genotype and sample quality control, phasing, imputation, ancestry inference, realistic LD and genome-wide discovery. The 20 markers and single simulated phenotype represent a small illustration. The run also does not assess sensitivity to heterogeneous cohort effects or demonstrate the calibration or power of the procedure across repeated simulations. The documentation describes selected steps in Sun et al. and the corresponding base R calculations, without asserting numerical equivalence to SNPTEST or METAL.
