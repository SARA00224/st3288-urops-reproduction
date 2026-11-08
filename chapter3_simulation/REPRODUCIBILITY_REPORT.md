# Chapter 3 reproducibility report

Verification date: 6 November 2026. The experiment was executed with R 4.4.1 (UCRT) on Windows 11 x64, using the standard R packages and the PNG graphics device. The verification process used the C locale. The generated [session record](results/sessionInfo.txt) documents the runtime.

## Audit and changes

The audit compared the existing simulation scripts, Chinese method notes, generated tables, diagnostic figure, and Chapter 3 result snapshot dated 4 November 2026. The audit also checked the reported sample sizes, preprocessing, SNP association and meta-analysis against the Methods subsections “Study participants”, “Genome-wide association study” and “Meta-analysis and statistical significance”. [Sun et al. (2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6697541/).

The previous launcher relied on a workstation path and the parent project marker. The revised entry locates its companion scripts relative to the invoked file and writes outputs beside it. The Windows wrapper now resolves Rscript through PATH. The seed already existed, and the entry now explicitly sets the default random number algorithms used by the reference run. These changes leave the data-generating equation, models, seed and parameters intact. The original Chinese explanations remain in the parent project, with the previous module README preserved as README_ZH.md.

## Commands and checks

The first run used a separate directory containing only the simulation scripts. It started with no generated inputs or output directories. From the root of that standalone directory, the command is:

~~~sh
Rscript --vanilla run_all.R
~~~

The entry was also run by passing the script path from outside its directory and by sourcing it from R. The parent convenience script was sourced and its three compatibility objects were checked. The working directory remained the same during sourcing.

| Check | Observed result |
|---|---|
| Population and genotype columns | 3,301 unique participant records across cohorts A (2,481) and B (820), with 20 SNP columns |
| Preprocessing | Raw and processed participant order matched. Log abundance, residual subtraction and INT quantiles agreed with the saved fields. |
| Output tables | Eight CSVs were generated with 3,301 raw rows, 3,301 processed rows, 20 variant rows, 14 covariate coefficient rows, 40 cohort association rows, 40 weight contribution rows, 20 pooled rows and six scale comparison rows. |
| Association calculations | The existing arithmetic checks of slope and SE against the linear model passed for all 40 fits. |
| Meta-analysis | All 20 beta, SE and P-value triples agreed with an independent calculation using weighted.mean(), inverse-variance SE and the normal tail probability at relative tolerance 1e-12. Allele labels and both cohort concordance flags were checked. |
| Existing output baseline | SHA-256 hashes matched for all eight CSVs, the preprocessing PNG and the saved configuration against the pre-change baseline. |
| Chapter 3 result snapshot | The five listed SNPs and their meta-analysis beta, SE and P-values matched the 4 November result snapshot at relative tolerance 1e-12. |

The positive-control meta-analysis estimate was beta 0.411066755949578, SE 0.0260942661712653 and two-sided P = 6.53489941238052e-56. The 0.30 generating coefficient applies to log abundance. The estimated association applies to the residual INT scale.

## Statistical methods translation and equation audit

The original METHODS_EXPLAINED_ZH.md remains in the parent project. Its SHA-256 hash matches the copy in the 4 November snapshot. The separate [English statistical methods document](docs/STATISTICAL_METHODS.md) follows the original eight-section order, retains the derivations and output mappings, and uses relative paths for the standalone package.

| Check | Finding and treatment |
|---|---|
| Data generation and natural log | The two cohort sizes, allele-frequency and genotype draws, covariate distributions, log-abundance coefficients and error SD agree with the script. The saved log abundance was reconstructed from the generating components. |
| Covariate regression and INT | The per-cohort fitted values, residuals, average ranks, plotting positions and normal quantiles were verified from the output tables. The offset and tie rule are simulation implementation conventions. |
| OLS | The slope, intercept, residual variance, SE and signed t statistic agree with the saved regressions. The original two-sided P-value expression uses one minus the t CDF. The equivalent lower-tail computation used by the R model summary avoids subtraction from one in the extreme tail. The translation retains the mathematical expression and clarifies its numerical evaluation. |
| Fixed-effect pooling | Allele-label checks, inverse-variance weights, normal-tail P-values and concordance fields agree with the implementation and the prior run verification. |
| Scale comparison | The raw and log models adjust for covariates jointly, and the INT model uses the preprocessed phenotype. The documentation explicitly attributes the differences to both model specification and phenotype scale. |
| Formula rendering and navigation | All 21 display formulas and 47 inline expressions in the English methods document parsed in strict mode with local KaTeX. The README formulas and relative links also passed. The files use GitHub Markdown math delimiters. Rendering on a published GitHub page has yet to be checked. |

No substantive discrepancy requiring a change to a model, parameter or reported result was found. The translation expands the original plain-text covariate contribution into a display equation and typesets the original formula notation. The numerical audit reran the same R workflow and reproduced the existing meta-analysis results.

## Outstanding checks

The final Chapter 3 manuscript was unavailable for a line-by-line review. Before submission, check its reported numbers, scale labels and caption against the verified results. The standalone bundle contains the English documentation, base R source, and synthetic outputs. Publication and the choice of a version-specific citation remain for the author. The verification covers the local Windows/R environment. Other operating systems, R versions and graphics devices were not tested, and binary graphics hashes may differ across environments.
