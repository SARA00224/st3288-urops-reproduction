# ST3288 UROPS computational reproduction

This companion package contains the two R experiments in Chapter 5, Computational Exercises and Reproducibility. The thesis describes the reproduction materials in Section 5.3 and Supplementary Section S3. Each experiment has its own inputs, scripts, documentation and results.

| Thesis section | Experiment | Documentation and outputs |
|---|---|---|
| 5.1 | Generate SNP–protein association statistics from simulated participant-level data. The workflow covers protein preprocessing, single-SNP regression and fixed-effect meta-analysis. | [Experiment guide](section5_1_simulation/README.md) · [Statistical methods](section5_1_simulation/docs/STATISTICAL_METHODS.md) · [Verification report](section5_1_simulation/REPRODUCIBILITY_REPORT.md) · [Meta-analysis results](section5_1_simulation/results/tables/meta_associations.csv) |
| 5.2 | Recalculate two selected protein–Crohn's disease Wald ratios from public GWAS summary statistics, validate the outcome proxy alleles and compare the estimates with the published Supplementary Table 7. | [Experiment guide](section5_2_mr/README.md) · [Verification report](section5_2_mr/REPRODUCIBILITY_REPORT.md) · [Published and reproduced results](section5_2_mr/results/published_vs_reproduced.csv) |

## Run both experiments

Extract the package into a writable directory. Open a terminal in the directory containing this README. The reference environment was R 4.4.1 on Windows. Section 5.1 requires the standard R packages and PNG graphics support. Section 5.2 also requires data.table, R.utils, yaml, jsonlite and TwoSampleMR. If needed, run these commands in an R session:

~~~r
install.packages(c("data.table", "R.utils", "yaml", "jsonlite"),
                 repos = "https://cloud.r-project.org")
install.packages("TwoSampleMR",
                 repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
~~~

Use an R library visible to the Rscript executable used below. Exact versions from the original runs are recorded in the experiment session files. Package installation instructions for additional platforms are in the [TwoSampleMR documentation](https://mrcieu.github.io/TwoSampleMR/).

From the package root:

~~~sh
Rscript --vanilla run_all.R
~~~

The launcher runs Section 5.1 and then Section 5.2 in separate R sessions. Once the packages are installed, the default run is offline. Section 5.1 regenerates its synthetic participants. Section 5.2 uses the six source association records and the cached public reference inputs shipped with the package. It recomputes allele alignment, LD, the Wald ratios and the publication comparison. The bundled results are examples from the verified run. Rerunning the scripts overwrites matching generated files and creates any needed output directories.

## Run one experiment or repeat the full GWAS search

To run Section 5.1 alone from the package root, use Rscript --vanilla section5_1_simulation/run_all.R. To run Section 5.2 alone, follow its [experiment guide](section5_2_mr/README.md).

The Section 5.2 excerpt has only the association records needed for the selected calculations. To repeat the original full-file absence checks for the instrument SNPs, obtain the full source files listed in the [Section 5.2 guide](section5_2_mr/README.md), put them under section5_2_mr/data/raw and verify their checksums. Then run:

~~~sh
Rscript --vanilla run_all.R --full
~~~

The original publisher workbook is available through the [Zheng et al. article](https://www.nature.com/articles/s41588-020-0682-6). The project includes a [record-level extraction](section5_2_mr/sources/published_workbook_records.csv) with sheet names, row numbers and source cells. Study-specific data terms apply to third-party inputs.

## Verification

The combined launcher was verified from a clean writable directory with R 4.4.1 and the recorded packages available on the R library search path. Fourteen regenerated files matched the reference checksums, including the eight Section 5.1 CSVs, its preprocessing figure and five Section 5.2 result files.

## References


- Sun BB et al. *Genomic atlas of the human plasma proteome*. Nature 558, 73–79 (2018). [DOI](https://doi.org/10.1038/s41586-018-0175-2).
- Zheng J et al. *Phenome-wide Mendelian randomization mapping the influence of the plasma proteome on complex diseases*. Nature Genetics 52, 1122–1131 (2020). [DOI](https://doi.org/10.1038/s41588-020-0682-6).
