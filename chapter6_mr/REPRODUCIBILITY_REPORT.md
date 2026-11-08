# Chapter 6 reproducibility report

## Audit and inputs

The audit used the existing local Supplementary_Tables.xlsx downloaded with Zheng et al. (2020), the two original Sun protein TSVs, the original ieu-a-12 VCF, and the previously cached Ensembl and 1000 Genomes reference records. The workbook SHA-256 matched the earlier provenance record. The three compressed GWAS inputs matched the recorded file sizes and MD5 values. No input GWAS version was replaced and no new online association or LD request was made. [Publisher article and Supplementary Tables](https://www.nature.com/articles/s41588-020-0682-6).

The original workbook was read directly. S1 records were located by protein, probe, Sun study, original SNP, cis status and Tier 1. S7 rows were independently matched by protein, probe, study, original SNP, cis status, phenotype, outcome ID 12, outcome N = 51,874 and WR method. Each filter returned a unique record. The reported MR fields are cells X60:AA60 for IL23R and X97:AA97 for RSPO3, with beta, SE and P in Y, Z and AA. The original exposure records are S1 rows 1021 and 1712. The workbook extraction preserves the missing historical Outcome_SNP and LD_r2 fields. S9 row 30 is the separate two-instrument IL23R IVW result and was excluded from the comparison. The freshly extracted values match the existing project reference export. A final Chapter 6 manuscript was unavailable for direct comparison.

The full Crohn's disease VCF reports a conflicting StudyType header label (Continuous). The outcome is treated as case–control using the study metadata and case and control counts, with the supplied ES used on its log-odds scale. The archived response URLs identify the pairwise LD population as 1000GENOMES:phase_3:EUR and the coordinate build as GRCh37. Exact original retrieval times and genotype-request URLs are unavailable in the saved metadata. The reference JSON files preserve variant IDs and phased genotypes.

## Environment and execution

R 4.4.1 on Windows x64, with data.table 1.18.2.1, R.utils 2.13.0, yaml 2.3.12, jsonlite 2.0.0 and TwoSampleMR 0.7.13. See the [session record](results/steps/05_sessionInfo.txt). The runs used the C locale and the existing project package library. The full run used:

~~~powershell
$env:MR_INPUT_MODE = 'full'
Rscript --vanilla papers/01_plasma_proteome_MR/real_pair_reproduction/run_mr.R
~~~

A separate copy containing only the code, cached sources and selected raw association lines was run from its own root:

~~~powershell
$env:MR_INPUT_MODE = 'excerpt'
Rscript --vanilla run_mr.R
~~~

Both runs completed all five stages. The gzip reader emitted temporary file access warnings during the full run. It used the successful file operation and returned complete results. The full run verified zero direct hits in the available VCF for both instruments. The excerpt run rechecked allele alignment and LD with the same inputs and reproduced the six checked CSVs byte for byte. These were 02_target_exposure, 04_ld_validation, 04_mr_input, 05_mr_results, 05_same_proxy_sensitivity, and published_vs_reproduced. Main and sensitivity MR outputs also matched the pre-edit outputs byte for byte.

The changes replace workstation paths with paths relative to the entry script, add the excerpt input mode and emit a comparison based on the fresh workbook extraction. The numbered analysis, proxy candidates and Wald estimator were retained. The full local input files and original workbook were preserved. The transferable bundle includes selected association records and cached public reference data. It excludes the full GWAS downloads and publisher workbook. Obtain the workbook from the publisher for independent inspection. Confirm source redistribution terms before depositing the bundle publicly.

The packaged ZIP was also extracted and rerun. A first extraction into the sandbox-managed temporary directory encountered a permission error when the gzip reader renamed an intermediate file. Extraction into a writable workspace directory succeeded, and the same six comparison CSVs matched the full-input run byte for byte. Run the bundle from a writable directory.

## Numerical verification

The manually computed Wald ratio and first-order SE matched the package result. The absolute differences from the workbook for beta were 4.44089209850063e−15 and 2.77555756156289e−16 for IL23R and RSPO3. Their SE differences were 1.90438094876955e−08 and 6.69544444609205e−08. These are observed differences. The workbook and cached inputs do not establish their exact cause.

| Protein | Published beta | Reproduced beta | Published SE | Reproduced SE | Published P | Reproduced P | Absolute P difference |
|---|---:|---:|---:|---:|---:|---:|---:|
| IL23R | 1.83896428571429 | 1.83896428571429 | 0.0707459333371429 | 0.0707459523809524 | 5.80114354300137e−149 | 5.80220033460784e−149 | 1.05679160646551e−152 |
| RSPO3 | 0.242606296296296 | 0.242606296296296 | 0.0446565997122222 | 0.0446566666666667 | 5.55049876784864e−08 | 5.55075222543042e−08 | 2.53457581780886e−12 |

## Remaining checks

Publication values and instrument identities are verified. The comparison is marked partially verified because S7 does not identify the historical outcome SNP or proxy LD for either row. Exact historical input precision and the cause of the small SE differences remain unresolved. A current API-to-cache comparison, fresh package installation, cross-platform execution and final thesis Table 6.1 verification were not performed.
