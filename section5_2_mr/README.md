# Section 5.2 selected MR reproduction

## Objective and scope

This NUS UROP exercise reproduces the single-instrument Wald ratio calculations for IL23R and RSPO3 in relation to Crohn's disease. It starts from the cis Tier 1 instruments reported by Zheng et al. (2020) and uses public protein and disease GWAS summary statistics. The scope covers association lookup, allele alignment, LD proxy validation and numerical comparison. Genome-wide instrument selection, phenome-wide MR and colocalization are outside the workflow. [Zheng et al., Nature Genetics, including the Supplementary Tables download](https://www.nature.com/articles/s41588-020-0682-6).

## Sources and record locations

| Input | Provenance |
|---|---|
| Instruments and published MR results | Publisher's Supplementary Tables workbook. Table S1, sheet TableS1.Instruments, rows 1021 (IL23R) and 1712 (RSPO3). Table S7, sheet TableS7.cis-only-coloc, rows 60 and 97, columns X–AA (method, beta, SE and P). The [direct workbook extraction](sources/published_workbook_records.csv) records the workbook SHA-256, matching identifiers and source cells. |
| IL23R exposure | Sun et al. (2018), INTERVAL, probe IL23R.5088.175.3. GCST90241623_GRCh37.tsv.gz, decompressed line 3,995,932. [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90241623) |
| RSPO3 exposure | Sun et al. (2018), INTERVAL, probe RSPO3.13094.75.3. GCST90242572_GRCh37.tsv.gz, decompressed line 8,739,277. [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90242572) |
| Crohn's disease outcome | Liu et al. (2015), International IBD Genetics Consortium. OpenGWAS ieu-a-12, 17,897 cases and 33,977 controls. The two proxy records are in ieu-a-12.vcf.gz, decompressed lines 2,943 and 60,055. [OpenGWAS dataset](https://opengwas.io/datasets/ieu-a-12) |
| LD and direction | Cached Ensembl GRCh37 pairwise LD responses and 1000 Genomes Phase 3 reference genotypes. The signed LD check uses the 503 EUR samples identified in the panel. The [LD output](results/steps/04_ld_validation.csv) gives both exact pairwise URLs and EUR population parameters. |

Line numbers are one-based and include headers. The [record locator](sources/input_record_locations.csv) includes coordinates, roles and the two protein-proxy records used only in the existing sensitivity analysis. The [full-input manifest](sources/full_input_checksums.csv) identifies the exact compressed GWAS files by bytes and MD5. Input excerpts retain the original headers and selected lines. The protein metadata MD5 values identify the full downloads.

## Selected associations

Both pairs use a distinct disease-side proxy. EA is the effect allele, and beta values below are on the aligned scales used in the ratio.

| Protein | Original SNP, EA | Protein beta (SE) | Outcome SNP, EA | Disease beta (SE) | Ensembl EUR LD r² |
|---|---|---:|---|---:|---:|
| IL23R | rs11581607, A | −0.42 (0.0491) | rs80174646, T | −0.772365 (0.0297133) | 0.920472 |
| RSPO3 | rs2489623, C | 0.27 (0.0245) | rs2503322, G | 0.0655037 (0.0120573) | 0.938008 |

The full VCF was searched by original rsID and by GRCh37 chromosome, position and allele pair. Both searches returned zero direct matches for both tools. The S7 Outcome_SNP and LD_r2 cells are missing in the matched records. S1 lists rs11209026 in the IL23R instrument-validation Proxy_SNP column. This is a separate field and does not document the disease association used in the S7 ratio. The IL23R S7 row also contains a separate multi-cis result, which is reported as an IVW analysis in Table S9 row 30.

## Procedure

1. **Define pairs and read exposures.** Match the probe and instrument by GRCh37 coordinates and allele pair. Align to the Table S1 exposure allele. RSPO3 is recorded as A/C with beta −0.27 in the raw TSV, so reorientation to C/A gives +0.27. The IL23R exposure direction is retained.
2. **Read outcomes and validate LD.** Parse VCF FORMAT fields by name. ES is the ALT-allele log-odds coefficient, SE its standard error and LP the negative base-10 log P-value. Verify that the proposed target and proxy alleles are positively correlated using signed dosage and phased haplotype correlations. The signed haplotype correlations for A↔T and C↔G are approximately +0.95941 and +0.96857. Squared LD alone cannot determine effect direction. Both disease ALT alleles already match the validated orientation.
3. **Calculate and compare.** Use the original protein coefficient and the aligned disease proxy coefficient without LD scaling. Compute the Wald ratio $b_Y/b_X$ and first-order SE $s_Y/|b_X|$, which treats the exposure coefficient as fixed for the SE calculation. Compare the hand calculation with TwoSampleMR::mr_wald_ratio and use a two-sided normal test. Compare with the workbook only after estimating the ratio. Beta is on the Crohn's log-odds scale per unit of the GWAS-transformed protein phenotype.

## Running the analysis

R 4.4.1 was used with data.table, R.utils, yaml, jsonlite and TwoSampleMR. Exact versions are in the [R session record](results/steps/05_sessionInfo.txt). Install these packages into an R library on the search path before running. Package installation can require internet access. With packages and inputs in place, the workflow runs offline and needs no API token.

For a full-input run from the combined package root, place the original files under section5_2_mr/data/raw and use full mode:

~~~sh
Rscript --vanilla section5_2_mr/run_mr.R
~~~

To run only Section 5.2 from its section5_2_mr directory, use the excerpt mode. In PowerShell:

~~~powershell
$env:MR_INPUT_MODE = 'excerpt'
Rscript --vanilla run_mr.R
~~~

In a POSIX shell:

~~~sh
MR_INPUT_MODE=excerpt Rscript --vanilla run_mr.R
~~~

The excerpt contains the exact six selected GWAS records and their headers, plus the cached reference inputs in sources/. It reproduces the downstream estimates and LD checks. A full-file search is required to establish that the target is absent across the available disease VCF. To repeat the full search from the standalone root, place the three original compressed files and two protein metadata YAMLs in data/raw, check the file names and MD5 values against the manifest, and use MR_INPUT_MODE=full. The script defaults to full mode when the variable is unset.

run_mr.R executes 01_define_pairs.R through 05_run_mr_compare.R in order and creates output directories. Files under data/ named 01_design.rds through 05_results.rds are generated state. results/steps contains the successive matching, harmonization, LD and estimate tables. Use [04_mr_input.csv](results/steps/04_mr_input.csv) for the aligned input records and [published_vs_reproduced.csv](results/published_vs_reproduced.csv) for the publication comparison.

## Results and comparison

| Protein–phenotype | Published beta | Reproduced beta | Absolute beta difference | Published SE | Reproduced SE | Absolute SE difference |
|---|---:|---:|---:|---:|---:|---:|
| IL23R–Crohn's disease | 1.83896428571429 | 1.83896428571429 | 4.44e−15 | 0.0707459333371429 | 0.0707459523809524 | 1.90438e−08 |
| RSPO3–Crohn's disease | 0.242606296296296 | 0.242606296296296 | 2.78e−16 | 0.0446565997122222 | 0.0446566666666667 | 6.69544e−08 |

The published values were extracted directly from the original workbook. Beta agrees at the displayed precision. The SE and P differences are retained in the CSV. The same-proxy sensitivity output changes the exposure SNP and remains separate from the primary results.

## Limitations

The comparison verifies two published single-instrument numerical results. The historical outcome proxies remain unresolved. Candidate selection preceded this audit and was informed by the published values. Signed LD was checked independently using reference genotypes. The exercise does not independently establish the original causal conclusions. See the [verification report](REPRODUCIBILITY_REPORT.md) for the exact P-value comparison and remaining source questions.
