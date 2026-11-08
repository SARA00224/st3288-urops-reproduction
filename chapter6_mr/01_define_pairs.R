# Run the numbered scripts in order. Objects remain visible in Environment.
m <- local({
  paths <- Filter(is.character, lapply(sys.frames(), function(frame) frame$ofile))
  script <- if (length(paths)) tail(paths, 1L)[[1L]] else
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(script) == 1L) dirname(normalizePath(script, mustWork = TRUE)) else getwd()
})
dir.create(file.path(m, "data"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(m, "results/steps"), recursive=TRUE, showWarnings=FALSE)

# 1. Read source rows exported directly from the paper's supplementary workbook.
# Provenance, source sheet and Excel row are preserved under sources/.
paper_instruments <- read.csv(file.path(m,"sources/paper_instruments.csv"),
                              stringsAsFactors=FALSE,check.names=FALSE)
paper_results <- read.csv(file.path(m,"sources/paper_mr_results.csv"),
                          stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(paper_instruments)==2,nrow(paper_results)==2)
paper_instruments <- paper_instruments[match(c("IL23R","RSPO3"),paper_instruments$Exposure),]
paper_results <- paper_results[match(c("IL23R","RSPO3"),paper_results$Exposure),]

# 2. Define the analysis. Do not include published MR answers in pair_design.
pair_design <- data.frame(
  protein=paper_instruments$Exposure, probe=paper_instruments$Platform_id,
  target_snp=paper_instruments$SNP, chr=paper_instruments$CHR_SNP,
  target_pos=paper_instruments$POS_SNP,
  target_ea=paper_instruments$Effect_allele,
  target_oa=paper_instruments$Other_allele,
  outcome="Crohn disease",outcome_id="ieu-a-12",
  n_case=17897,n_control=33977,build="GRCh37",ancestry="European",
  stringsAsFactors=FALSE)
input_mode <- match.arg(Sys.getenv("MR_INPUT_MODE", "full"), c("full", "excerpt"))
input_dir <- file.path(m, "data", if (input_mode == "full") "raw" else "excerpts")
pair_design$protein_file <- file.path(input_dir,
  c("GCST90241623_GRCh37.tsv.gz", "GCST90242572_GRCh37.tsv.gz"))
pair_design$metadata_file <- paste0(pair_design$protein_file, "-meta.yaml")
outcome_file <- file.path(input_dir, "ieu-a-12.vcf.gz")
required <- c(pair_design$protein_file, pair_design$metadata_file, outcome_file)
if (!all(file.exists(required))) stop("Missing input files: ",
  paste(basename(required[!file.exists(required)]), collapse = ", "))
message("Input mode: ", input_mode)

# 3. A curated shortlist from the previous investigation, NOT a genome-wide LD search.
# Outcome beta and published MR results are not used by the validation in step 04.
proxy_candidates <- data.frame(
  protein=c("IL23R","RSPO3"),proxy_snp=c("rs80174646","rs2503322"),
  proxy_pos=c(67708155L,127457260L),proxy_ea=c("T","G"),proxy_oa=c("G","A"))
ld_r2_threshold <- 0.8
saveRDS(list(pair_design=pair_design,proxy_candidates=proxy_candidates,
             outcome_file=outcome_file,ld_r2_threshold=ld_r2_threshold,input_mode=input_mode),
        file.path(m,"data/01_design.rds"))
pair_design_export <- pair_design
relative_input <- file.path("data", basename(input_dir))
pair_design_export$protein_file <- file.path(relative_input, basename(pair_design$protein_file))
pair_design_export$metadata_file <- file.path(relative_input, basename(pair_design$metadata_file))
write.csv(pair_design_export, file.path(m, "results/steps/01_pair_design.csv"), row.names = FALSE)
print(pair_design[,c("protein","target_snp","outcome_id")],row.names=FALSE)
message("01 complete. Inspect pair_design, paper_instruments and paper_results.")

