# Numbered scripts expose intermediate objects directly in Environment.
m <- local({
  paths <- Filter(is.character, lapply(sys.frames(), function(frame) frame$ofile))
  script <- if (length(paths)) tail(paths, 1L)[[1L]] else
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(script) == 1L) dirname(normalizePath(script, mustWork = TRUE)) else getwd()
})
dir.create(file.path(m,"data"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(m,"results/steps"),recursive=TRUE,showWarnings=FALSE)
tmp_storage <- file.path(m,"data/.tmp")
dir.create(tmp_storage,recursive=TRUE,showWarnings=FALSE)

stopifnot(requireNamespace("R.utils",quietly=TRUE),requireNamespace("data.table",quietly=TRUE),requireNamespace("yaml",quietly=TRUE))
if(!file.exists(file.path(m,"data/01_design.rds"))) stop("Run 01_define_pairs.R first.")
design_state <- readRDS(file.path(m,"data/01_design.rds"))
pair_design <- design_state$pair_design;proxy_candidates <- design_state$proxy_candidates
paper_instruments <- read.csv(file.path(m,"sources/paper_instruments.csv"),stringsAsFactors=FALSE)

# 1. Inspect the raw protein fields and check the cohort, probe and coordinate build.
metadata_audit <- data.frame();source_audit <- data.frame()
protein_preview <- list();protein_records <- list();target_exposure <- list()
for(i in seq_len(nrow(pair_design))) {
  x <- pair_design[i,];pc <- proxy_candidates[proxy_candidates$protein==x$protein,]
  meta <- yaml::read_yaml(x$metadata_file)
  stopifnot(meta$genome_assembly==x$build,meta$samples[[1]]$sample_size==3301,
            any(grepl(x$probe,unlist(meta$trait_description),fixed=TRUE)),
            "European" %in% unlist(meta$samples[[1]]$sample_ancestry_category))
  metadata_audit <- rbind(metadata_audit,data.frame(protein=x$protein,gwas_id=meta$gwas_id,
    probe=x$probe,build=meta$genome_assembly,n=meta$samples[[1]]$sample_size))
  message("Reading ",x$protein," protein summary...")
  protein_raw <- data.table::fread(x$protein_file,
    select=c("variant_id","chromosome","base_pair_location","effect_allele",
             "other_allele","beta","standard_error","log(P)","p_value"),data.table=FALSE,tmpdir=tmp_storage)
  protein_preview[[x$protein]] <- head(protein_raw,10)

  # 2. Extract the target and candidate-proxy positions, retaining true raw beta.
  # Proxy protein records are only needed by the separate same-proxy comparison in 05.
  keep <- protein_raw$chromosome==x$chr &
    protein_raw$base_pair_location %in% c(x$target_pos,pc$proxy_pos)
  z <- protein_raw[keep,,drop=FALSE]
  z <- data.frame(protein=x$protein,raw_variant_id=as.character(z$variant_id),
    chr=z$chromosome,pos=z$base_pair_location,effect_allele=toupper(z$effect_allele),
    other_allele=toupper(z$other_allele),beta=z$beta,se=z$standard_error,
    p=z$p_value,n=meta$samples[[1]]$sample_size)
  protein_records[[i]] <- z

  # 3. Identify the target by GRCh37 chr/position AND its allele pair.
  er <- z[z$pos==x$target_pos &
    ((z$effect_allele==x$target_ea & z$other_allele==x$target_oa) |
     (z$effect_allele==x$target_oa & z$other_allele==x$target_ea)),]
  stopifnot(nrow(er)==1)
  er$target_snp <- x$target_snp
  er$raw_effect_allele <- er$effect_allele;er$raw_beta <- er$beta
  # 4. Same allele: keep beta; swapped alleles: reverse beta. SE and P unchanged.
  er$sign_multiplier <- ifelse(er$effect_allele==x$target_ea,1,-1)
  er$beta <- er$raw_beta*er$sign_multiplier
  er$effect_allele <- x$target_ea;er$other_allele <- x$target_oa
  target_exposure[[i]] <- er
  source_audit <- rbind(source_audit,data.frame(protein=x$protein,
    rows_read=nrow(protein_raw),rows_retained=nrow(z)))
  rm(protein_raw);gc(verbose=FALSE)
}
protein_records <- do.call(rbind,protein_records)
target_exposure <- do.call(rbind,target_exposure)
stopifnot(all(is.finite(protein_records$beta)),all(protein_records$se>0),
          all(protein_records$p>=0 & protein_records$p<=1))
# Source-value comparison does not supply or overwrite the extracted beta.
exposure_check <- data.frame(protein=target_exposure$protein,
  beta_read=target_exposure$beta,paper_beta=paper_instruments$Beta,
  se_read=target_exposure$se,paper_se=paper_instruments$SE)
saveRDS(list(protein_records=protein_records,target_exposure=target_exposure,
  protein_preview=protein_preview,metadata_audit=metadata_audit,source_audit=source_audit,
  exposure_check=exposure_check),file.path(m,"data/02_protein.rds"))
write.csv(target_exposure,file.path(m,"results/steps/02_target_exposure.csv"),row.names=FALSE)
print(target_exposure[,c("protein","target_snp","raw_effect_allele","raw_beta",
                         "effect_allele","beta","se")],row.names=FALSE)
message("02 complete: protein fields -> target extraction -> allele-direction alignment.")

