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

stopifnot(requireNamespace("R.utils",quietly=TRUE),requireNamespace("data.table",quietly=TRUE))
if(!file.exists(file.path(m,"data/02_protein.rds"))) stop("Run 02_read_gwas.R first.")
design_state <- readRDS(file.path(m,"data/01_design.rds"))
protein_state <- readRDS(file.path(m,"data/02_protein.rds"))
pair_design <- design_state$pair_design;proxy_candidates <- design_state$proxy_candidates
target_exposure <- protein_state$target_exposure

# 1. Read the disease VCF header: this is GWAS summary data, not individual genotypes.
con <- gzfile(design_state$outcome_file,"rt");vcf_header <- character()
repeat {
  line <- readLines(con,n=1,warn=FALSE)
  if(!length(line)) stop("Missing VCF column header.")
  vcf_header <- c(vcf_header,line)
  if(startsWith(line,"#CHROM")) break
}
close(con)
stopifnot(any(grepl("assembly=HG19/GRCh37",vcf_header,fixed=TRUE)))
vcf_raw <- data.table::fread(design_state$outcome_file,skip="#CHROM",data.table=FALSE,tmpdir=tmp_storage)
stopifnot(ncol(vcf_raw)==10,names(vcf_raw)[10]=="ieu-a-12")
names(vcf_raw)[1] <- "chr"
disease_preview <- head(vcf_raw,10)

# 2. Assign effect allele=ALT and beta=ES; decode FORMAT by field NAME.
disease_data <- data.frame(SNP=vcf_raw$ID,chr=as.integer(vcf_raw$chr),pos=vcf_raw$POS,
  effect_allele=toupper(vcf_raw$ALT),other_allele=toupper(vcf_raw$REF),
  beta=NA_real_,se=NA_real_,LP=NA_real_,eaf=NA_real_,n=NA_real_)
field_map <- c(ES="beta",SE="se",LP="LP",AF="eaf",SS="n")
for(layout in unique(vcf_raw$FORMAT)) {
  idx <- which(vcf_raw$FORMAT==layout)
  field_names <- strsplit(layout,":",fixed=TRUE)[[1]]
  values <- data.table::tstrsplit(vcf_raw[[10]][idx],":",fixed=TRUE)
  if(!all(c("ES","SE","LP") %in% field_names)) stop("Missing disease beta/SE/P fields.")
  for(field in intersect(names(field_map),field_names)) {
    raw_values <- values[[match(field,field_names)]]
    raw_values[raw_values=="."] <- NA_character_
    disease_data[idx,field_map[[field]]] <- as.numeric(raw_values)
  }
}
# 3. LP stores -log10(P). ES already stores log odds, so do NOT log-transform ES.
disease_data$p <- 10^(-disease_data$LP)
stopifnot(all(is.finite(disease_data$beta)),all(disease_data$se>0),
          all(disease_data$p>=0 & disease_data$p<=1),
          all(disease_data$n[!is.na(disease_data$n)]==51874))
# StudyType=Continuous in this VCF header conflicts with the original case-control design.

# 4. Look for each original target, then retain a distinct candidate if it is absent.
candidate_outcome <- list();match_audit <- list();nearby_records <- list()
for(i in seq_len(nrow(pair_design))) {
  x <- pair_design[i,];pc <- proxy_candidates[proxy_candidates$protein==x$protein,]
  id_hit <- which(disease_data$SNP==x$target_snp)
  coord_hit <- which(disease_data$chr==x$chr & disease_data$pos==x$target_pos &
    ((disease_data$effect_allele==x$target_ea & disease_data$other_allele==x$target_oa) |
     (disease_data$effect_allele==x$target_oa & disease_data$other_allele==x$target_ea)))
  match_audit[[i]] <- data.frame(protein=x$protein,target_snp=x$target_snp,
    disease_id_hits=length(id_hit),disease_coordinate_hits=length(coord_hit))
  if(length(id_hit)+length(coord_hit)>0) stop("Direct target now available: review proxy plan.")
  cr <- disease_data[disease_data$chr==x$chr & disease_data$pos==pc$proxy_pos &
    disease_data$SNP==pc$proxy_snp &
    ((disease_data$effect_allele==pc$proxy_ea & disease_data$other_allele==pc$proxy_oa) |
     (disease_data$effect_allele==pc$proxy_oa & disease_data$other_allele==pc$proxy_ea)),]
  stopifnot(nrow(cr)==1)
  cr$protein <- x$protein;candidate_outcome[[i]] <- cr
  # Nearby records are inspectable in Environment, not separate output files.
  nr <- disease_data[disease_data$chr==x$chr & abs(disease_data$pos-x$target_pos)<=5000,]
  nr$protein <- x$protein;nr$distance_to_target <- nr$pos-x$target_pos
  nearby_records[[i]] <- nr[order(abs(nr$distance_to_target)),]
}
candidate_outcome <- do.call(rbind,candidate_outcome)
match_audit <- do.call(rbind,match_audit);nearby_records <- do.call(rbind,nearby_records)
# A candidate is not approved by closeness or disease significance: validate LD in 04.
disease_matches <- merge(match_audit,candidate_outcome,by="protein",sort=FALSE)
disease_matches$status <- if (identical(design_state$input_mode, "excerpt"))
  "Cached candidate; target absence checked in the full input run" else
  "Target absent; candidate awaiting independent LD/direction validation"
saveRDS(list(target_exposure=target_exposure,candidate_outcome=candidate_outcome,
  disease_matches=disease_matches,match_audit=match_audit,
  disease_preview=disease_preview,vcf_header=vcf_header),
  file.path(m,"data/03_matches.rds"))
write.csv(disease_matches,file.path(m,"results/steps/03_disease_matches.csv"),row.names=FALSE)
print(disease_matches[,c("protein","target_snp","disease_id_hits","SNP","beta","se")],
      row.names=FALSE)
message("03 complete: VCF fields -> disease estimates -> target search -> candidate matches.")

