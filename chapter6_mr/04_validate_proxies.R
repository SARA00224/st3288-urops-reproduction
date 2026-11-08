# Run the numbered scripts in order. Objects remain visible in Environment.
m <- local({
  paths <- Filter(is.character, lapply(sys.frames(), function(frame) frame$ofile))
  script <- if (length(paths)) tail(paths, 1L)[[1L]] else
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(script) == 1L) dirname(normalizePath(script, mustWork = TRUE)) else getwd()
})
dir.create(file.path(m, "data"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(m, "results/steps"), recursive=TRUE, showWarnings=FALSE)

stopifnot(requireNamespace("jsonlite",quietly=TRUE),requireNamespace("data.table",quietly=TRUE))
if(!file.exists(file.path(m,"data/03_matches.rds"))) stop("Run 03_match_instruments.R first.")
design_state <- readRDS(file.path(m,"data/01_design.rds"))
matches <- readRDS(file.path(m,"data/03_matches.rds"))
pair_design <- design_state$pair_design;proxy_candidates <- design_state$proxy_candidates
# Use the raw external responses, NOT the previous pass/fail summary.
ld_api <- jsonlite::fromJSON(file.path(m,"sources/ensembl_proxy_checks.json"),simplifyVector=FALSE)
panel <- data.table::fread(file.path(m,"sources/1000G_phase3_sample_panel.tsv"),header=TRUE,fill=TRUE,data.table=FALSE)
eur_ids <- panel$sample[panel$super_pop=="EUR"]
stopifnot(length(eur_ids)==503)
ld_validation <- list();mr_input <- list();haplotype_counts <- list();dosage_tables <- list()
mapping_validation <- list()

for(i in seq_len(nrow(pair_design))) {
  x <- pair_design[i,];pc <- proxy_candidates[proxy_candidates$protein==x$protein,]
  # Verify cached build/coordinate/allele mappings for both distinct SNPs.
  for(kind in c("MAP","PROXY_MAP")) {
    mp <- ld_api[[paste0(x$protein,"_",kind)]]
    expected_snp <- if(kind=="MAP") x$target_snp else pc$proxy_snp
    expected_pos <- if(kind=="MAP") x$target_pos else pc$proxy_pos
    expected_alleles <- if(kind=="MAP") c(x$target_ea,x$target_oa) else c(pc$proxy_ea,pc$proxy_oa)
    stopifnot(mp$data$name==expected_snp)
    found <- FALSE
    for(mm in mp$data$mappings) {
      alleles <- strsplit(mm$allele_string,"/",fixed=TRUE)[[1]]
      if(mm$assembly_name=="GRCh37" && as.character(mm$seq_region_name)==as.character(x$chr) &&
         mm$start==expected_pos && all(expected_alleles %in% alleles)) found <- TRUE
    }
    stopifnot(found)
    mapping_validation[[length(mapping_validation)+1]] <- data.frame(
      protein=x$protein,SNP=expected_snp,chr=x$chr,pos=expected_pos,
      build="GRCh37",mapping_passed=found,source_url=mp$url)
  }
  api <- ld_api[[paste0(x$protein,"_LD")]]
  a <- api$data[[1]]
  stopifnot(a$variation1==x$target_snp,a$variation2==pc$proxy_snp,
            a$population_name=="1000GENOMES:phase_3:EUR")
  api_r2 <- as.numeric(a$r2)
  stopifnot(api_r2>=design_state$ld_r2_threshold)
  tg <- jsonlite::fromJSON(file.path(m,"sources",paste0(x$target_snp,"_genotypes.json")))$genotypes
  pg <- jsonlite::fromJSON(file.path(m,"sources",paste0(pc$proxy_snp,"_genotypes.json")))$genotypes
  tg$sample_id <- sub("^.*:","",tg$sample);pg$sample_id <- sub("^.*:","",pg$sample)
  j <- merge(tg[,c("sample_id","genotype")],pg[,c("sample_id","genotype")],
             by="sample_id",suffixes=c("_target","_proxy"))
  j <- j[j$sample_id %in% eur_ids,]
  stopifnot(nrow(j)==503,!anyDuplicated(j$sample_id))
  # Counts 0/1/2 for the proposed corresponding effect alleles.
  target_hap <- strsplit(j$genotype_target,"|",fixed=TRUE)
  proxy_hap <- strsplit(j$genotype_proxy,"|",fixed=TRUE)
  target_dosage <- numeric(nrow(j));proxy_dosage <- numeric(nrow(j))
  for(k in seq_len(nrow(j))) {
    stopifnot(length(target_hap[[k]])==2,length(proxy_hap[[k]])==2,
              all(target_hap[[k]] %in% c(x$target_ea,x$target_oa)),
              all(proxy_hap[[k]] %in% c(pc$proxy_ea,pc$proxy_oa)))
    target_dosage[k] <- sum(target_hap[[k]]==x$target_ea)
    proxy_dosage[k] <- sum(proxy_hap[[k]]==pc$proxy_ea)
  }
  signed_dosage_r <- cor(target_dosage,proxy_dosage)
  # Paired phased haplotypes: compute D, signed r and r^2 independently.
  ht <- unlist(target_hap,use.names=FALSE);hp <- unlist(proxy_hap,use.names=FALSE)
  pA <- mean(ht==x$target_ea);pB <- mean(hp==pc$proxy_ea)
  pAB <- mean(ht==x$target_ea & hp==pc$proxy_ea)
  D <- pAB-pA*pB
  signed_haplotype_r <- D/sqrt(pA*(1-pA)*pB*(1-pB))
  haplotype_r2 <- signed_haplotype_r^2
  stopifnot(signed_dosage_r>0,signed_haplotype_r>0,
            haplotype_r2>=design_state$ld_r2_threshold)
  # Different SNPs are paired through LD phase, not by equating their allele letters.
  er <- matches$target_exposure[matches$target_exposure$protein==x$protein,]
  yr <- matches$candidate_outcome[matches$candidate_outcome$protein==x$protein,]
  stopifnot(nrow(er)==1,nrow(yr)==1)
  if(yr$effect_allele==pc$proxy_ea && yr$other_allele==pc$proxy_oa) sign_y <- 1 else
    if(yr$effect_allele==pc$proxy_oa && yr$other_allele==pc$proxy_ea) sign_y <- -1 else
      stop("Unexpected proxy allele pair; no silent strand flip.")
  mr_input[[i]] <- data.frame(
    protein=x$protein,outcome=x$outcome,outcome_id=x$outcome_id,
    target_snp=x$target_snp,target_ea=x$target_ea,target_oa=x$target_oa,
    beta_X=er$beta,se_X=er$se,p_X=er$p,n_X=er$n,
    outcome_proxy_snp=pc$proxy_snp,proxy_ea=pc$proxy_ea,proxy_oa=pc$proxy_oa,
    raw_beta_Y=yr$beta,sign_multiplier_Y=sign_y,beta_Y=sign_y*yr$beta,
    se_Y=yr$se,p_Y=yr$p,proxy_eaf=if(sign_y==1)yr$eaf else 1-yr$eaf,
    n_Y=yr$n,n_case=x$n_case,n_control=x$n_control,
    API_r2=api_r2,haplotype_r2=haplotype_r2,signed_dosage_r=signed_dosage_r,
    pairing="Target exposure + validated outcome proxy; no r scaling")
  ld_validation[[i]] <- data.frame(protein=x$protein,target_snp=x$target_snp,
    proxy_snp=pc$proxy_snp,target_ea=x$target_ea,proxy_ea=pc$proxy_ea,
    EUR_n=nrow(j),API_r2=api_r2,D=D,signed_haplotype_r=signed_haplotype_r,
    haplotype_r2=haplotype_r2,signed_dosage_r=signed_dosage_r,
    dosage_r2=signed_dosage_r^2,threshold=design_state$ld_r2_threshold,
    passed=TRUE,source_url=api$url)
  hc <- as.data.frame(table(target_allele=ht,proxy_allele=hp),stringsAsFactors=FALSE)
  hc$protein <- x$protein;haplotype_counts[[i]] <- hc
  dosage_tables[[i]] <- data.frame(protein=x$protein,sample_id=j$sample_id,
                                  target_dosage=target_dosage,proxy_dosage=proxy_dosage)
}
mr_input <- do.call(rbind,mr_input);ld_validation <- do.call(rbind,ld_validation)
haplotype_counts <- do.call(rbind,haplotype_counts)
dosage_data <- do.call(rbind,dosage_tables)
mapping_validation <- do.call(rbind,mapping_validation)
saveRDS(list(mr_input=mr_input,ld_validation=ld_validation,
  haplotype_counts=haplotype_counts,dosage_data=dosage_data,
  mapping_validation=mapping_validation),file.path(m,"data/04_validation.rds"))
write.csv(mr_input,file.path(m,"results/steps/04_mr_input.csv"),row.names=FALSE)
write.csv(ld_validation,file.path(m,"results/steps/04_ld_validation.csv"),row.names=FALSE)
print(ld_validation[,c("protein","API_r2","haplotype_r2","signed_dosage_r","passed")],row.names=FALSE)
message("04 complete. Inspect mr_input, ld_validation, haplotype_counts and dosage_data.")


