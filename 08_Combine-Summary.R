##############################################################
# SCRIPT 8: Combine summary statistics from meta-analysis
##############################################################
# Setup
rm(list=ls())
library(tidyverse)

##############################################################
# Meta-Analysis Results
load_meta <- function(folder){
  results <- read_tsv(file=paste0("../IL6_Data/Processing/Meta_Output/METAANALYSIS_",file,".tbl"))
  results <- results %>% 
    mutate(pFDR = p.adjust(`P-value`, method="fdr")) %>% 
    dplyr::select(MarkerName, Effect, StdErr,
                  `P-value`, pFDR, Direction,
                  HetISq, HetPVal)
  return(results)
}

# IL-6 EWAS - Base model
results <- load_meta(file="il6_base")
colnames(results) <- c(
  "cpg", "base_meta_beta", "base_meta_se",
  "base_meta_p", "base_meta_padj_fdr",
  "base_het_dir", "base_het_i2", "base_het_p"
)
il6_results <- results

# IL-6 EWAS - Smoking sensitivity
results <- load_meta(file="il6_ext1")
colnames(results) <- c(
  "cpg", "ext1_meta_beta", "ext1_meta_se",
  "ext1_meta_p", "ext1_meta_padj_fdr",
  "ext1_het_dir", "ext1_het_i2", "ext1_het_p"
)
il6_results <- full_join(il6_results, results)

# IL-6 EWAS - IDOLext sensitivity
results <- load_meta(file="il6_ext3")
colnames(results) <- c(
  "cpg", "ext3_meta_beta", "ext3_meta_se",
  "ext3_meta_p", "ext3_meta_padj_fdr",
  "ext3_het_dir", "ext3_het_i2", "ext3_het_p"
)
il6_results <- full_join(il6_results, results)

# IL-6 EWAS - CRP sensitivity
results <- load_meta(file="il6_ext2")
colnames(results) <- c(
  "cpg", "ext2_meta_beta", "ext2_meta_se",
  "ext2_meta_p", "ext2_meta_padj_fdr",
  "ext2_het_dir", "ext2_het_i2", "ext2_het_p"
)

# CRP EWAS - base model
results <- load_meta(file="crp_base")
colnames(results) <- c(
  "cpg", "crp_base_meta_beta", "crp_base_meta_se",
  "crp_base_meta_p", "crp_base_meta_padj_fdr",
  "crp_base_het_dir", "crp_base_het_i2", "crp_base_het_p"
)
il6_results <- full_join(il6_results, results)

# CRP EWAS - IL-6 sensitivity
results <- load_meta(file="crp_ext2")
colnames(results) <- c(
  "cpg", "crp_ext2_meta_beta", "crp_ext2_meta_se",
  "crp_ext2_meta_p", "crp_ext2_meta_padj_fdr",
  "crp_ext2_het_dir", "crp_ext2_het_i2", "crp_ext2_het_p"
)
il6_results <- full_join(il6_results, results)

##############################################################
# Cohort Specific Results
# Function
read_cohort <- function(cohort, model){
  results <- read_tsv(file=paste0("../IL6_Data/Processing/Meta_Input/OUT_",
                model, "-", toupper(cohort), ".tsv")) %>% 
    dplyr::select(cpg, N, beta, se, p)
  return(results)
}

# LLS
# IL6 EWAS - Base model
results <- read_cohort("lls", "il6_base")
colnames(results) <- c(
  "cpg", "base_N_LLS", "base_beta_LLS", "base_SE_LLS",
  "base_p_LLS", "base_padj_fdr_LLS")
il6_results <- left_join(il6_results, results)

# IL6 EWAS - Smoking
results <- read_cohort("lls", "il6_ext1")
colnames(results) <- c(
  "cpg", "ext1_N_LLS", "ext1_beta_LLS", "ext1_SE_LLS",
  "ext1_p_LLS", "ext1_padj_fdr_LLS")
il6_results <- left_join(il6_results, results)

# IL6 EWAS - CRP
results <- read_cohort("lls", "il6_ext2")
colnames(results) <- c(
  "cpg", "ext2_N_LLS", "ext2_beta_LLS", "ext2_SE_LLS",
  "ext2_p_LLS", "ext2_padj_fdr_LLS")
il6_results <- left_join(il6_results, results)

# IL6 EWAS - IDOLext
results <- read_cohort("lls", "il6_ext3")
colnames(results) <- c(
  "cpg", "ext3_N_LLS", "ext3_beta_LLS", "ext3_SE_LLS",
  "ext3_p_LLS", "ext3_padj_fdr_LLS")
il6_results <- left_join(il6_results, results)

# CRP EWAS - Base
results <- read_cohort("lls", "crp_base")
colnames(results) <- c(
  "cpg", "crp_base_N_LLS", "crp_base_beta_LLS", "crp_base_SE_LLS",
  "crp_base_p_LLS", "crp_base_padj_fdr_LLS")
il6_results <- left_join(il6_results, results)

# CRP EWAs - IL-6
results <- read_cohort("lls", "crp_ext2")
colnames(results) <- c(
  "cpg", "crp_ext2_N_LLS", "crp_ext2_beta_LLS", "crp_ext2_SE_LLS",
  "crp_ext2_p_LLS", "crp_ext2_padj_fdr_LLS")
il6_results <- left_join(il6_results, results)

# KORA
# IL6 EWAS - Base
results <- read_cohort("kora", "il6_base")
colnames(results) <- c(
  "cpg", "base_N_KORA", "base_beta_KORA", "base_SE_KORA",
  "base_p_KORA", "base_padj_fdr_KORA")
il6_results <- left_join(il6_results, results)

# IL6 EWAS - Smoking
results <- read_cohort("kora", "il6_ext1")
colnames(results) <- c(
  "cpg", "ext1_N_KORA", "ext1_beta_KORA", "ext1_SE_KORA",
  "ext1_p_KORA", "ext1_padj_fdr_KORA")
il6_results <- left_join(il6_results, results)

# IL6 EWAS - CRP
results <- read_cohort("kora", "il6_ext2")
colnames(results) <- c(
  "cpg", "ext2_N_KORA", "ext2_beta_KORA", "ext2_SE_KORA",
  "ext2_p_KORA", "ext2_padj_fdr_KORA")
il6_results <- left_join(il6_results, results)

# IL6 EWAS - IDOLext
results <- read_cohort("kora", "il6_ext3")
colnames(results) <- c(
  "cpg", "ext3_N_KORA", "ext3_beta_KORA", "ext3_SE_KORA",
  "ext3_p_KORA", "ext3_padj_fdr_KORA")
il6_results <- left_join(il6_results, results)

# CRP EWAS - Base
results <- read_cohort("kora", "crp_base")
colnames(results) <- c(
  "cpg", "crp_base_N_KORA", "crp_base_beta_KORA", "crp_base_SE_KORA",
  "crp_base_p_KORA", "crp_base_padj_fdr_KORA")
il6_results <- left_join(il6_results, results)

# CRP EWAS - IL-6
results <- read_cohort("kora", "crp_ext2")
colnames(results) <- c(
  "cpg", "crp_ext2_N_KORA", "crp_ext2_beta_KORA", "crp_ext2_SE_KORA",
  "crp_ext2_p_KORA", "crp_ext2_padj_fdr_KORA")
il6_results <- left_join(il6_results, results)

# NTR
# IL6 EWAS - Base
results <- read_cohort("ntr", "il6_base")
colnames(results) <- c(
  "cpg", "base_N_NTR", "base_beta_NTR", "base_SE_NTR",
  "base_p_NTR", "base_padj_fdr_NTR")
il6_results <- left_join(il6_results, results)

# IL6 EWAS - Smoking
results <- read_cohort("ntr", "il6_ext1")
colnames(results) <- c(
  "cpg", "ext1_N_NTR", "ext1_beta_NTR", "ext1_SE_NTR",
  "ext1_p_NTR", "ext1_padj_fdr_NTR")
il6_results <- left_join(il6_results, results)

# IL6 EWAS - CRP
results <- read_cohort("ntr", "il6_ext2")
colnames(results) <- c(
  "cpg", "ext2_N_NTR", "ext2_beta_NTR", "ext2_SE_NTR",
  "ext2_p_NTR", "ext2_padj_fdr_NTR")
il6_results <- left_join(il6_results, results)

# IL6 EWAS - IDOLext
results <- read_cohort("ntr", "il6_ext3")
colnames(results) <- c(
  "cpg", "ext3_N_NTR", "ext3_beta_NTR", "ext3_SE_NTR",
  "ext3_p_NTR", "ext3_padj_fdr_NTR")
il6_results <- left_join(il6_results, results)

# CRP EWAS - Base
results <- read_cohort("ntr", "crp_base")
colnames(results) <- c(
  "cpg", "crp_base_N_NTR", "crp_base_beta_NTR", "crp_base_SE_NTR",
  "crp_base_p_NTR", "crp_base_padj_fdr_NTR")
il6_results <- left_join(il6_results, results)

# CRP EWAS - IL6
results <- read_cohort("ntr", "crp_ext2")
colnames(results) <- c(
  "cpg", "crp_ext2_N_NTR", "crp_ext2_beta_NTR", "crp_ext2_SE_NTR",
  "crp_ext2_p_NTR", "crp_ext2_padj_fdr_NTR")
il6_results <- left_join(il6_results, results)

##############################################################
# Methylation Statistics
load("../IL6_Data/Cohort_Results/ALL_ln-il6_base.Rdata")

limma_base <- limma_base %>% 
  dplyr::select(cpg, starts_with("cpg_mean"), starts_with("cpg_sd"))

il6_results <- left_join(il6_results, limma_base,by="cpg")

# Annotation - hg19
manifest_hg19 <- read_tsv("../Shared_Data/Manifests/HM450.hg19.manifest.tsv.gz")

anno <- manifest_hg19 %>% dplyr::select(
    cpg = probeID, cpg_chr_hg19 = CpG_chrm, cpg_start_hg19 = CpG_beg,
    cpg_end_hg19 = CpG_end, cpg_strand = probe_strand, gene_HGNC) %>% 
  mutate(cpg_chr_hg19 = substr(cpg_chr_hg19,4,5))
anno <- anno %>% dplyr::filter(cpg %in% limma_base$cpg)

# Annotation - hg38
manifest_hg38 <- read_tsv("../Shared_Data/Manifests/HM450.hg38.manifest.tsv.gz")

manifest_hg38 <- manifest_hg38 %>% dplyr::select(cpg = Probe_ID,
    cpg_chr_hg38 = CpG_chrm, cpg_start_hg38 = CpG_beg, cpg_end_hg38 = CpG_end) %>% 
  mutate(cpg_chr_hg38 = substr(cpg_chr_hg38,4,5))

anno <- left_join(anno, manifest_hg38,by="cpg")

# Annotation - ROADMAP
manifest_roadmap <- read_tsv("../Shared_Data/Manifests/HM450.hg19.REMC.chromHMM.tsv.gz")

manifest_roadmap <- manifest_roadmap %>% 
  dplyr::select(cpg = probeID, E062)

anno <- left_join(anno, manifest_roadmap, by="cpg")

il6_results <- left_join(il6_results, anno)

##############################################################
# Saving Results
save(il6_results, file="../IL6_Data/Results/IL6_EWAS_Results-Full.Rdata")

# Suppl. Tab. 1
write_csv(il6_results %>% arrange(cpg_chr_hg19,
                                  cpg_start_hg19), 
          file="../IL6_Data/Tables/ST1.csv")

##############################################################

