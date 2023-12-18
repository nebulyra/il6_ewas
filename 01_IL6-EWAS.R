##############################################################
# SCRIPT 1: An Epigenome-Wide Association Study (EWAS) between 
# IL-6 and whole blood DNA methylation (DNAm) levels in the 
# Leiden Longevity Study (LLS)
##############################################################
# Setup
rm(list=ls())
library(SummarizedExperiment)
library(haven)
library(tidyverse)
library(ggpubr)
library(DNAmArray)
library(minfi)
library(limma)

# Data Import

# A SummarizedExperiment object containing the DNA methylation data is loaded into the environment
load("../Shared_Data/LLS_betas.Rdata")

# IL-6 data is also loaded in to the environment
file_name <- "LLS_IOP1-450kDNAmethylation-AdiponectinLeptinIl6_20201009.sav"
lls_cytokines <- read_sav(paste0("../Shared_Data/", file_name))

##############################################################
# Data preprocessing
# Remove NA values
lls_cytokines <- lls_cytokines %>% filter(!is.na(IL6))

# Remove values below LOD
lls_cytokines <- lls_cytokines %>% filter(IL6>0)

# Check Normality
ggdensity(lls_cytokines$IL6)
ggqqplot(lls_cytokines$IL6)
shapiro.test(lls_cytokines$IL6)

# Natural log transformation
lls_cytokines$log_il6 <- log(lls_cytokines$IL6)

# 3IQR outlier removal
lower_limit <- quantile(lls_cytokines$log_il6, na.rm=TRUE)[[2]] - 
  (3 * IQR(lls_cytokines$log_il6, na.rm=TRUE))
upper_limit <- quantile(lls_cytokines$log_il6, na.rm=TRUE)[[4]] +
  (3 * IQR(lls_cytokines$log_il6, na.rm=TRUE))

lls_cytokines <- lls_cytokines %>% 
  mutate(
    log_il6 = ifelse(
      lls_cytokines$log_il6 < lower_limit | lls_cytokines$log_il6 > upper_limit, 
      NA, log_il6),
    IL6 = ifelse(
      lls_cytokines$log_il6 < lower_limit | lls_cytokines$log_il6 > upper_limit, 
      NA, IL6)
  )

# Check Normality
ggdensity(lls_cytokines$log_il6)
ggqqplot(lls_cytokines$log_il6)
shapiro.test(lls_cytokines$log_il6)

##############################################################
# Phenotype data extraction
# IDs
targets <- data.frame(uuid = rownames(colData(betas)),
                      run_id = colData(betas)$run_id)

# Add sex coded as 0=male, 1=female
targets <- targets %>% 
  mutate(
    sex = factor(
      case_when(
        colData(betas)$sex == "male" ~ 0,
        colData(betas)$sex == "female" ~ 1
      ), labels = c("Male", "Female")))

# Add smoking coded as 0=never, 1=ex, 2=current
targets <- targets %>% 
  mutate(
    smoking = factor(
      case_when(
        colData(betas)$smoking == "non-smoker" ~ 0,
        colData(betas)$smoking == "former-smoker" ~ 1,
        colData(betas)$smoking == "current smoker" ~ 2
      ), labels = c("Never", "Former", "Current")))

# hsCRP data
targets$crp <- colData(betas)$hscrp

# IDOL cell counts
counts <- read_tsv("../Shared_Data/LLS_countsUCB.tsv")
colnames(counts) <- c("run_id", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
targets <- left_join(targets, counts, by = "run_id")

# IDOL extended cell counts
load("../Shared_Data/LLS_deconvolution.RData")
counts_ext <- LLS.prop %>% rownames_to_column(var="uuid")
counts_ext <- counts_ext %>% dplyr::select(
  uuid,
  Mono_ext = Mono, 
  NK_ext = NK, 
  Bas_ext = Bas, 
  Bmem_ext = Bmem, 
  Bnv_ext = Bnv, 
  CD4mem_ext = CD4mem, 
  CD4nv_ext = CD4nv,
  CD8mem_ext = CD8mem, 
  CD8nv_ext = CD8nv, 
  Eos_ext = Eos, 
  Neu_ext = Neu,
  Treg_ext = Treg
)
targets <- left_join(targets, counts_ext, by="uuid")

# Technical covariates
targets$plate <- colData(betas)$sample_plate
targets$row <- colData(betas)$sentrix_position

# Recode row as numeric
targets <- targets %>% mutate(row = as.numeric(substr(row, 3, 3)))

# Join data
targets <- left_join(targets, lls_cytokines, by = "uuid")
targets <- zap_labels(targets)

##############################################################
# Exclusion criteria
# Sex chromosomes
betas <- betas[!(seqnames(betas) %in% c("chrX", "chrY")),]

# Encode Blacklist Regions
load("../Shared_Data/ENCODE_Blacklist/ENCODEBlacklist_CpGomit-450K.RData")
betas <- betas[!rownames(betas) %in% cpg_blacklist,]

# Zhou probes
maskProbes <- read_tsv("../Shared_Data/Manifests/HM450.hg19.manifest.tsv.gz")
maskProbes <- maskProbes[maskProbes$MASK_general == TRUE,]$probeID
betas <- betas[!rownames(betas) %in% maskProbes,]

# 3IQR outlier Removal
iqr_dnam <- apply(assay(betas), 1, function(x){
  iqr <- IQR(x, na.rm = TRUE)
  q1 <- quantile(x, na.rm=TRUE)[[2]]
  q3 <- quantile(x, na.rm=TRUE)[[4]]
  x <- ifelse((x <= q1 - (3*iqr) | x >= q3 + (3*iqr)), NA, x)
})
assay(betas) <- t(iqr_dnam)

# CpGs with >5% missingness
perc_na <- rowSums(is.na(iqr_dnam))*100/ncol(iqr_dnam)
betas <- betas[,perc_na <= 95]

# Samples with >5% missingness
perc_na <- colSums(is.na(iqr_dnam))*100/nrow(iqr_dnam)
betas <- betas[perc_na <= 95,]

# Subset phenotype data based on exclusions
targets <- targets %>% 
  filter(uuid %in% rownames(colData(betas)))

# Check order
rownames(targets) <- targets$uuid
identical(rownames(targets), rownames(colData(betas)))

# Add phenotype data to SummarizedExperiment
colData(betas) <- DataFrame(targets)

##############################################################
# Data preparation for analysis
# Subset complete data for base model variables
# Base model
betas <- subset(betas, 
                select = rowSums(is.na(colData(betas)[, c(1:13)])) == 0)


# Smoking Sensitivity Analysis
betas_smoke <- subset(betas, 
                      select = rowSums(is.na(colData(betas)[, c(1:14)])) == 0)

# CRP Sensitivity Analysis
betas_crp <- subset(betas, 
                    select = rowSums(is.na(colData(betas)[, c(1:13,15)])) == 0)

# IDOl Extended Sensitivity Analysis
betas_idol <- subset(betas, 
                     select = rowSums(is.na(colData(betas)[, c(1:13,16)])) == 0)

##############################################################
# EWAS
# Base Model
# Fit model
metadata(betas)$formula <- ~log_il6 + age + sex + CD8T + CD4T + NK + Bcell + Mono + plate + row
design <- model.matrix(metadata(betas)$formula, data=colData(betas))
fit <- lmFit(assay(betas), design)

# Store results
beta <- fit$coefficients[, 2]
SE <- fit$stdev.unscaled[, 2] * fit$sigma
t <- beta / SE
p <- 2 * pt(-abs(t), fit$df.residual)
N <- ncol(design) + fit$df.residual
cpg_mean <- apply(assay(betas), 1, mean)
cpg_se <- apply(assay(betas), 1, sd)

limma_base <- data.frame(cpg = rownames(fit$coefficients), 
                         cpg_mean, cpg_se, beta, SE, p, N)

# Save
save(limma_base, file='../IL6_Data/Cohort_Results/LLS/LLS_ln-il6_base.Rdata')

# Smoking Sensitivity Analysis
# Fit model
metadata(betas_smoke)$formula <- ~log_il6 + age + sex + smoking + CD8T + CD4T + NK + Bcell + Mono + plate + row
design <- model.matrix(metadata(betas_smoke)$formula, data=colData(betas_smoke))
fit <- lmFit(assay(betas_smoke), design)

# Store results
beta <- fit$coefficients[, 2]
SE <- fit$stdev.unscaled[, 2] * fit$sigma
t <- beta / SE
p <- 2 * pt(-abs(t), fit$df.residual)
N <- ncol(design) + fit$df.residual

limma_base <- data.frame(cpg = rownames(fit$coefficients), 
                         beta, SE, p, N)

# Save results
save(limma_base, file='../IL6_Data/Cohort_Results/LLS/LLS_ln-il6_ext1.Rdata')

# CRP Sensitivity Analysis
# Fit model
metadata(betas_crp)$formula <- ~log_il6 + age + sex + crp + CD8T + CD4T + NK + Bcell + Mono + plate + row
design <- model.matrix(metadata(betas_crp)$formula, data=colData(betas_crp))
fit <- lmFit(assay(betas_crp), design)

# Store results
beta <- fit$coefficients[, 2]
SE <- fit$stdev.unscaled[, 2] * fit$sigma
t <- beta / SE
p <- 2 * pt(-abs(t), fit$df.residual)
N <- ncol(design) + fit$df.residual

limma_base <- data.frame(cpg = rownames(fit$coefficients), 
                         beta, SE, p, N)

# Save results
save(limma_base, file='../IL6_Data/Cohort_Results/LLS/LLS_ln-il6_ext2.Rdata')

# IDOL Extended Sensitivity Analysis
# Fit model
metadata(betas_idol)$formula <- ~log_il6 + age + sex + plate + row + 
  Bas_ext + Bmem_ext + Bnv_ext + 
  CD4mem_ext + CD4nv_ext + CD8mem_ext + CD8nv_ext + 
  Eos_ext + Mono_ext + NK_ext + Treg_ext 
design <- model.matrix(metadata(betas_idol)$formula, 
                       data=colData(betas_idol))
fit <- lmFit(assay(betas_idol), design)

# Store results
beta <- fit$coefficients[, 2]
SE <- fit$stdev.unscaled[, 2] * fit$sigma
t <- beta / SE
p <- 2 * pt(-abs(t), fit$df.residual)
N <- ncol(design) + fit$df.residual

limma_base <- data.frame(cpg = rownames(fit$coefficients), 
                         beta, SE, p, N)

# Save results
save(limma_base, file='../IL6_Data/Cohort_Results/LLS/LLS_ln-il6_ext3.Rdata')

##############################################################






























