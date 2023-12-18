##############################################################
# SCRIPT 2: An Epigenome-Wide Association Study (EWAS) between 
# CRP and whole blood DNA methylation (DNAm) levels in the 
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

# Data import
# A SummarizedExperiment object containing the DNA methylation data is loaded into the environment
load("../Shared_Data/LLS_betas.Rdata")

##############################################################
# Data preprocessing
# Store CRP
targets <- data.frame(
  uuid = colData(betas)$uuid,
  run_id = colData(betas)$run_id,
  crp = colData(betas)$hscrp)

# Remove values below LOD
targets <- targets %>% filter(crp>0)
  mutate(
    crp = ifelse(targets$crp == 0, NA, crp)
  )

# Check Normality
ggdensity(targets$crp)
ggqqplot(targets$crp)
shapiro.test(targets$crp)

# Natural Log Transformation
targets$log_crp <- log(targets$crp)

# 3IQR outlier Removal
lower_limit <- quantile(targets$crp, na.rm=TRUE)[[2]] - 
  (3 * IQR(targets$crp, na.rm=TRUE))
upper_limit <- quantile(targets$crp, na.rm=TRUE)[[4]] +
  (3 * IQR(targets$crp, na.rm=TRUE))
targets <- targets %>% 
  mutate(
    log_crp = ifelse(
      targets$crp < lower_limit | targets$crp > upper_limit, 
      NA, log_crp),
    crp = ifelse(
      targets$crp < lower_limit | targets$crp > upper_limit, 
      NA, crp)
  )

# Check Normality
ggdensity(targets$log_crp)
ggqqplot(targets$log_crp)
shapiro.test(targets$log_crp)

##############################################################
# Covariates
file_name <- "LLS_IOP1-450kDNAmethylation-AdiponectinLeptinIl6_20201009.sav"
lls_cytokines <- read_sav(paste0("../Shared_Data/", file_name))

# Sex coded as 0=male, 1=female
targets <- targets %>% 
  mutate(
    sex = factor(
      case_when(
        colData(betas)$sex == "male" ~ 0,
        colData(betas)$sex == "female" ~ 1
      ), labels = c("Male", "Female")))

# Smoking Status coded as 0=never, 1=ex, 2=current
targets <- targets %>% 
  mutate(
    smoking = factor(
      case_when(
        colData(betas)$smoking == "non-smoker" ~ 0,
        colData(betas)$smoking == "former-smoker" ~ 1,
        colData(betas)$smoking == "current smoker" ~ 2
      ), labels = c("Never", "Former", "Current")))

# IL-6 data
lls_cytokines <- lls_cytokines %>% filter(!is.na(IL6) & IL6>0)

# Cell counts
counts <- read_tsv("../Shared_Data/LLS_countsUCB.tsv")
colnames(counts) <- c("run_id", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")

# Join
targets <- left_join(targets, counts, by = "run_id")

# Technical covariates
targets$plate <- colData(betas)$sample_plate
targets$row <- colData(betas)$sentrix_position

# Row as numeric
targets <- targets %>% mutate(row = as.numeric(substr(row, 3, 3)))

# Join
targets <- left_join(targets, lls_cytokines, by = "uuid")
targets <- zap_labels(targets)

##############################################################
# Exclusion Criteria
# Sex Chromosomes
betas <- betas[!(seqnames(betas) %in% c("chrX", "chrY")),]

# Encode Blacklist Regions
load("../Shared_Data/ENCODE_Blacklist/ENCODEBlacklist_CpGomit-450K.RData")
betas <- betas[!rownames(betas) %in% cpg_blacklist,]

# Zhou probes
maskProbes <- read_tsv("../Shared_Data/Manifests/HM450.hg19.manifest.tsv.gz")
maskProbes <- maskProbes[maskProbes$MASK_general == TRUE,]$probeID
betas <- betas[!rownames(betas) %in% maskProbes,]

# 3IQR Outlier Removal
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

# Merging
targets <- targets %>% filter(uuid %in% rownames(colData(betas)))

# Check order
rownames(targets) <- targets$uuid
identical(rownames(targets), rownames(colData(betas)))

# Add phentype data to SummarizedExperiment object
colData(betas) <- DataFrame(targets)

##############################################################
# Analysis data preparation
# Base Model
betas <- subset(betas, select = rowSums(is.na(colData(betas)[, c(1:13)])) == 0)

# IL-6 Sensitivity Analysis
betas_il6 <- subset(betas, select = rowSums(is.na(colData(betas)[, c(1:13,15)])) == 0)

##############################################################
# EWAS
# Base Model
# Fit model
metadata(betas)$formula <- ~log_crp + age + sex + CD8T + CD4T + NK + Bcell + Mono + plate + row
design <- model.matrix(metadata(betas)$formula, data=colData(betas))
fit <- lmFit(assay(betas), design)

# Store results
beta <- fit$coefficients[, 2]
SE <- fit$stdev.unscaled[, 2] * fit$sigma
t <- beta / SE
p <- 2 * pt(-abs(t), fit$df.residual)
N <- ncol(design) + fit$df.residual

limma_base <- data.frame(cpg = rownames(fit$coefficients), 
                         beta, SE, p, N)

# Save results
save(limma_base, file='../IL6_Data/Cohort_Results/LLS/LLS_ln-crp_base.Rdata')

# IL-6 Sensitivity Analysis
# Fit model  
metadata(betas_il6)$formula <- ~log_crp + age + sex + il6 + CD8T + CD4T + NK + Bcell + Mono + plate + row
design <- model.matrix(metadata(betas_il6)$formula, data=colData(betas_il6))
fit <- lmFit(assay(betas_il6), design)

# Store results
beta <- fit$coefficients[, 2]
SE <- fit$stdev.unscaled[, 2] * fit$sigma
t <- beta / SE
p <- 2 * pt(-abs(t), fit$df.residual)
N <- ncol(design) + fit$df.residual

limma_base <- data.frame(cpg = rownames(fit$coefficients), 
                         beta, SE, p, N)

# Save results
save(limma_base, file='../IL6_Data/Cohort_Results/LLS/LLS_ln-crp_ext2.Rdata')

##############################################################










