##############################################################
# SCRIPT 3: Loading and saving cohort-level results
##############################################################
# Setup
rm(list=ls())
library(tidyverse)

##############################################################
# IL6 EWAS - Base Model
# LLS
load("../IL6_Data/Cohort_Results/LLS/LLS_ln-il6_base.Rdata")

limma_base <- limma_base %>% dplyr::select(cpg, 
    cpg_mean_lls = cpg_mean, cpg_se_lls = cpg_se,
    beta_lls = beta, SE_lls = SE, p_lls = p, N_lls = N)

# KORA
df <- read_csv2("../IL6_Data/Cohort_Results/KORA/KORA_ln-il6_base.csv")

df <- df %>% dplyr::select(cpg = cg, 
    cpg_mean_kora = mean_cpg, cpg_se_kora = sd_cpg,
    beta_kora = beta, SE_kora = se, p_kora = p, N_kora = N_obs)

limma_base <- full_join(limma_base, df, by = "cpg")

# NTR
df <- read_csv("../IL6_Data/Cohort_Results/NTR/NTR_ln-il6_base.csv")
stats <- read_csv("../IL6_Data/Cohort_Results/NTR/NTR_CpG-meanSD.csv")
df <- left_join(df, stats, by=c("CpG" = "cpg")) 

df <- df %>% dplyr::select(cpg = CpG, 
    cpg_mean_ntr = mean, cpg_sd_ntr = sd,
    beta_ntr = beta, SE_ntr = SE, p_ntr = p, N_ntr = N)

limma_base <- full_join(limma_base, df, by = "cpg")

# Save 
save(limma_base, file="../IL6_Data/Cohort_Results/ALL_ln-il6_base.Rdata")

##############################################################
# IL6 EWAS - Smoking Sensitivity Analysis
# LLS
load("../IL6_Data/Cohort_Results/LLS/LLS_ln-il6_ext1.Rdata")
limma_base <- limma_base %>% dplyr::select(cpg, 
    beta_lls = beta, SE_lls = SE, p_lls = p, N_lls = N)

# KORA
df <- read_csv2("../IL6_Data/Cohort_Results/KORA/KORA_ln-il6_ext1.csv")
df <- df %>% dplyr::select(cpg = cg, 
    beta_kora = beta, SE_kora = se, p_kora = p, N_kora = N_obs)

limma_base <- full_join(limma_base, df, by = "cpg")

# NTR
df <- read_csv("../IL6_Data/Cohort_Results/NTR/NTR_ln-il6_ext1.csv")

df <- df %>% dplyr::select(cpg = CpG, 
    beta_ntr = beta, SE_ntr = SE, p_ntr = p, N_ntr = N)

limma_base <- full_join(limma_base, df, by = "cpg")

# Save 
save(limma_base, file="../IL6_Data/Cohort_Results/ALL_ln-il6_ext1.Rdata")

##############################################################
# IL6 EWAS - CRP Sensitivity Analysis
# LLS
load("../IL6_Data/Cohort_Results/LLS/LLS_ln-il6_ext2.Rdata")

limma_base <- limma_base %>% dplyr::select(cpg, 
    beta_lls = beta, SE_lls = SE, p_lls = p, N_lls = N)

# KORA
df <- read_csv2("../IL6_Data/Cohort_Results/KORA/KORA_ln-il6_ext2.csv")

df <- df %>% dplyr::select(cpg = cg, 
    beta_kora = beta, SE_kora = se, p_kora = p, N_kora = N_obs)

limma_base <- full_join(limma_base, df, by = "cpg")

# NTR
df <- read_csv("../IL6_Data/Cohort_Results/NTR/NTR_ln-il6_ext2.csv")

df <- df %>% dplyr::select(cpg = CpG, 
    beta_ntr = beta, SE_ntr = SE, p_ntr = p, N_ntr = N)

limma_base <- full_join(limma_base, df, by = "cpg")

# Save
save(limma_base, file="../IL6_Data/Cohort_Results/ALL_ln-il6_ext2.Rdata")

##############################################################
# IL6 EWAS - IDOL Extended Sensitivity Analysis
# LLS
load("../IL6_Data/Cohort_Results/LLS/LLS_ln-il6_ext3.Rdata")

limma_base <- limma_base %>% dplyr::select(cpg, 
    beta_lls = beta, SE_lls = SE, p_lls = p, N_lls = N)

# KORA
df <- read_csv2("../IL6_Data/Cohort_Results/KORA/KORA_ln-il6_ext3.csv")

df <- df %>% dplyr::select(cpg = cg, 
    beta_kora = beta, SE_kora = se, p_kora = p, N_kora = N_obs)

limma_base <- full_join(limma_base, df, by = "cpg")

# NTR
df <- read_csv("../IL6_Data/Cohort_Results/NTR/NTR_ln-il6_ext3.csv")

df <- df %>% dplyr::select(cpg = CpG, 
    beta_ntr = beta, SE_ntr = SE, p_ntr = p, N_ntr = N)

limma_base <- full_join(limma_base, df, by = "cpg")

# Save
save(limma_base, file="../IL6_Data/Cohort_Results/ALL_ln-il6_ext3.Rdata")

##############################################################
# CRP EWAS - Base Model
# LLS

load("../IL6_Data/Cohort_Results/LLS/LLS_ln-crp_base.Rdata")

limma_base <- limma_base %>% dplyr::select(cpg, 
    beta_lls = beta, SE_lls = SE, p_lls = p, N_lls = N)

# KORA
df <- read_csv2("../IL6_Data/Cohort_Results/KORA/KORA_ln-crp_base.csv")

df <- df %>%dplyr::select(cpg = cg, 
    beta_kora = beta, SE_kora = se, p_kora = p, N_kora = N_obs)

limma_base <- full_join(limma_base, df, by = "cpg")

# NTR
df <- read_csv("../IL6_Data/Cohort_Results/NTR/NTR_ln-crp_base.csv")

df <- df %>% dplyr::select(cpg = CpG, 
    beta_ntr = beta, SE_ntr = SE, p_ntr = p, N_ntr = N)

limma_base <- full_join(limma_base, df, by = "cpg")

# Save 
save(limma_base, file="../IL6_Data/Cohort_Results/ALL_ln-crp_base.Rdata")

##############################################################
# CRP EWAS - IL6 Sensitivity Analysis
# Leiden Longevity Study (LLS) 
load("../IL6_Data/Cohort_Results/LLS/LLS_ln-crp_ext2.Rdata")

limma_base <- limma_base %>% dplyr::select(cpg, 
    beta_lls = beta, SE_lls = SE, p_lls = p, N_lls = N)

# KORA
df <- read_csv2("../IL6_Data/Cohort_Results/KORA/KORA_ln-crp_ext2.csv")

df <- df %>% dplyr::select(cpg = cg, 
    beta_kora = beta, SE_kora = se, p_kora = p, N_kora = N_obs)

limma_base <- full_join(limma_base, df, by = "cpg")

# NTR
df <- read_csv("../IL6_Data/Cohort_Results/NTR/NTR_ln-crp_ext2.csv")

df <- df %>% dplyr::select(cpg = CpG, 
    beta_ntr = beta, SE_ntr = SE, p_ntr = p, N_ntr = N)

limma_base <- full_join(limma_base, df, by = "cpg")

# Save 
save(limma_base, file="../IL6_Data/Cohort_Results/ALL_ln-crp_ext2.Rdata")

##############################################################


