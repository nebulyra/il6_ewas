##############################################################
# SCRIPT 5: Bacon - Correcting for bias and inflation in the
# test statistics
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(bacon)

##############################################################
# IL6 EWAS - Base Model
# Load Data
load('../IL6_Data/Processing/ALL-2-PRUNE_ln-il6_base.Rdata')
limma_base[1:5, 1:5]

# Bacon function
run_bacon <- function(cohort){
  df <- limma_base %>% 
    dplyr::select(cpg,
      beta = paste0("beta_", cohort),
      se = paste0("SE_", cohort),
      N = paste0("N_", cohort)
    ) %>% 
    filter(!is.na(beta) & !is.na(se))
  
  cpg <- df$cpg
  num <- df$N
  es <- as.numeric(df$beta)
  se <- as.numeric(df$se)
  
  bc <- bacon(NULL, es, se)
  print(bc)
  
  tstat <- tstat(bc)
  p <- pval(bc) 
  coef <- es(bc)
  se <- se(bc)
  
  bc <- bacon(NULL, coef, se)
  print(bc)
  
  assign(paste0("bacon_", cohort), data.frame(
    cpg, coef, se, p, tstat, num
  ), envir=.GlobalEnv)
  
  name_list <- c(
    "cpg", 
    paste0("beta_", cohort), 
    paste0("SE_", cohort),
    paste0("p_", cohort),
    paste0("tstat_", cohort),
    paste0("N_", cohort))
  
  assign(paste0("bacon_",cohort),
         setNames(get(paste0("bacon_",cohort)), name_list),
         envir=.GlobalEnv)
}

# Run for each cohort
run_bacon("lls")
run_bacon("kora")
run_bacon("ntr")

# Merge
bacon_df <- full_join(bacon_lls, bacon_kora, by = "cpg")
bacon_df <- full_join(bacon_df, bacon_ntr, by = "cpg")

# Save
save(bacon_df, file='../IL6_Data/Processing/ALL-3-BACON_ln-il6_base.Rdata')

##############################################################
# IL6 EWAS - Smoking Sensitivity Analysis
# Load data
load('../IL6_Data/Processing/ALL-2-PRUNE_ln-il6_ext1.Rdata')

# Run bacon
run_bacon("lls")
run_bacon("kora")
run_bacon("ntr")

# Merge
bacon_df <- full_join(bacon_lls, bacon_kora, by = "cpg")
bacon_df <- full_join(bacon_df, bacon_ntr, by = "cpg")

# Save
save(bacon_df, file='../IL6_Data/Processing/ALL-3-BACON_ln-il6_ext1.Rdata')

##############################################################
# IL6 EWAS - CRP Sensitivity Analysis
# Load data
load('../IL6_Data/Processing/ALL-2-PRUNE_ln-il6_ext2.Rdata')

# Bacon
run_bacon("lls")
run_bacon("kora")
run_bacon("ntr")

# Merge
bacon_df <- full_join(bacon_lls, bacon_kora, by = "cpg")
bacon_df <- full_join(bacon_df, bacon_ntr, by = "cpg")

# Save
save(bacon_df, file='../IL6_Data/Processing/ALL-3-BACON_ln-il6_ext2.Rdata')

##############################################################
# IL6 EWAS - IDOL Extended Cell Counts Sensitivity Analysis
# Load data
load('../IL6_Data/Processing/ALL-2-PRUNE_ln-il6_ext3.Rdata')

# Bacon
run_bacon("lls")
run_bacon("kora")
run_bacon("ntr")

# Merge
bacon_df <- full_join(bacon_lls, bacon_kora, by = "cpg")
bacon_df <- full_join(bacon_df, bacon_ntr, by = "cpg")

# Save
save(bacon_df, file='../IL6_Data/Processing/ALL-3-BACON_ln-il6_ext3.Rdata')

##############################################################
# CRP EWAS - Base Model
# Load Data
load('../IL6_Data/Processing/ALL-2-PRUNE_ln-crp_base.Rdata')

# Bacon
run_bacon("lls")
run_bacon("kora")
run_bacon("ntr")

# Merge
bacon_df <- full_join(bacon_lls, bacon_kora, by = "cpg")
bacon_df <- full_join(bacon_df, bacon_ntr, by = "cpg")

# Save
save(bacon_df, file='../IL6_Data/Processing/ALL-3-BACON_ln-crp_base.Rdata')

##############################################################
# CRP EWAS - IL6 Sensitivity Analysis
# Load Data
load('../IL6_Data/Processing/ALL-2-PRUNE_ln-crp_ext2.Rdata')

# Bacon
run_bacon("lls")
run_bacon("kora")
run_bacon("ntr")

# Merge
bacon_df <- full_join(bacon_lls, bacon_kora, by = "cpg")
bacon_df <- full_join(bacon_df, bacon_ntr, by = "cpg")

# Save
save(bacon_df, file='../IL6_Data/Processing/ALL-3-BACON_ln-crp_ext2.Rdata')

##############################################################