##############################################################
# SCRIPT 22: Two-sample MR
##############################################################
# Setup
rm(list=ls())
library(TwoSampleMR)
library(tidyverse)
library(ieugwasr)
library(jsonlite)

##############################################################
# Read in data
# Read in coloc results and save colocalized CpGs
res <- read_csv(file='../IL6_Data/Tables/ST8.csv')
res <- res %>% filter(pass == TRUE)

# Look at investigated CpGs
length(unique(res$cpg))
cpg_list <- unique(res$cpg)

# IL-6 GWAS data
il6_dat <- read_tsv('../IL6_Data/GWAS/IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz')
hg19 <- read_tsv('../Data/snp_postion_lift_B36_B37_key.txt.gz')
hg19 <- hg19 %>% select(SNP, pos=BP_37)

# Format
il6_dat <- il6_dat %>% 
  select(
    SNP, chr=CHR,  
    beta.outcome = beta,
    se.outcome = se,
    samplesize.outcome = N,
    pval.outcome = p,
    eaf.outcome = EAF,
    effect_allele.outcome = EA,
    other_allele.outcome = OA) %>% 
  mutate(
    outcome = 'IL6',
    id.outcome = 'IL6')

# Merge
il6_dat <- left_join(il6_dat, hg19, by='SNP')

##############################################################
# DNAm --> IL-6
# Loop through CpGs
for(i in cpg_list){
  print(paste0('Looking at ', i))
  
  # Exposure data for CpGs (mQTLs) - downloaded from GoDMC 
  exposure_dat <- fromJSON(
      paste0('http://api.godmc.org.uk/v0.1/assoc_meta/cpg/',
             i))$assoc_meta
  
  # Save
  write_csv(exposure_dat,
            file=paste0('../IL6_Data/MR/cis_mQTLs/',i,'.csv'))
    
  # Format
  exposure_dat <- exposure_dat %>% 
    filter(a1 != 'I' & a1 != 'D') %>% #Remove indels
    filter(cistrans == TRUE) %>%  # Keep only cis-mQTLs
    select(SNP = rsid, 
           beta = beta_a1, se,
           effect_allele = a1, other_allele = a2, 
           Phenotype = cpg, 
           position = snp, samplesize, pval)
  
  print(paste0(i, ' has ', 
               nrow(exposure_dat), ' cis-mQTLs in GoDMC...'))
  
  if(nrow(exposure_dat)!=0){
    exposure_dat <- format_data(exposure_dat, type = 'exposure')
    exposure_dat <- exposure_dat %>% arrange(pval.exposure)
      
    # Keep only mQTLs in IL-6 data
    exposure_dat <- exposure_dat %>% 
      filter(SNP %in% il6_dat$SNP) 
    print(paste0(nrow(exposure_dat), 
                   ' of these are in the IL6 data...'))
      
    # Clumping
    if(nrow(exposure_dat)!=0){
      exposure_dat <- clump_data(exposure_dat)
      print(paste0(nrow(exposure_dat), 
                     ' of these are independent...'))
      print(exposure_dat$SNP)
        
      # Outcome data from GWAS
      # Keep only matching SNPs
      outcome_dat <- il6_dat %>% 
        filter(SNP %in% exposure_dat$SNP)
      print(paste0(nrow(outcome_dat), 
                     ' of these are in the IL6 summary data...'))
        
      print(exposure_dat)
      print(outcome_dat)
      
      # Harmonise data
      if(nrow(outcome_dat)!=0){
        dat <- harmonise_data(exposure_dat, outcome_dat)
          
        # Perform SMR
        res <- mr(dat)
        print(res)
          
        if(i == cpg_list[1]){
          res_full <- res
        } else {
          res_full <- rbind(res, res_full)
        }
      }
    }
  }
}

# Adjust p-values
res_full$padj <- p.adjust(res_full$pval, method='fdr')
res_full <- res_full %>% arrange(padj) %>% 
  select(exposure, b, se, pval, padj)
print(res_full)

# Save
save(res_full, file='../IL6_Data/DNAm_to_IL6.Rdata')
write_csv(res_full, files='../IL6_Data/Tables/ST11.csv')

##############################################################
# IL-6 --> DNAm
# Exposure data
exposure_dat <- read_csv('../IL-6 Data/GWAS/il6_exposure.csv')

# Save list of SNPs
snplist <- paste0('chr', exposure_dat$chr.exposure, ':',
                  exposure_dat$pos.exposure)
snplist

# Outcome data - trans-mQTLs from GoDMC
file_list <- list.files('../IL6_Data/MR/17/', pattern = '17')

# Scan function
mqtl_scan <- function(x) {
  print(x)
  mqtl_data <- read_tsv(paste0('../IL6_Data/MR/17/',x))
  mqtl_data <- separate(mqtl_data, MarkerName, c('snp', 'cpg'), '_')
  mqtl_data$snp <- sub("(:[^:]+):.*", "\\1", mqtl_data$snp)
  mqtl_data <- mqtl_data %>% filter(snp %in% snplist)
  save(mqtl_data, file=paste0('../IL6_Data/MR/trans_mQTLs/',x,'.Rdata'))
}

# Run scan
for(i in file_list){
  mqtl_scan(i)
}

# Clumping
exposure_dat <- clump_data(exposure_dat)
print(paste0('IL6 has ', nrow(exposure_dat), 
             ' independent instrumental SNPs'))

# Map SNPs
key <- exposure_dat %>% 
  dplyr::select(SNP, chr.exposure, pos.exposure) %>% 
  mutate(snp = paste0('chr', chr.exposure, ':', pos.exposure))
key <- key %>% dplyr::select(SNP, snp)

# Load trans mQTLs
file_list <- list.files('../Test/IL6_Base/MR/trans_mQTLs/')

# Format outcome data
outcome_dat <- data.frame()
outcome_dat <- lapply(file_list, function(x) {
  load(paste0('../IL6_Data/MR/trans_mQTLs/',x))
  outcome <- right_join(key, mqtl_data, by = 'snp')
  outcome <- outcome %>% 
    dplyr::select(SNP,
                  beta.outcome = Effect,
                  se.outcome = StdErr,
                  samplesize.outcome = TotalSampleSize,
                  pval.outcome = Pvalue,
                  eaf.outcome = Freq1, 
                  effect_allele.outcome = Allele1,
                  other_allele.outcome = Allele2,
                  outcome = cpg, 
                  id.outcome = cpg) %>% 
    mutate(effect_allele.outcome = toupper(effect_allele.outcome),
           other_allele.outcome = toupper(other_allele.outcome))
  return(outcome)
})
outcome_dat <- bind_rows(outcome_dat)

# Save only CpGs of interest
outcome_dat <- outcome_dat %>% filter(outcome %in% cpg_list)

# Save
save(outcome_dat, file = paste0('../IL6_Data/MR/Outcome/IL6.RData'))

# Harmonise
dat <- harmonise_data(
  exposure_dat, outcome_dat)

# SMR
res_ivw <- mr(dat, method = c('mr_ivw'))

# Adjust p-values
res_ivw$padj <- p.adjust(res_ivw$pval, method = 'fdr')
print(res_ivw %>% arrange(padj) %>% 
        select(id.outcome, b, se, pval, padj))

# Save
save(res_ivw, file='../IL6_Data/IL6_to_DNAm.Rdata')
write_csv(res_ivw, file='../IL6_Data/Tables/ST12.csv')

##############################################################




