##############################################################
# SCRIPT 21: SuSiE results
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(coloc)

##############################################################
# Save coloc out filenames
files <- list.files('../IL6_Data/Susie/Coloc_out/',
                    pattern = ".Rdata")

# For each CpG gene pair
for(i in files){
  # Save cpg and gene names
  cpg <- strsplit(i, "_")[[1]][2]
  gene <- strsplit(i, "_")[[1]][3]
  
  # Store list of CpGs tested
  if(i == files[1]){
    cpg_list <- cpg
  } else {
    cpg_list <- c(cpg_list, cpg)
  }
  
  # Load SuSiE results
  load(i)
  
  # Test H4 > 0.9
  res <- sensitivity(susie.res, "H4 > 0.9", plot.manhattans = FALSE)
  pass <- all((res %>% filter(p12 > 2.6E-06))$pass == TRUE)
  row <- (res %>% filter(p12 > 2.6E-06))[1,]
  
  # Save results  
  df <- data.frame(cpg = cpg,
                   gene = gene,
                   pass = pass,
                   PP.H0.abf = row$PP.H0.abf,
                   PP.H1.abf = row$PP.H1.abf,
                   PP.H2.abf = row$PP.H2.abf,
                   PP.H3.abf = row$PP.H3.abf,
                   PP.H4.abf = row$PP.H4.abf)
  
  # Merge
  if(i == files[1]){
    out <- df
  } else {
    out <- rbind(out, df)
  }
}

# Save results
write_csv(out, file='../IL6_Data/Tables/ST8.csv')

##############################################################