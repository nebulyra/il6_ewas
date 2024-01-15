##############################################################
# SCRIPT 10: Identification of distinct loci
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(ggrepel)
library(GenomicRanges)
library(ggpubr)
library(DNAmArray)

# Load data
load("../IL6_Data/Results/IL6_EWAS_Results-Top.Rdata")

##############################################################
# Data processing
il6_top <- il6_top %>% 
  dplyr::select(chromosome = cpg_chr_hg19,
    start = cpg_start_hg19, padj = base_meta_padj_fdr) %>% 
  filter(!is.na(chromosome) & chromosome!='chrX')

# 5% significance identifier
il6_top <- il6_top %>% mutate(
  padj = ifelse(is.na(padj), 1, padj))

il6_top <- il6_top %>% mutate(
  crit = ifelse(padj <= 0.05, 1, 0))

# Remove NAs
il6_top <- il6_top %>% filter(!is.na(crit) &
                                !is.na(chromosome))

# Arrange by position
il6_top <- il6_top %>% arrange(chromosome, start) %>% 
  dplyr::select(chromosome, start, crit)

##############################################################
# DMRfinder
# Initialize
chromosome=1:22
MAXIMUM_REGION_LENGTH = 1000
mismatches = 3
chr_list <- 1:22

# Find distinct loci
for(x in chr_list){
  
  chr1 = il6_top[il6_top[,1]==x,]
  chr1 <- chr1 %>% arrange(start)
  chr.final = data.frame(
    coord = chr1$start,
    crit = chr1$crit
  )
  
  last_coordinate = length( chr.final$crit )
  next_coordinate = 0
  
  for (i in 1:(last_coordinate-1)) {
    if ( i>=next_coordinate ) {
      if (chr.final$crit[ i ]==1) {
        start_location = chr.final$coord[ i ]
        last_visited_crit_loc = start_location
        sum_of_ones = 1
        number_of_items = 1
        
        # start crawling loop
        for (j in (i+1):last_coordinate ) {
          if (chr.final$coord[ j ] > (last_visited_crit_loc + MAXIMUM_REGION_LENGTH)) { break }
          if((number_of_items-sum_of_ones)>mismatches) { break }   #Number of mismatches
          number_of_items = number_of_items + 1
          if (chr.final$crit[j]==1) { 
            last_visited_crit_loc = chr.final$coord[ j ]
            sum_of_ones = sum_of_ones + 1 
          }
        }
        
        # now check the result
        if (sum_of_ones>=3) {
          last_one=i+number_of_items-1
          for (k in (i+number_of_items-1):1) {
            if ( chr.final$crit[k] == 0 ) {
              last_one = last_one - 1
              number_of_items = number_of_items - 1
            }
            else {
              break
            }
          }
          cat(x, ';',start_location,";",chr.final$coord[last_one],";",sum_of_ones/number_of_items,"\n")
          next_coordinate = last_one + 1
        }
      }
    }
  }
}
