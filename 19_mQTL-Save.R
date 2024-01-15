##############################################################
# SCRIPT 19: Save mQTL data
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)

# Load data
load('../IL6_Data/Results/IL6_EWAS_Results-Top.Rdata')

##############################################################
# Setup
# List files
files <- list.files(path="../IL6_Data/mQTLs/out/", pattern='.txt')
length(files)
head(files)

# Load genome
genome <- BSgenome.Hsapiens.UCSC.hg19
all_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh37

##############################################################
# Save for each CpG
for(i in files){
  # Read in mQTL data
  mqtl_dat <- read_tsv(file=paste0('../IL6_Data/mQTL/',i),
                       col_names = c('MarkerName', 
                                   'Allele1', 'Allele2',
                                   'Freq1', 'FreqSE',
                                   'Effect', 'StdErr', 
                                   'Pvalue', 'Direction',
                                   'HetISq', 'HetChiSq',
                                   'HetDf', 'HetPVal',
                                   'EffectARE', 'StdErrARE',
                                   'PvalueARE', 'tausq', 
                                   'StdErrMRE', 'PvalueMRE',
                                   'TotalSampleSize'))
  # Filer out insertions and deletions
  mqtl_dat <- mqtl_dat %>% 
      dplyr::filter(Allele1 %in% c('a', 'c', 'g', 't') &
           Allele2 %in% c('a', 'c', 'g', 't'))

  # Format
  mqtl_dat <- mqtl_dat %>% 
      separate(
        MarkerName, into = c('CHR', 'POS', 'ProbeID'), sep=':') %>% 
      mutate(
        CHR = substr(CHR, 4, 5),
        POS = as.numeric(POS),
        ProbeID = substr(ProbeID, 5, 15)
      ) %>% 
      dplyr::select(
        ProbeID,
        BETA = Effect,
        SE = StdErr,
        N = TotalSampleSize,
        PVAL = Pvalue,
        MAF = Freq1,
        A1 = Allele1,
        A2 = Allele2,
        CHR, 
        POS
      ) %>% 
      mutate(
        BETA = as.numeric(BETA),
        SE = as.numeric(SE),
        N = as.numeric(N),
        PVAL = as.numeric(PVAL),
        MAF = 1-as.numeric(MAF),
        A1 = toupper(A1),
        A2 = toupper(A2),
        CHR = as.character(CHR),
        POS = as.numeric(POS),
        CHRPOS = paste0(CHR, ':', POS)
      ) %>% 
      arrange(CHR, POS) %>% unique()

  # Check
  print(paste0("There are ", nrow(mqtl_dat), " mQTLs for ", i))

  # GRanges
  print(paste0("Creating mQTL GRanges... "))
  positions <- GPos(seqnames=mqtl_dat$CHR, 
                    pos=mqtl_dat$POS)
  positions
  
  # Find overlaps between genome and mQTLs
  print("Finding overlaps between hg19 SNPs and mQTLs...")
  my_snps <- snpsByOverlaps(all_snps, 
                            positions, 
                            genome=genome)
  my_snps

  # Save alternate alleles
  alt <- unlist(lapply(my_snps$alt_alleles, function(x) {
        x[1]
      }))
  
  # Save as data frame
  my_snps <- data.frame(CHR = seqnames(my_snps),
                        POS = start(my_snps), 
                        SNP = my_snps$RefSNP_id,
                        A1 = my_snps$ref_allele,
                        A2 = alt)



  my_snps <- my_snps %>% mutate(
        CHR = as.numeric(CHR),
        POS = as.numeric(POS),
        CHRPOS = paste0(CHR,':',POS)) %>% 
        arrange(CHR, POS) %>% 
        distinct(CHRPOS, .keep_all=TRUE)
  my_snps <- my_snps %>% dplyr::select(SNP, CHRPOS)
  
  # Merge
  mqtl_dat <- left_join(mqtl_dat, my_snps, by = 'CHRPOS')

  mqtl_dat <- mqtl_dat %>%
       dplyr::select(ProbeID, SNP, BETA, SE, N, 
                     PVAL, MAF, A1, A2, CHR, POS)

  # Save
  save(mqtl_dat, 
       file=paste0('../IL6_Data/mQTLs/out',substr(i,1,10),'.Rdata'))
}

##############################################################

