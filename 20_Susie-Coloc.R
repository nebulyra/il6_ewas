##############################################################
# SCRIPT 20: SuSiE colocalization
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(SummarizedExperiment)
library(GenomicRanges)
library(coloc)
library(susieR)
library(LDlinkR)
library(Biostrings)
library(XVector)
library(edgeR)

# Load data
# eQTLgen
eqtl_dat <- read_tsv('../IL6_Data/eQTL/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz')
eqtl_maf <- read_tsv('../IL6_Data/eQTL/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz')
# IL-6 GWAS results (gwas_dat)
load('../IL6_Data/GWAS/all_snps.Rdata')
load('../IL6_Data/GWAS/IL6_gwas_dat.Rdata')
# Key for hg18 to hg19 from IL-6 GWAS
hg19 <- read_tsv('../IL6_Data/GWAS/snp_postion_lift_B36_B37_key.txt.gz')
# IL-6 EWAS results
load('../IL6_Data/Results/IL6_EWAS_Results-Top.Rdata')
# BIOS DNAm - betas
load('../../BIOS/eQTM/data/methData_Betas_BIOS_Freeze2_unrelated.RData')
# BIOS RNAseq - counts
load('../../BIOS/eQTM/data/rnaSeqData_ReadCounts_BIOS_Freeze2_unrelated_GRCh38.RData')
# Gene CpG pairs (res)
load('../IL6_Data/Susie/gene_cpg_pairs.Rdata')

##############################################################
# IL-6 GWAS data
# Keep hg18
gwas_dat <- gwas_dat %>% dplyr::select(-POS)
hg19 <- hg19 %>% dplyr::select(SNP, POS=BP_37)
gwas_dat <- left_join(gwas_dat, hg19, by='SNP')
head(gwas_dat)

##############################################################
# IL-6 CpGs
cpg_list <- il6_top$cpg

##############################################################
# eQTL data
# MAF
eqtl_maf <- eqtl_maf %>% 
  dplyr::select(SNP, 
                MAF = AlleleB_all)

# Create df
eqtl_dat <- eqtl_dat %>% 
  dplyr::select(ProbeID = GeneSymbol,
                SNP,
                CHR = SNPChr,
                POS = SNPPos,
                Zscore,
                N = NrSamples,
                PVAL = Pvalue,
                A1 = AssessedAllele,
                A2 = OtherAllele) %>% 
  dplyr::mutate(
    sigma_sqrd = 1/(1+(Zscore^2/N)),
    BETA = Zscore*sqrt(sigma_sqrd/N),
    SE = sqrt(sigma_sqrd/N)
  ) %>% 
  dplyr::select(
    ProbeID, SNP, CHR, POS, BETA, SE, N, PVAL, A1, A2)

# Merge
eqtl_dat <- left_join(eqtl_dat, eqtl_maf, by='SNP')

##############################################################
# BIOS data processing - RNAseq
# Drop sex chromosomes
counts  <- dropSeqlevels(counts,  c("X", "Y"), pruning.mode = "coarse")

# Remove missing flowcell data
count.coldata <- colData(counts)
idx <- which(is.na(count.coldata$flowcell_num) == TRUE)
count.coldata <- count.coldata[-idx,]
counts <- counts[,-idx]

# Remove lowly expressed genes
counts <- counts[rowSums(assays(counts)$data > 0) > 0.5 * ncol(counts), ]

# Log2 CPM
log.cpm <- DGEList(counts = assays(counts)$data)
log.cpm <- calcNormFactors(log.cpm)
log.cpm <- cpm(log.cpm, log = T)

# RIN transformation
RIN <- function(x) {
  y <- rank(x, NA)
  y <- ppoints(y)[y]
  y <- qnorm(y)
  x[!is.na(x)] <- y
  x
}
RIN.counts <- t(apply(log.cpm, 1, RIN))

##############################################################
# Gene names
ens2gene <- cinaR::grch38
ens2gene <- ens2gene %>% dplyr::select(gene_ens = ensgene,
                                       gene_sym = symbol)

# Merge with gene-Cpg pair data (res)
res <- unique(left_join(res, ens2gene, by='gene_ens'))

##############################################################
# Colocalization analysis
# Save previous
prev_coloc <- list.files(path='../IL6_Data/Coloc_out/')

# Print
print(paste0('We are looking at ', nrow(res), 
             ' CpGs-gene pairs...'))

for(i in 1:nrow(res)){
  # Save CpG and gene pair names
  cpg <- res$cpg[i]
  gene <- res$gene_sym[i]
  
  # Save output file name
  coloc_filename <- paste0(
    'expression_', cpg, '_', gene, '_ext.Rdata'
  )
  
  # If not already run
  if(coloc_filename %in% prev_coloc){
    print('Coloc already run... skipping...')
  } else {
    # Save ensembl ID
    gene_ens <- res$gene_ens[i]
    
    # Print progress
    print(paste0('Looking at CpG ', cpg, 
                 " and gene ", gene,
                 ', which are ', i, 
                 ' out of ', nrow(eqtm),
                 " CpG-gene pairs..."))
    
    # Save mQTL file names
    files <- list.files(path='../IL6_Data/mQTLs/out')
    
    # Save mQTL file name for this CpG
    look <- paste0(cpg, '.Rdata')
    
    # If there is no mQTL data for this CpG in GoDMC
    if(!look %in% files){
      print(paste0('mQTL data not available for ',cpg))
    } else {
      # Store output manhattan plots
      png(paste0('../IL6_Data/Susie/Plots/expression_',cpg,'_',gene,'.png'))
      
      # Load mQTL data
      load(paste0('../IL6_Data/mQTLs/out/',cpg,'.Rdata'))
      
      # If there is no mQTL data
      if(nrow(mqtl_dat)==0){
        print(paste0('mQTL data not available for ',cpg))
      } else {
        # Get CpG SD from BIOS
        betas_df <- as.data.frame(t(assay(betas))) %>% 
            dplyr::select(all_of(cpg))
        cpg_sd <- sd(betas_df[,1], na.rm=TRUE)
        print(paste0('SD of CpG ', cpg, ' is ', cpg_sd))
          
        # Get gene SD from BIOS
        print(paste0('Gene is', gene,
                       ' (', gene_ens, ')'))
          
        rin <- as.data.frame(t(RIN.counts))
        rin_sub <- rin %>% dplyr::select(all_of(gene_ens))
        gene_sd <- sd(rin_sub[,1], na.rm=TRUE)
        print(paste0('SD of gene ', gene, ' is ', gene_sd))
          
        # Make GRange object for +-1000kb from each IL-6 CpG
        bed <- il6_top %>% 
          filter(cpg == res$cpg[i]) %>% 
          select(ProbeID = cpg,
                 CHR = cpg_chr_hg19,
                 POS = cpg_start_hg19) %>% 
          mutate(CHR = sub('chr', '', CHR),
                 START = POS - 1000000,
                 STOP = POS + 1000000) %>% 
          select(-POS)
        print(bed)
        
        # Save mQTLs within 1000kb of CpG
        if(length(bed$STOP) == 1){
          mqtl <- mqtl_dat %>% 
            dplyr::filter(ProbeID == res$cpg[i] &
                          CHR == bed$CHR &
                          POS >= bed$START &
                          POS <= bed$STOP) %>% 
            arrange(CHR, POS)
          print(dim(mqtl))
            
        # Save eQTLs withint 1000kb of CpG for this gene
        eqtl <- eqtl_dat %>% 
          dplyr::filter(ProbeID == gene &
                        CHR == bed$CHR &
                        POS >= bed$START &
                        POS <= bed$STOP) %>% 
          arrange(CHR, POS)
          print(dim(eqtl))
        
        # Make into list    
        data_genome <- list(eqtl, mqtl)
            
        # Subset only SNPs with mQTL and eQTL data
        eqtl <- eqtl %>% 
          filter(SNP %in% mqtl$SNP) %>% 
          arrange(SNP)
        mqtl <- mqtl %>% 
          filter(SNP %in% eqtl$SNP) %>% 
          arrange(SNP)
            
        # Get LD if more than 1 SNP
        if(length(eqtl$SNP)>=2){
          # LD filename
          look_long <- paste0(cpg,'_',gene,'.Rdata')
          
          # Previously saved LD files
          LD_files <- list.files(path='../IL6_Data/LD/')
          
          # If LD already saved, load it    
          if(look_long %in% LD_files){
            load(paste0('../IL6_Data/LD/',cpg,'_',gene,'.Rdata'))
          # If LD not already saved, fetch it
          } else {
            LD <- LDmatrix(eqtl$SNP, token='7c9554d72b01')
            names <- LD$RS_number
          
          # If LD matrix was available      
          if(nrow(LD) > 1){
            # Format
            LD <- as.matrix(LD[,2:ncol(LD)])
            row.names(LD) <- names
            
            # R2      
            LD <- apply(LD, 2, sqrt)
            
            # Save
            save(LD, 
                 file=paste0('../IL6_Data/LD/',cpg,'_',gene,'.Rdata'))
          }
          }
          
          # Remove NAs from LS
          logi <- rowSums(is.na(LD)) == 0
          LD <- LD[logi,]
          
          if(nrow(LD) > 1){
            # Output unavailable SNPs
            print(eqtl %>% filter(!SNP %in% colnames(LD)))
            print(mqtl %>% filter(!SNP %in% colnames(LD)))
            
            # Subset and calculate variances
            eqtl <- eqtl %>% 
              filter(SNP %in% colnames(LD)) %>% 
              mutate(var = SE * SE) %>% 
              arrange(POS)
            eqtl_var <- (0.15/gene_sd)^2
            print(paste0('Variance for eQTLs: ', eqtl$var))
            
            mqtl <- mqtl %>% 
              filter(SNP %in% colnames(LD)) %>% 
              mutate(var = SE * SE) %>% 
              arrange(POS)
            mqtl_var <- (0.15/cpg_sd)^2
            print(paste0('Variance for mQTLs: ', mqtl$var))
                  
            # SNP match check
            mqtl$SNP == eqtl$SNP
            mqtl$A1 == eqtl$A1
            mqtl$MAF <- eqtl$MAF
                  
            # Make lists
            eqtl <- list(
                    pvalues = eqtl$PVAL,
                    beta = eqtl$BETA,
                    varbeta = eqtl$var,
                    N = max(eqtl$N),
                    MAF = eqtl$MAF,
                    LD = LD,
                    type = "quant",
                    sdY=gene_sd,
                    snp = eqtl$SNP,
                    position = eqtl$POS
                  )
                  
              mqtl <- list(
                    pvalues = mqtl$PVAL,
                    beta = mqtl$BETA,
                    varbeta = mqtl$var,
                    N = max(mqtl$N),
                    MAF = mqtl$MAF,
                    LD = LD,
                    type = "quant",
                    sdY=cpg_sd,
                    snp = mqtl$SNP,
                    position = mqtl$POS
                  )
               
              # SuSiE check datasets
              check_dataset(eqtl, req="LD")
              check_dataset(mqtl, req='LD')
                    
              # Run SuSiE - eQTL data
              S3 = runsusie(
                eqtl,
                prior_variance=eqtl_var,
                estimate_prior_variance=FALSE)
              print(summary(S3))
              
              # Run SuSiE - mQTL data
              S4 = runsusie(
                mqtl,
                prior_variance = mqtl_var,
                estimate_prior_variance=FALSE)
              print(summary(S4))
              
              # Colocalization
              susie.res=coloc.susie(S3,S4)
              print(susie.res$summary)
              
              # Save
              save(susie.res, 
                   file=paste0('../IL6_Data/Susie/Coloc_out/expression_', 
                               cpg, '_', gene, '_ext.Rdata'))
              
              # Check if H4 > 0.9
              sensitivity(susie.res,
                          "H4 > 0.9",
                          dataset1=eqtl,
                          dataset2=mqtl)
          }
        }
        }
      }
      dev.off()
    }
  }
}

##############################################################


