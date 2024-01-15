##############################################################
# SCRIPT 16: Identify Genes
##############################################################
# Setup
rm(list=ls())
library(edgeR)
library(cinaR)
library(limma)
library(tidyverse)
library(GenomicFeatures)
library(AnnotationHub)
library(SummarizedExperiment)

# Load data
load('../IL6_Data/Results/IL6_EWAS_Results-Top.Rdata')

##############################################################
# Data on gene locations
# Make txdb
txdb <- makeTxDbFromEnsembl(organism="Homo sapiens",
                            release=75,
                            circ_seqs=NULL,
                            server="ensembldb.ensembl.org",
                            username="anonymous", 
                            password=NULL, 
                            port=0L,
                            tx_attrib=NULL)

# GRanges object
gene_range <- unlist(cdsBy(txdb, by = "gene"))
gene_range <- unique(gene_range)
gene_range

# Save GRanges
save(gene_range, file='../Shared_Data/geneRanges-75.Rdata')

##############################################################
# Data on CpG locations
# Make GRanges
cpg_ranges <- GRanges(seqnames=sub('chr', '', 
                                   il6_top$cpg_chr_hg19), 
                      IRanges(start=il6_top$cpg_start_hg19,
                              width=2,
                              names=il6_top$cpg))
cpg_ranges

##############################################################
# Genes within 100kb
# Initialize
distance.matrix <- matrix(NaN, 
                          nrow = length(gene_range), 
                          ncol = length(cpg_ranges), 
                          dimnames = list(names(gene_range), 
                                          names(cpg_ranges)))

# Store distances
for(i in 1:nrow(il6_top)){
  cpg <- cpg_ranges[i]
  calc.distance <- distance(cpg, gene_range, ignore.strand=T)
  distance.matrix[,i] <- calc.distance
}
dim(distance.matrix)

# Process
distance.matrix <- as.data.frame(distance.matrix)
colnames(distance.matrix) <- names(cpg_ranges)
distance.matrix$gene <- names(gene_range)

# List of genes
for(i in names(cpg_ranges)){
  df <- (distance.matrix %>% dplyr::select(gene, all_of(i)))
  df <- df[,2]<=100000
  out <- data.frame(
    gene=df$gene,
    cpg=i
  )
  
  if(i==names(cpg_ranges)[1]){
    res <- out
  } else {
    res <- rbind(res, out)
  }
}

# Save list
save(res, file='../IL6_Data/Susie/gene_cpg_pairs.Rdata')

##############################################################
  