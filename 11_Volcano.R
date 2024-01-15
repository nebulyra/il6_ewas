##############################################################
# SCRIPT 11: Volcano Plot
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(ggrepel)
library(GenomicRanges)
library(ggpubr)
library(DNAmArray)

# Load data
load("../IL6_Data/Results/IL6_EWAS_Results-Full.Rdata")
load("../IL6_Data/Results/IL6_EWAS_Results-Top.Rdata")

# Filter out significant CpGs that were heterogenous or failed sensitivity analysis
cpg_list <- il6_top$cpg
length(cpg_list)

il6_results <- il6_results %>% 
  filter(cpg %in% cpg_list | base_meta_padj_fdr > 0.05)
dim(il6_results)

##############################################################
# Figure 1c: Volcano plot
# Initialize
sig_limit <- max(il6_top$base_meta_p)

min <- as.numeric(-max(
  abs(il6_results$base_meta_beta),
  na.rm=TRUE)) - 0.005

min <- as.numeric(max(
  abs(il6_results$base_meta_beta),
  na.rm=TRUE)) + 0.005

p_max <- as.numeric(-log10(min(il6_results$base_meta_p, na.rm=TRUE))) + 2

# Plot
il6_results %>% 
  ggplot(aes(
    x=base_meta_beta,
    y=-log10(base_meta_p)
  )) +
  geom_hline(
    yintercept=-log10(sig_limit),
    linetype="dashed",
    size=0.5,
    color="#ACAFAF"
  ) +
  xlim(min, max) +
  geom_point(
    color=ifelse(il6_results$base_meta_p<=sig_limit, 
                   "#237194", "#c1c0bc"),
    size=0.8, 
  ) +
  ggtitle("") +
  ylab(bquote(-log[10]~"p")) +
  xlab("Effect Size") +
  theme(
    axis.text=element_text(
      size=9, 
      color="#1B2021"),
    axis.title=element_text(
      size=11, 
      hjust=0.5, 
      color="#1B2021"),
    plot.title=element_text(
      size=16, 
      hjust=0.5,
      face="bold", 
      color="#548687"),
    panel.background=element_rect(
      fill="white"),
    panel.border=element_rect(
      color="#1B2021",
      fill=NA),
    panel.grid.major=element_line(
      color="grey95"),
    panel.grid.minor=element_line(
      color="grey95"),
    plot.background=element_rect(
      fill="white"))

##############################################################