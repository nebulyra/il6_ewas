##############################################################
# SCRIPT 13: CRP Investigations
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(ggrepel)
library(GenomicRanges)
library(ggpubr)

# Load data
load("../IL6_Data/Results/IL6_EWAS_Results-Full.Rdata")
load("../IL6_Data/Results/IL6_EWAS_Results-Top.Rdata")

##############################################################
# Adjusting IL-6 effects for CRP
# Save effects
il6_top <- il6_top %>% 
  mutate(crp_effect=abs(base_meta_beta-ext2_meta_beta),
         crp_perc=abs(crp_effect*100/base_meta_beta))

# Look at large effects
il6_top %>% 
  arrange(desc(crp_perc)) %>% 
  select(cpg, crp_effect, crp_perc) %>% 
  head(20)

# Figure 1d
ggplot(il6_top, 
       aes(x=base_meta_beta, 
           y=ext2_meta_beta)) +
  geom_hline(aes(yintercept=0),
             color="#1B2021") +
  geom_vline(aes(xintercept=0),
             color="#1B2021") +
  geom_abline(intercept=0, 
              slope=1, 
              color='grey40', 
              linetype='dashed') +
  geom_point(color="#237194") +
  stat_cor(p.accuracy=0.001, 
           r.accuracy=0.01) +
  xlab('Effect Size (Base)') +
  ylab('Effect Size (CRP)') +
  ggtitle('') +
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
      fill="white"),
    legend.position="none")

# Adjust p-values
il6_top$padj_crp <- p.adjust(il6_top$ext2_meta_p, 
                             method='fdr')

# Look at FDR 5% CpGs
paste0(
  nrow(il6_top %>% filter(padj_crp<=0.05)), 
  ' of our ', 
  nrow(il6_top), 
  ' CpGs are significant when looked up',
  ' in the CRP sensitivity analysis')

# Suppl. Tab. 4
write_csv(il6_top,
          file='../IL6_Data/Tables/ST4.csv')

##############################################################
# CRP EWAS adjusted for IL-6
# Load 1,649 CpGs previously associated with CRP
load("../IL6_Data/EWAS/CRP_CpGs.Rdata")

# Save these in our results
il6_crp <- il6_results %>% filter(cpg %in% crp_trait$cpg)
nrow(il6_crp)

# Look at effects of adjusting for IL-6
il6_crp <- il6_crp %>% 
  mutate(il6_effect=abs(crp_base_meta_beta-crp_ext2_meta_beta),
         il6_perc=abs(il6_effect*100/crp_base_meta_beta))

# Look at strongest effects
il6_crp %>% 
  arrange(desc(il6_perc)) %>% 
  select(cpg, il6_effect, il6_perc) %>% 
  head(20)

# Suppl. Fig. 1
il6_crp %>% 
  ggplot(aes(x=crp_base_meta_beta, 
             y=crp_ext2_meta_beta)) +
  geom_hline(aes(yintercept=0),
             color="#1B2021") +
  geom_vline(aes(xintercept=0),
             color="#1B2021") +
  geom_abline(intercept=0, 
              slope=1, 
              color='grey40', 
              linetype='dashed') +
  geom_point(color="#237194") +
  stat_cor(p.accuracy=0.001, 
           r.accuracy=0.01) +
  xlab('Effect Size (Base)') +
  ylab('Effect Size (IL-6)') +
  ggtitle('') +
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
      fill="white"),
    legend.position="none")

# Suppl. Tab. 5
write_csv(il6_crp,
          file='../IL6_Data/Tables/ST5.csv')

##############################################################