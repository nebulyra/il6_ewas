##############################################################
# SCRIPT 9: Sensitivity Analyses
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(ggrepel)
library(GenomicRanges)
library(ggpubr)

##############################################################
# Setup
load("../IL6_Data/Results/IL6_EWAS_Results-Full.Rdata")

# Filter out CpGs with heterogeneity of effects between cohorts
il6_results <- il6_results %>% 
  filter(base_het_i2 <= 90)

# Save top hits for sensitivity analyses
il6_top <- il6_results %>% filter(base_meta_padj_fdr <= 0.05)
dim(il6_top)

##############################################################
# Extended Cell Count Sensitivity Analysis
# Save effect size changes
il6_top <- il6_top %>% 
  mutate(cell_effect = abs(base_meta_beta - ext3_meta_beta),
         cell_perc = abs(cell_effect * 100 / base_meta_beta))

# Look at top affected CpGs
il6_top %>% 
  arrange(desc(cell_perc)) %>% 
  select(cpg, smoke_effect, smoke_perc) %>% 
  head(20)
  
# Save adjusted p-values
il6_top$padj_cell <- p.adjust(il6_top$ext3_meta_p, 
                               method='fdr')

# Look at remaining CpGs
paste0(
  nrow(il6_top %>% filter(padj_cell <= 0.05)), 
  ' of our ', 
  nrow(il6_top), 
  ' CpGs are significant when looked up in',
  ' the cell counts sensitivity analysis')

# Look at removed CpG
il6_top %>% 
  filter(padj_cell > 0.05) %>% 
  select(cpg, cell_effect, cell_perc)

##############################################################
# Smoking Sensitivity Analysis
# Save effect size changes
il6_top <- il6_top %>% 
  mutate(smoke_effect = abs(base_meta_beta - ext1_meta_beta),
         smoke_perc = abs(smoke_effect * 100 / base_meta_beta))

# Look at top affected CpGs
il6_top %>% 
  arrange(desc(smoke_perc)) %>% 
  select(cpg, smoke_effect, smoke_perc) %>% 
  head(20)

# Save adjusted p-values
il6_top$padj_smoke <- p.adjust(il6_top$ext1_meta_p, 
                              method='fdr')

# Look at remaining CpGs
paste0(
  nrow(il6_top %>% filter(padj_smoke <= 0.05)), 
  ' of our ', 
  nrow(il6_top), 
  ' CpGs are significant when looked up in',
  ' the smoking sensitivity analysis')

# Look at removed CpGs
il6_top %>% 
  filter(padj_smoke > 0.05) %>% 
  select(cpg, smoke_effect, smoke_perc)

##############################################################
# Figures
# Figure 1a: The IL-6 effect before and after adjustment for cell counts 

il6_top %>% 
  ggplot(aes(x=base_meta_beta, 
             y=ext3_meta_beta)) +
  geom_hline(aes(yintercept=0),
             color="#1B2021") +
  geom_vline(aes(xintercept=0),
             color="#1B2021") +
  geom_abline(intercept=0, 
              slope=1, 
              color='grey40', 
              linetype='dashed') +
  geom_point(color=ifelse(il6_top$padj_smoke>0.05, 
                            "#c1c0bc", "white"),
             alpha=ifelse(il6_top$padj_smoke>0.05, 1, 0)) +
  geom_point(color=ifelse(il6_top$padj_cell>0.05, 
                            "#B3AF8F", "white"),
             alpha=ifelse(il6_top$padj_cell>0.05, 1, 0)) +
  geom_point(color=ifelse(il6_top$padj_cell<=0.05 & 
                            il6_top$padj_smoke<=0.05, 
                            "#237194", "white"),
             alpha=ifelse(il6_top$padj_cell<=0.05 & 
                            il6_top$padj_smoke<=0.05, 1, 0)) +
  stat_cor(data=il6_top %>% 
                filter(padj_cell<=0.05 & padj_smoke<=0.05),
           color="#237194",
           p.accuracy=0.001, 
           r.accuracy=0.01) +
  xlab('Effect Size (Base)') +
  ylab('Effect Size (IDOLext)') +
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

# Figure 1b: The IL-6 effect before and after adjustment for smoking

il6_top %>% 
  ggplot(aes(x=base_meta_beta, 
             y=ext1_meta_beta)) +
  geom_hline(aes(yintercept=0),
             color="#1B2021") +
  geom_vline(aes(xintercept=0),
             color="#1B2021") +
  geom_abline(intercept=0, 
              slope=1, 
              color='grey40', 
              linetype='dashed') +
  geom_point(color=ifelse(il6_top$padj_cell>0.05, 
                          "#c1c0bc", "white"),
             alpha=ifelse(il6_top$padj_cell>0.05, 1, 0)) +
  geom_point(color=ifelse(il6_top$padj_smoke>0.05, 
                          "#B3AF8F", "white"),
             alpha=ifelse(il6_top$padj_smoke>0.05, 1, 0)) +
  geom_point(color=ifelse(il6_top$padj_smoke<=0.05 & 
                            il6_top$padj_cell<=0.05, 
                          "#237194", "white"),
             alpha=ifelse(il6_top$padj_smoke<=0.05 & 
                            il6_top$padj_cell<=0.05, 1, 0)) +
  stat_cor(data=il6_top %>% 
             filter(padj_smoke<=0.05 & padj_cell<=0.05),
           color="#237194",
           p.accuracy=0.001, 
           r.accuracy=0.01) +
  xlab('Effect Size (Base)') +
  ylab('Effect Size (Smoking)') +
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

##############################################################
# Save Results
# Suppl. Tab. 2
write_csv(il6_top %>% arrange(cpg_chr_hg19,
                              cpg_start_hg19),
          file='../IL6_Data/Tables/ST2.csv')

# Filter
il6_top <- il6_top %>% 
  filter(padj_smoke <= 0.05 & padj_cell <= 0.05)
dim(il6_top)

# Save Rdata
save(il6_top, file="../IL6_Data/Results/IL6_EWAS_Results-Top.Rdata")

##############################################################