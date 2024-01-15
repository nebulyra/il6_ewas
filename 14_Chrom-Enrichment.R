##############################################################
# SCRIPT 14: Chromatin State Enrichment
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(ggrepel)
library(GenomicRanges)
library(ggpubr)
library(DNAmArray)
library(MASS)

# Load data
load("../IL6_Data/Results/IL6_EWAS_Results-Full.Rdata")
load("../IL6_Data/Results/IL6_EWAS_Results-Top.Rdata")

##############################################################
# Enrichment analysis
# States
chrom_states <- c('1_TssA', '2_TssAFlnk', '3_TxFlnk',
                  '4_Tx', '5_TxWk', '6_EnhG',
                  '7_Enh', '8_ZNF/Rpts', '9_Het',
                  '10_TssBiv', '11_BivFlnk', '12_EnhBiv',
                  '13_ReprPC', '14_ReprPCWk', '15_Quies')

# Logistic regression for each state (sig ~ state)
for(i in chrom_states){
  res_road <- il6_results %>% 
    # Indicators for 5% FDR and state
    mutate(sig=ifelse(il6_results$cpg %in% il6_top$cpg, 1, 0),
           state=ifelse(E062 == i, 1, 0)) %>% 
    mutate(state=ifelse(is.na(state), 0, state))
  
  # Logistic regression
  x <- glm(state ~ sig, 
           family=binomial, data=res_road)
  
  # Save results
  out <- c(coef(summary(x))[2,],
           exp(cbind(coef(x), confint.default(x)))[2,])
  names(out) <- c('logOR', 'SE', 'z', 'p', 
                  'OR', 'low_CI', 'upp_CI')
  out <- as.data.frame(t(out))
  out$Trait = i
  out <- out %>% 
    dplyr::select(Trait, OR, logOR, low_CI, upp_CI, z, p) 
  
  if(i==chrom_states[1]){
    res <- out
  } else {
    res <- rbind(res, out)
  }
}

# Suppl. Tab. 6
write_csv(res %>% arrange(p, OR),
          file="../IL6_Data/Tables/ST6.csv")

##############################################################
# Fig. 3a
# Labels
res$Trait <- c("(15) Quies", "(14) ReprPCWk", "(13) ReprPC", 
               "(12) EnhBiv", "(11) BivFlnk",
               "(10) TssBiv", "(9) Het", "(8) ZNF/Rpts", 
               "(7) Enh", "(6) EnhG", 
               "(5) TxWk", "(4) Tx",  "(3) TxFlnk", 
               "(2) TssAFlnk", "(1) TssA")

# Initialize
res <- res %>% 
  mutate(
    loglowCI=log(low_CI),
    loguppCI=log(upp_CI),
    padj=p.adjust(p, method='fdr')
  ) %>% 
  filter(OR<200)
res$fill <- ifelse(res$padj<0.05, "Enriched", "Not Enriched")
res$invlogOR <- -res$logOR

# Fig. 3a
png("../IL6_Data/Chromatin_States/OR-Plot_Chromatin.png", 
    width=1000, height=800)
res %>% 
  arrange(logOR) %>% 
  ggplot(aes(
    y=logOR, 
    x=reorder(Trait,-invlogOR), 
    fill=fill,
    ymin=loglowCI,
    ymax=loguppCI)) +
  geom_hline(yintercept=0, 
             size=1, 
             linetype='dashed') +
  geom_errorbar(width=0.5,
                size=1,
                position=position_dodge(width=0.9)) +
  geom_point(
    size=10,
    shape=21,
    stroke=1.2,
    position=position_dodge(width=0.9)) +
  scale_fill_manual(values=c('#237194', '#ACAFAF'), 
                    name="") +
  geom_vline(xintercept=0.4, 
             size=1) +
  ylab('log(OR)') + xlab('') + ylim(c(-5,5)) +
  coord_flip() +
  theme(
    axis.text=element_text(
      size=32, 
      color='#373334'),
    axis.title=element_text(
      size=32, 
      hjust=0.5, 
      color='#373334'),
    text=element_text(size=32),
    plot.title=element_text(
      size=16, 
      hjust=0.5,
      face='bold', 
      color='#66cec8'),
    panel.background=element_rect(
      fill='white'),
    panel.grid.major=element_line(
      color='grey95'),
    panel.grid.minor=element_line(
      color='grey95'),
    plot.background=element_rect(
      fill='white'))
dev.off()

##############################################################
# Save CpGs for HOMER
homer <- il6_top %>% 
  select(cpg, cpg_chr_hg19, cpg_start_hg19, 
         cpg_end_hg19, cpg_strand)

write_tsv(homer, file='..IL6_Data/Homer/IL6_HomerInput.tsv')

##############################################################