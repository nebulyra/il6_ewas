##############################################################
# SCRIPT 12: EWAS Enrichment
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(MASS)

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
# EWAS Catalog (2,059,897 associations) - Process Data
# Load data
ewas_res <- read_tsv('../IL6_Data/EWAS/ewascatalog-results.txt.gz')
ewas_stu <- read_tsv('../IL6_Data/EWAS/ewascatalog-studies.txt.gz')

# Merge into one data frame
ewas_cat <- left_join(ewas_res, ewas_stu, by='StudyID')

# Remove those with N<500
ewas_cat <- ewas_cat %>% filter(N>=500)

# Remove unpublished results (no PMID)
ewas_cat <- ewas_cat[!is.na(as.numeric(ewas_cat$PMID)),]

# Keep only whole blood or leukocytes
table(ewas_cat$Tissue)
ewas_cat <- ewas_cat %>% 
  filter(Tissue %in% c("blood", "Blood", 
                       "CD4+ T-cells, whole blood",
                       "CD4+ T-cells, Whole blood",
                       "Cord blood, whole blood",
                       "Leukocytes", "Peripheral blood",
                       "whole blood", "Whole blood",
                       "Whole Blood",
                       "Whole blood, breast tissue",
                       "Whole blood, CD4+ T-cells, CD14+ monocytes",
                       "Whole blood, CD4+ T cells",
                       "Whole blood, CD4+ T cells, CD14+ monocytes",
                       "Whole blood, cord blood",
                       "Whole blood, heel prick blood spots"
))

# Keep only adults
table(ewas_cat$Age)
ewas_cat <- ewas_cat %>% filter(Age %in% c(
  "Adults", 
  "Adults (18-65 years)", 
  "Geriatrics"))

##############################################################
# EWAS catalog - Key
write_csv(as.data.frame(table(ewas_cat$Trait)), 
          file='../IL6_Data/EWAS/EWAS_cat_Trait.csv')

# This key was then manually created to collapse similar traits
# e.g. BMI and body mass index --> BMI
# See Suppl. Tab. 3 for more details
# Load in key
trait_key <- read_csv('../IL6_Data/EWAS/EWAS_cat_Trait.csv')

# Merge
ewas_cat <- left_join(ewas_cat, trait_key, by = 'Trait')

# Keep only CpGs in our results
ewas_cat <- ewas_cat %>% filter(CpG %in% il6_results$cpg)

ewas_full <- ewas_cat %>% 
  dplyr::select(cpg = CpG, Trait = Trait_Simple)
dim(ewas_full)
length(unique(ewas_full$cpg))

##############################################################
# EWAS Atlas (617,084 associations) - Process Data
# Load data
ewas_atlas <- read_tsv('../IL6_Data/EWAS/EWAS_Atlas_associations.tsv')
ewas_stu <- read_tsv('../IL6_Data/EWAS/EWAS_Atlas_cohorts.tsv')

# Subset relevant variables
ewas_stu <- ewas_stu %>% 
  dplyr::select(study_ID, tissue, mean_age, sample_size) 

# Remove duplicates
ewas_stu <- ewas_stu[!duplicated(ewas_stu$study_ID),]

# Merge
ewas_atlas <- left_join(ewas_atlas, ewas_stu, by='study_ID')

# Remove if N<500
ewas_atlas <- ewas_atlas %>% 
  filter(sample_size >= 500)

# Remove unpublished (no PMID) or present in EWAS catalog
ewas_atlas <- ewas_atlas %>% 
  filter(!is.na(as.numeric(PMID))) %>% 
  filter(!PMID %in% ewas_cat$PMID)

# Keep only whole blood or leukocytes
table(ewas_atlas$tissue)
ewas_atlas <- ewas_atlas %>% 
  filter(tissue %in% c("blood", "blood spot",
                       "bloodspot", "buffy coat",
                       "bufy coat", "leukocyte",
                       "leukocytes", "peripheral blood",
                       "white blood cells", "whole blood"
))

# Keep only adults
ewas_atlas <- ewas_atlas %>% filter(mean_age >= 18)

##############################################################
# EWAS atlas - Key
write_csv(as.data.frame(table(ewas_atlas$trait)), 
          file='../IL6_Data/EWAS/EWAS_atlas_Trait.csv')

# This key was then manually created to collapse similar traits
# e.g. BMI and body mass index --> BMI
# See Suppl. Tab. 3 for more details
# Load in key
trait_key <- read_csv('../IL6_Data/EWAS/EWAS_atlas_Trait.csv')
ewas_atlas <- left_join(ewas_atlas, trait_key, by = 'trait')

# Keep only our CpGs
ewas_atlas <- ewas_atlas %>% 
  filter(probe_ID %in% il6_results$cpg)
dim(ewas_atlas)
length(unique(ewas_atlas$probe_ID))

# Merge
ewas_full <- rbind(ewas_full, ewas_atlas %>% 
                   dplyr::select(cpg=probe_ID, Trait=Trait_Simple))
dim(ewas_full)
length(unique(ewas_full$cpg))

##############################################################
# CRP EWAS
# Load data
crp_ewas <- read_csv('../IL6_Data/EWAS/CRP-EWAS_CpG-list.csv')

# Save only our CpGs
crp_ewas <- crp_ewas %>% 
  filter(ID %in% il6_results$cpg)

# Add trait variable
crp_ewas <- data.frame(cpg = crp_ewas$ID, Trait="CRP")

# Merge
ewas_full <- rbind(ewas_full, crp_ewas)
dim(ewas_full)
length(unique(ewas_full$cpg))

# Save CRP CpGs
crp_trait <- ewas_full %>% 
  filter(Trait == "CRP") %>% 
  dplyr::select(cpg)
nrow(crp_trait)

save(crp_trait, file='../IL6_Data/EWAS/CRP_CpGs.Rdata')

##############################################################
# Enrichment analysis
# Logistic regression of sig ~ trait if >100 CpGs found
trait_list <- unique(ewas_full$Trait)

for(i in trait_list){
  cpg_trait <- (ewas_full %>% filter(Trait == i))$cpg
  n <- length(cpg_trait)
  
  if(n>=100){
    df <- il6_results %>% 
      mutate(
        trait = ifelse(cpg %in% cpg_trait, 1, 0),
        sig = ifelse(cpg %in% il6_top$cpg, 1, 0)
      )
    nFound = sum(df$trait)
    nSig = sum((df %>% filter(sig==1))$trait)
    
    trait_df <- df %>% filter(trait==1)
    save(trait_df, 
         file=paste0('../IL6_Data/EWAS/Trait/', i, '.Rdata'))
    
    x <- glm(trait~sig, family=binomial, data=df)
    out <- c(coef(summary(x))[2,],
             exp(cbind(coef(x),
                       confint.default(x)))[2,])
    names(out) <- c('logOR', 'SE', 'z', 'p', 
                    'OR', 'low_CI', 'upp_CI')
    out <- as.data.frame(t(out))
    out$Trait = i
    out$nCpG = n
    out$nFound = nFound
    out$nSig = nSig
    out <- out %>% 
      dplyr::select(Trait, OR, logOR, low_CI, upp_CI, 
                    z, p, nCpG, nFound, nSig) 
    out$padj <- p.adjust(out$p, method='fdr')
    if(i==trait_list[1]){
      res <- out
    } else {
      res <- rbind(out, res)
    }
  }    
}

# Suppl. Tab. 3
write_csv(res %>% 
          arrange(p, OR),
          file="../IL6_Data/Tables/ST3.csv")

##############################################################
# Figure 2: Forest plot
# Initialize
res <- unique(res)
res <- res %>% 
  dplyr::select(Trait, OR, logOR,
                low_CI, upp_CI, p, 
                nCpG, nFound, nSig) %>% 
  mutate(
    loglowCI=log(low_CI),
    loguppCI=log(upp_CI),
    padj=p.adjust(p, method='fdr'))

# Look at top traits
top_trait <- (res %>% arrange(padj) 
              %>% head(20))$Trait
top_trait
res$label <- paste0(res$Trait, 
                    ' (', res$nSig, ', ', 
                    round(res$nSig/res$nFound*100,0), '%)')
res$invlogOR <- -res$logOR

# Plot
png("../IL6_Data/EWAS/OR-Plot-Enrichment.png", 
    width=1000, height=800)
res %>% 
  filter(Trait %in% top_trait) %>% 
  arrange(OR) %>% 
  ggplot(aes(
    y=logOR, 
    x=reorder(label,-invlogOR), 
    ymin=0,
    ymax=loguppCI+0.5)) +
  geom_errorbar(width=0.5,
                size=1,
                position=position_dodge(width=0.9)) +
  geom_point(
    size=9,
    shape=21,
    fill='#237194',
    stroke=1.2,
    position=position_dodge(width=0.9)) +
  ylab('log(OR)') + xlab('') +
  coord_flip() +
  theme(
    axis.text=element_text(
      size=32,
      color='#373334'),
    axis.title=element_text(
      size=32, 
      hjust=0.5, 
      color='#373334'),
    text=element_text(
      size=32),
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
