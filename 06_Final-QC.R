##############################################################
# SCRIPT 6: Final QC of EWAS summary data
##############################################################
# Setup
rm(list=ls())
library(ggpubr)
library(DNAmArray)
library(lattice)
library(ggrepel)
library(tidyverse)

# Load zhou probes
anno <- read_tsv("../Shared_Data/Manifests/HM450.hg19.manifest.tsv.gz")
anno <- anno %>% dplyr::select(
    cpg = probeID, chr = CpG_chrm, start = CpG_beg, MASK_general) %>% 
  mutate(chr = as.numeric(sub("chr", "", chr)))
maskProbes <- anno[anno$MASK_general == TRUE,]$cpg

# Save autosomal probes
autoProbes <- anno[anno$chr %in% paste0("chr", 1:22),]

# Load probes in ENCODE blacklist regions
load("../Shared_Data/ENCODE_Blacklist/ENCODEBlacklist_CpGomit-450K.RData")

##############################################################
# IL6 EWAS - Base Model
# Load data
load('../IL6_Data/Processing/ALL-3-BACON_ln-il6_base.Rdata')

# Save Results
# Function
save_results <- function(cohort){
  assign(paste0("df_", cohort),
         bacon_df %>% dplyr::select(cpg,
             beta = paste0("beta_", cohort),
             se = paste0("SE_", cohort),
             p = paste0("p_", cohort),
             N = paste0("N_", cohort)) %>% 
           mutate(beta = as.numeric(beta),
             se = as.numeric(se), 
             p = as.numeric(p),
             N = as.numeric(N)) %>% 
           filter(N>50), envir=.GlobalEnv)
}

# Save for each cohort
save_results("lls")
write_tsv(df_lls, file="../IL6_Data/Processing/Meta_Input/OUT_il6_base-LLS.tsv")
save_results("kora")
write_tsv(df_kora, file="../IL6_Data/Processing/Meta_Input/OUT_il6_base-KORA.tsv")
save_results("ntr")
write_tsv(df_ntr, file="../IL6_Data/Processing/Meta_Input/OUT_il6_base-NTR.tsv")

##############################################################
# Quality Control Plots
# QQ plots 
# Function
draw_qq <- function(cohort){
  n <- nrow(get(paste0("df_",cohort)))
  exp_p <- -log10((rank(get(paste0("df_", cohort))$p,
                        ties.method="first") - 0.5) / n)
  obs_p <- -log10(get(paste0("df_", cohort))$p)
  df <- data.frame(e=exp_p, o=obs_p)
  
  assign(paste0("plot_", cohort), df %>% 
            ggplot() +
            geom_abline(aes(intercept=0, slope=1), color="#D62839") +
            geom_point(aes(x=e, y=o), alpha=0.5, size=1, color="#1B2021") +
            ylim(0,NA) + xlim(0,NA) + ggtitle(toupper(cohort)) +
            xlab("Expected") + ylab("Observed") +
  theme(axis.text = element_text(size=9, color="#1B2021"),
    axis.title = element_text(size=11,hjust=0.5, color="#1B2021"),
    plot.title = element_text(size=16, hjust=0.5,face="bold", color = "#548687"),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(color="#1B2021",fill=NA),
    panel.grid.major = element_line(color="grey95"),
    panel.grid.minor = element_line(color="grey95"),
    plot.background = element_rect(fill="white")),
  envir = .GlobalEnv)
}

# Draw for each cohort
draw_qq("lls")
draw_qq("kora")
draw_qq("ntr")

ggarrange(plot_lls, plot_kora, plot_ntr, plot_lls,
          ncol = 2, nrow = 2)

# Volcano Plots
# Scale
min <- as.numeric(-max(c(abs(bacon_df$beta_lls), 
                         abs(bacon_df$beta_kora),
                         abs(bacon_df$beta_ntr)),na.rm=TRUE) - 0.005)

max <- as.numeric(max(c(abs(bacon_df$beta_lls), 
                        abs(bacon_df$beta_kora),
                        abs(bacon_df$beta_ntr)),na.rm=TRUE) + 0.005)

p_max <- as.numeric(-log10(min(c(bacon_df$p_lls,
                                 bacon_df$p_kora,
                                 bacon_df$p_ntr),na.rm=TRUE))) + 2

# Function
draw_volcano <- function(cohort){
  sig_df <- get(paste0("df_", cohort)) %>% 
    mutate(padj = p.adjust(p, method="fdr")) %>% 
    filter(padj <= 0.05)
  
  if(nrow(sig_df) >= 1){
    sig_limit <- max(sig_df$p)
  } else {
    sig_limit <- 10E-07
  }
  
  assign(paste0("plot_", cohort), 
         get(paste0('df_', cohort)) %>% 
           ggplot(aes(x = beta, y = -log10(p))) +
           geom_hline(yintercept = -log10(sig_limit), linetype = "dashed",
             size = 0.5, color = "#ACAFAF") +
           geom_vline(xintercept = -0.05, linetype = "dashed",
             size = 0.5, color = "#ACAFAF") +
           geom_vline(xintercept = 0.05, linetype = "dashed",
                      size = 0.5, color = "#ACAFAF") +
           geom_point(color = ifelse(get(paste0('df_', cohort))$p <= sig_limit, 
                          "#D62839", "#ACAFAF"), size = 0.8) +
           xlim(min, max) + ylim(0, p_max) + ggtitle(toupper(cohort)) +
           ylab(bquote(-log[10]~"p")) + xlab("Effect Size") +
  theme(axis.text = element_text(size=9, color="#1B2021"),
    axis.title = element_text(size=11, hjust=0.5, color="#1B2021"),
    plot.title = element_text(size=16, hjust=0.5,face="bold", color = "#548687"),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(color="#1B2021",fill=NA),
    panel.grid.major = element_line(color="grey95"),
    panel.grid.minor = element_line(color="grey95"),
    plot.background = element_rect(fill="white")),
  envir = .GlobalEnv)
}

# Draw for each cohort
draw_volcano("lls")
draw_volcano("kora")
draw_volcano("ntr")

ggarrange(plot_lls, plot_kora, plot_ntr, plot_lls,
  ncol = 2, nrow = 2)

# Manhattan plots
draw_man <- function(cohort){
  sig_df <- get(paste0("df_", cohort)) %>% 
    mutate(padj = p.adjust(p, method="fdr")) %>% 
    filter(padj <= 0.05)
  
  if(nrow(sig_df) >= 1){ sig_limit <- max(sig_df$p)
  } else { sig_limit <- 10E-07 }

  limma_man <- get(paste0("df_", cohort)) %>% 
    arrange(chr, start)
  
  limma_man$start_cum <- NA
  s <- 0
  nbp <- c()
  
  for(i in unique(limma_man$chr)){
    nbp[i] <- max(limma_man[limma_man$chr == i,]$start)
    limma_man[limma_man$chr == i, "startcum"] <- limma_man[limma_man$chr == i, 
                                                           "start"] + s
    s <- s + nbp[i]
  }
  
  axis.set <- limma_man %>% group_by(chr) %>% 
    summarize(center = (max(startcum) + min(startcum)) / 2)
  
  assign(paste0('plot_', cohort),
         ggplot(limma_man, aes(x=startcum, y=-log10(p),
             color=as.factor(chr))) +
         geom_point(size=0.8) + 
         geom_hline(yintercept = -log10(sig_limit),
             color = "#D62839", size = 0.5, linetype = "dashed") +
         scale_x_continuous(label = axis.set$chr, breaks = axis.set$center) +
         scale_y_continuous(expand = c(0,0), limits = c(0, p_max)) +
         scale_size_continuous(range = c(0.5, 3)) +
         scale_color_manual(values = rep(c("#ACAFAF", "#548687"), 11)) +
         xlab("") + ylab(bquote(-log[10]~"p")) + ggtitle(toupper(cohort)) +
         theme(legend.position = "none",
              panel.background = element_rect(fill="white"),
              panel.border = element_rect(color="#1B2021", fill=NA),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.text.y = element_text(size = 9, color = '#1B2021'),
              axis.text.x = element_text(angle = 90, size = 4, 
                           vjust = 0.5, color = '#1B2021'),
              axis.title = element_text(size=11, color = '#1B2021'),
              plot.title = element_text(size = 12, hjust=0.5,
                          color = '#548687', face = 'bold')),
              envir=.GlobalEnv)
}

# Draw for each cohort
draw_man("lls")
draw_man("kora")
draw_man("ntr")

ggarrange(plot_lls, plot_kora, plot_ntr, plot_lls,
  ncol = 2, nrow = 2)

# Boxplots
# Effect sizes
bacon_df <- bacon_df %>% 
  mutate(
    beta_LLS = as.numeric(beta_lls),
    beta_KORA = as.numeric(beta_kora),
    beta_NTR = as.numeric(beta_ntr)
  )

beta_df <- pivot_longer(
  data = bacon_df,
  cols = c("beta_LLS", "beta_KORA", "beta_NTR"),
  names_to = "cohort",
  names_prefix = "beta_",
  values_to = "beta",
  values_drop_na = TRUE
)

beta_df %>% 
  ggplot(aes(x = cohort, y = beta, fill = cohort)) +
  scale_fill_manual(values = c("#e2455f", "#eaa87d", "#66cec8")) +
  geom_boxplot(alpha = 0.8) +
  ggtitle("Effect size distributions by cohort") +
  ylab("\u03b2") + xlab("Cohort") +
  theme( 
    panel.background = element_rect(
      fill="white"),
    panel.border = element_rect(
      color="#1B2021", fill=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "#1B2021"),
    axis.text.x = element_text(size = 12, color = "#1B2021"),
    axis.title = element_text(size=11, color = "#1B2021"),
    plot.title = element_text(size = 16, hjust=0.5,
                              color = "#548687", face = "bold"))

# SEs
bacon_df <- bacon_df %>% 
  mutate(
    SE_LLS = as.numeric(SE_lls),
    SE_KORA = as.numeric(SE_kora),
    SE_NTR = as.numeric(SE_ntr)
  )

SE_df <- pivot_longer(
  data = bacon_df,
  cols = c("SE_LLS", "SE_KORA", "SE_NTR"),
  names_to = "cohort",
  names_prefix = "SE_",
  values_to = "SE",
  values_drop_na = TRUE
)

SE_df %>% 
  ggplot(aes(x = cohort, y = SE, fill = cohort)) +
  scale_fill_manual(values = c("#e2455f", "#eaa87d", "#66cec8")) +
  geom_boxplot(alpha = 0.6) +
  ggtitle("SE distributions by cohort") +
  ylab("SE") + xlab("Cohort") +
  theme( 
    panel.background = element_rect(
      fill="white"),
    panel.border = element_rect(
      color="#1B2021", fill=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "#1B2021"),
    axis.text.x = element_text(size = 12, color = "#1B2021"),
    axis.title = element_text(size=11, color = "#1B2021"),
    plot.title = element_text(size = 16, hjust=0.5,
                              color = "#548687", face = "bold"))

##############################################################
# IL6 EWAS - Smoking Sensitivity Analysis
# Load data
load("../IL6_Data/Processing/ALL-3-BACON_ln-il6_ext1.Rdata")

# Save Results
save_results("lls")
write_tsv(df_lls, file="../IL6_Data/Processing/Meta_Input/OUT_il6_ext1-LLS.tsv")
save_results("kora")
write_tsv(df_kora, file="../IL6_Data/Processing/Meta_Input/OUT_il6_ext1-KORA.tsv")
save_results("ntr")
write_tsv(df_ntr, file="../IL6_Data/Processing/Meta_Input/OUT_il6_ext1-NTR.tsv")

##############################################################
# Quality Control Plots
# QQ plots 
draw_qq("lls")
draw_qq("kora")
draw_qq("ntr")

ggarrange(plot_lls, plot_kora, plot_ntr, plot_lls,
          ncol = 2, nrow = 2)

# Volcano Plots
# Scales
min <- as.numeric(-max(c(abs(bacon_df$beta_lls), 
                         abs(bacon_df$beta_kora),
                         abs(bacon_df$beta_ntr)),na.rm=TRUE) - 0.005)

max <- as.numeric(max(c(abs(bacon_df$beta_lls), 
                        abs(bacon_df$beta_kora),
                        abs(bacon_df$beta_ntr)),na.rm=TRUE) + 0.005)

p_max <- as.numeric(-log10(min(c(bacon_df$p_lls,
                                 bacon_df$p_kora,
                                 bacon_df$p_ntr),na.rm=TRUE))) + 2

# Draw for each cohort
draw_volcano("lls")
draw_volcano("kora")
draw_volcano("ntr")

ggarrange(plot_lls, plot_kora, plot_ntr, plot_lls,
  ncol = 2, nrow = 2)

# Manhattan plots
draw_man("lls")
draw_man("kora")
draw_man("ntr")

ggarrange(plot_lls, plot_kora, plot_ntr, plot_lls,
  ncol = 2, nrow = 2)

# Boxplots
# Effect sizes
bacon_df <- bacon_df %>% 
  mutate(
    beta_LLS = as.numeric(beta_lls),
    beta_KORA = as.numeric(beta_kora),
    beta_NTR = as.numeric(beta_ntr)
  )

beta_df <- pivot_longer(
  data = bacon_df,
  cols = c("beta_LLS", "beta_KORA", "beta_NTR"),
  names_to = "cohort",
  names_prefix = "beta_",
  values_to = "beta",
  values_drop_na = TRUE
)

beta_df %>% 
  ggplot(aes(x = cohort, y = beta, fill = cohort)) +
  scale_fill_manual(values = c("#e2455f", "#eaa87d","#66cec8")) +
  geom_boxplot(alpha = 0.8) +
  ggtitle("Effect size distributions by cohort") +
  ylab("\u03b2") + xlab("Cohort") +
  theme( 
    panel.background = element_rect(
      fill="white"),
    panel.border = element_rect(
      color="#1B2021", fill=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "#1B2021"),
    axis.text.x = element_text(size = 12, color = "#1B2021"),
    axis.title = element_text(size=11, color = "#1B2021"),
    plot.title = element_text(size = 16, hjust=0.5,
                              color = "#548687", face = "bold"))

# SEs
bacon_df <- bacon_df %>% 
  mutate(
    SE_LLS = as.numeric(SE_lls),
    SE_KORA = as.numeric(SE_kora),
    SE_NTR = as.numeric(SE_ntr)
  )

SE_df <- pivot_longer(
  data = bacon_df,
  cols = c("SE_LLS", "SE_KORA", "SE_NTR"),
  names_to = "cohort",
  names_prefix = "SE_",
  values_to = "SE",
  values_drop_na = TRUE
)

SE_df %>% 
  ggplot(aes(x = cohort, y = SE, fill = cohort)) +
  scale_fill_manual(values = c("#e2455f", "#eaa87d", "#66cec8")) +
  geom_boxplot(alpha = 0.6) +
  ggtitle("SE distributions by cohort") +
  ylab("SE") + xlab("Cohort") +
  theme( 
    panel.background = element_rect(
      fill="white"),
    panel.border = element_rect(
      color="#1B2021", fill=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "#1B2021"),
    axis.text.x = element_text(size = 12, color = "#1B2021"),
    axis.title = element_text(size=11, color = "#1B2021"),
    plot.title = element_text(size = 16, hjust=0.5,
                              color = "#548687", face = "bold"))

##############################################################
# IL6 EWAS - CRP Sensitivity Analysis
# Load data
load("../IL6_Data/Processing/ALL-3-BACON_ln-il6_ext2.Rdata")

# Save Results
save_results("lls")
write_tsv(df_lls, file="../IL6_Data/Processing/Meta_Input/OUT_il6_ext2-LLS.tsv")
save_results("kora")
write_tsv(df_kora, file="../IL6_Data/Processing/Meta_Input/OUT_il6_ext2-KORA.tsv")
save_results("ntr")
write_tsv(df_ntr, file="../IL6_Data/Processing/Meta_Input/OUT_il6_ext2-NTR.tsv")

##############################################################
# Quality Control Plots
# QQ plots 
draw_qq("lls")
draw_qq("kora")
draw_qq("ntr")

ggarrange(plot_lls, plot_kora, plot_ntr, plot_lls,
          ncol = 2, nrow = 2)

# Volcano Plots
# Scales
min <- as.numeric(-max(c(abs(bacon_df$beta_lls), 
                         abs(bacon_df$beta_kora),
                         abs(bacon_df$beta_ntr)),na.rm=TRUE) - 0.005)

max <- as.numeric(max(c(abs(bacon_df$beta_lls), 
                        abs(bacon_df$beta_kora),
                        abs(bacon_df$beta_ntr)),na.rm=TRUE) + 0.005)

p_max <- as.numeric(-log10(min(c(bacon_df$p_lls,
                                 bacon_df$p_kora,
                                 bacon_df$p_ntr),na.rm=TRUE))) + 2

# Draw for each cohort
draw_volcano("lls")
draw_volcano("kora")
draw_volcano("ntr")

ggarrange(plot_lls, plot_kora, plot_ntr, plot_lls,
  ncol = 2, nrow = 2)

# Manhattan plots
draw_man("lls")
draw_man("kora")
draw_man("ntr")

ggarrange(plot_lls, plot_kora, plot_ntr, plot_lls,
  ncol = 2, nrow = 2)

# Boxplots
# Effect sizes
bacon_df <- bacon_df %>% 
  mutate(
    beta_LLS = as.numeric(beta_lls),
    beta_KORA = as.numeric(beta_kora),
    beta_NTR = as.numeric(beta_ntr)
  )

beta_df <- pivot_longer(
  data = bacon_df,
  cols = c("beta_LLS", "beta_KORA", "beta_NTR"),
  names_to = "cohort",
  names_prefix = "beta_",
  values_to = "beta",
  values_drop_na = TRUE
)

beta_df %>% 
  ggplot(aes(x = cohort, y = beta, fill = cohort)) +
  scale_fill_manual(values = c("#e2455f", "#eaa87d", "#66cec8")) +
  geom_boxplot(alpha = 0.8) +
  ggtitle("Effect size distributions by cohort") +
  ylab("\u03b2") + xlab("Cohort") +
  theme( 
    panel.background = element_rect(
      fill="white"),
    panel.border = element_rect(
      color="#1B2021", fill=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "#1B2021"),
    axis.text.x = element_text(size = 12, color = "#1B2021"),
    axis.title = element_text(size=11, color = "#1B2021"),
    plot.title = element_text(size = 16, hjust=0.5,
                              color = "#548687", face = "bold"))

# SEs
bacon_df <- bacon_df %>% 
  mutate(
    SE_LLS = as.numeric(SE_lls),
    SE_KORA = as.numeric(SE_kora),
    SE_NTR = as.numeric(SE_ntr)
  )

SE_df <- pivot_longer(
  data = bacon_df,
  cols = c("SE_LLS", "SE_KORA", "SE_NTR"),
  names_to = "cohort",
  names_prefix = "SE_",
  values_to = "SE",
  values_drop_na = TRUE
)

SE_df %>% 
  ggplot(aes(x = cohort, y = SE, fill = cohort)) +
  scale_fill_manual(values = c("#e2455f", "#eaa87d", "#66cec8")) +
  geom_boxplot(alpha = 0.6) +
  ggtitle("SE distributions by cohort") +
  ylab("SE") + xlab("Cohort") +
  theme( 
    panel.background = element_rect(
      fill="white"),
    panel.border = element_rect(
      color="#1B2021", fill=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "#1B2021"),
    axis.text.x = element_text(size = 12, color = "#1B2021"),
    axis.title = element_text(size=11, color = "#1B2021"),
    plot.title = element_text(size = 16, hjust=0.5,
                              color = "#548687", face = "bold"))

##############################################################
# IL6 EWAS - IDOL Extended Sensitivity Analysis
# Load data
load("../IL6_Data/Processing/ALL-3-BACON_ln-il6_ext3.Rdata")

# Save Results
save_results("lls")
write_tsv(df_lls, file="../IL6_Data/Processing/Meta_Input/OUT_il6_ext3-LLS.tsv")
save_results("kora")
write_tsv(df_kora, file="../IL6_Data/Processing/Meta_Input/OUT_il6_ext3-KORA.tsv")
save_results("ntr")
write_tsv(df_lls, file="../IL6_Data/Processing/Meta_Input/OUT_il6_ext3-NTR.tsv")

##############################################################
# Quality Control Plots
# QQ plots 
draw_qq("lls")
draw_qq("kora")
draw_qq("ntr")

ggarrange(plot_lls, plot_ntr, plot_kora, plot_lls,
          ncol = 2, nrow = 2)

# Volcano Plots
# Make scale
min <- as.numeric(-max(c(abs(bacon_df$beta_lls), 
                         abs(bacon_df$beta_kora),
                         abs(bacon_df$beta_ntr)),na.rm=TRUE) - 0.005)

max <- as.numeric(max(c(abs(bacon_df$beta_lls), 
                        abs(bacon_df$beta_kora),
                        abs(bacon_df$beta_ntr)),na.rm=TRUE) + 0.005)

p_max <- as.numeric(-log10(min(c(bacon_df$p_lls,
                                 bacon_df$p_kora,
                                 bacon_df$p_ntr),na.rm=TRUE))) + 2

# Draw for each cohort
draw_volcano("lls")
draw_volcano("kora")
draw_volcano("ntr")

ggarrange(plot_lls, plot_ntr, plot_kora, plot_lls,
  ncol = 2, nrow = 2)

# Manhattan plots
draw_man("lls")
draw_man("kora")
draw_man("ntr")

ggarrange(plot_lls, plot_ntr, plot_kora, plot_lls,
  ncol = 2, nrow = 2)

# Boxplots
# Effect sizes  
bacon_df <- bacon_df %>% 
  mutate(
    beta_LLS = as.numeric(beta_lls),
    beta_KORA = as.numeric(beta_kora),
    beta_NTR = as.numeric(beta_ntr)
  )

beta_df <- pivot_longer(
  data = bacon_df,
  cols = c("beta_LLS", "beta_KORA", "beta_NTR"),
  names_to = "cohort",
  names_prefix = "beta_",
  values_to = "beta",
  values_drop_na = TRUE
)

beta_df %>% 
  ggplot(aes(x = cohort, y = beta, fill = cohort)) +
  scale_fill_manual(values = c("#e2455f", "#eaa87d", "#66cec8")) +
  geom_boxplot(alpha = 0.8) +
  ggtitle("Effect size distributions by cohort") +
  ylab("\u03b2") + xlab("Cohort") +
  theme( 
    panel.background = element_rect(
      fill="white"),
    panel.border = element_rect(
      color="#1B2021", fill=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "#1B2021"),
    axis.text.x = element_text(size = 12, color = "#1B2021"),
    axis.title = element_text(size=11, color = "#1B2021"),
    plot.title = element_text(size = 16, hjust=0.5,
                              color = "#548687", face = "bold"))

# SEs
bacon_df <- bacon_df %>% 
  mutate(
    SE_LLS = as.numeric(SE_lls),
    SE_KORA = as.numeric(SE_kora),
    SE_NTR = as.numeric(SE_ntr)
  )

SE_df <- pivot_longer(
  data = bacon_df,
  cols = c("SE_LLS", "SE_KORA", "SE_NTR"),
  names_to = "cohort",
  names_prefix = "SE_",
  values_to = "SE",
  values_drop_na = TRUE
)

SE_df %>% 
  ggplot(aes(x = cohort, y = SE, fill = cohort)) +
  scale_fill_manual(values = c("#e2455f", "#eaa87d", "#66cec8")) +
  geom_boxplot(alpha = 0.6) +
  ggtitle("SE distributions by cohort") +
  ylab("SE") + xlab("Cohort") +
  theme( 
    panel.background = element_rect(
      fill="white"),
    panel.border = element_rect(
      color="#1B2021", fill=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "#1B2021"),
    axis.text.x = element_text(size = 12, color = "#1B2021"),
    axis.title = element_text(size=11, color = "#1B2021"),
    plot.title = element_text(size = 16, hjust=0.5,
                              color = "#548687", face = "bold"))

##############################################################
# CRP EWAS - Base Model
# Load data
load("../IL6_Data/Processing/ALL-3-BACON_ln-crp_base.Rdata")

# Save Results
save_results("lls")
write_tsv(df_lls, file="../IL6_Data/Processing/Meta_Input/OUT_crp_base-LLS.tsv")
save_results("kora")
write_tsv(df_kora, file="../IL6_Data/Processing/Meta_Input/OUT_crp_base-KORA.tsv")
save_results("ntr")
write_tsv(df_ntr, file="../IL6_Data/Processing/Meta_Input/OUT_crp_base-NTR.tsv")

##############################################################
# Quality Control Plots
# QQ plots 
draw_qq("lls")
draw_qq("kora")
draw_qq("ntr")

ggarrange(plot_lls, plot_ntr, plot_kora, plot_lls,
          ncol = 2, nrow = 2)

# Volcano Plots
# Scale
min <- as.numeric(-max(c(abs(bacon_df$beta_lls), 
                         abs(bacon_df$beta_kora),
                         abs(bacon_df$beta_ntr)),na.rm=TRUE) - 0.005)

max <- as.numeric(max(c(abs(bacon_df$beta_lls), 
                        abs(bacon_df$beta_kora),
                        abs(bacon_df$beta_ntr)),na.rm=TRUE) + 0.005)

p_max <- as.numeric(-log10(min(c(bacon_df$p_lls,
                                 bacon_df$p_kora,
                                 bacon_df$p_ntr),na.rm=TRUE))) + 2

# Draw for each cohort
draw_volcano("lls")
draw_volcano("kora")
draw_volcano("ntr")

ggarrange(plot_lls, plot_ntr, plot_kora, plot_lls,
  ncol = 2, nrow = 2)

# Manhattan plots
draw_man("lls")
draw_man("kora")
draw_man("ntr")

ggarrange(plot_lls, plot_ntr, plot_kora, plot_lls, 
  ncol = 2, nrow = 2)

# Boxplots
# Effect sizes
bacon_df <- bacon_df %>% 
  mutate(
    beta_LLS = as.numeric(beta_lls),
    beta_KORA = as.numeric(beta_kora),
    beta_NTR = as.numeric(beta_ntr)
  )

beta_df <- pivot_longer(
  data = bacon_df,
  cols = c("beta_LLS", "beta_KORA", "beta_NTR"),
  names_to = "cohort",
  names_prefix = "beta_",
  values_to = "beta",
  values_drop_na = TRUE
)

beta_df %>% 
  ggplot(aes(x = cohort, y = beta, fill = cohort)) +
  scale_fill_manual(values = c("#e2455f", "#eaa87d", "#66cec8")) +
  geom_boxplot(alpha = 0.8) +
  ggtitle("Effect size distributions by cohort") +
  ylab("\u03b2") + xlab("Cohort") +
  theme( 
    panel.background = element_rect(
      fill="white"),
    panel.border = element_rect(
      color="#1B2021", fill=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "#1B2021"),
    axis.text.x = element_text(size = 12, color = "#1B2021"),
    axis.title = element_text(size=11, color = "#1B2021"),
    plot.title = element_text(size = 16, hjust=0.5,
                              color = "#548687", face = "bold"))

# SEs
bacon_df <- bacon_df %>% 
  mutate(
    SE_LLS = as.numeric(SE_lls),
    SE_KORA = as.numeric(SE_kora),
    SE_NTR = as.numeric(SE_ntr)
  )

SE_df <- pivot_longer(
  data = bacon_df,
  cols = c("SE_LLS", "SE_KORA", "SE_NTR"),
  names_to = "cohort",
  names_prefix = "SE_",
  values_to = "SE",
  values_drop_na = TRUE
)

SE_df %>% 
  ggplot(aes(x = cohort, y = SE, fill = cohort)) +
  scale_fill_manual(values = c("#e2455f", "#eaa87d", "#66cec8")) +
  geom_boxplot(alpha = 0.6) +
  ggtitle("SE distributions by cohort") +
  ylab("SE") + xlab("Cohort") +
  theme( 
    panel.background = element_rect(
      fill="white"),
    panel.border = element_rect(
      color="#1B2021", fill=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "#1B2021"),
    axis.text.x = element_text(size = 12, color = "#1B2021"),
    axis.title = element_text(size=11, color = "#1B2021"),
    plot.title = element_text(size = 16, hjust=0.5,
                              color = "#548687", face = "bold"))

##############################################################
# CRP EWAS - IL6 Sensitivity Analysis
# Load data
load("../IL6_Data/Processing/ALL-3-BACON_ln-crp_ext2.Rdata")

# Save Results
save_results("lls")
write_tsv(df_lls, file="../IL6_Data/Processing/Meta_Input/OUT_crp_ext2-LLS.tsv")
save_results("kora")
write_tsv(df_kora, file="../IL6_Data/Processing/Meta_Input/OUT_crp_ext2-KORA.tsv")
save_results("ntr")
write_tsv(df_ntr, file="../IL6_Data/Processing/Meta_Input/OUT_crp_ext2-NTR.tsv")

##############################################################
# Quality Control Plots
# QQ plots 
draw_qq("lls")
draw_qq("kora")
draw_qq("ntr")

ggarrange(plot_lls, plot_ntr, plot_kora, plot_lls,
          ncol = 2, nrow = 2)

# Volcano Plots
# Scales
min <- as.numeric(-max(c(
  abs(bacon_df$beta_lls), 
  abs(bacon_df$beta_kora),
  abs(bacon_df$beta_ntr)),
  na.rm=TRUE) - 0.005)

max <- as.numeric(max(c(
  abs(bacon_df$beta_lls), 
  abs(bacon_df$beta_kora),
  abs(bacon_df$beta_ntr)),
  na.rm=TRUE) + 0.005)

p_max <- as.numeric(-log10(min(c(
  bacon_df$p_lls,
  bacon_df$p_kora,
  bacon_df$p_ntr),na.rm=TRUE))) + 2

# Draw for each cohort
draw_volcano("lls")
draw_volcano("kora")
draw_volcano("ntr")

ggarrange(plot_lls, plot_ntr, plot_kora, plot_lls,  
  ncol = 2, nrow = 2)

# Manhattan plots
draw_man("lls")
draw_man("kora")
draw_man("ntr")


ggarrange(plot_lls, plot_ntr, plot_kora, plot_lls,  
  ncol = 2, nrow = 2)

# Boxplots
# Effect sizes
bacon_df <- bacon_df %>% 
  mutate(
    beta_LLS = as.numeric(beta_lls),
    beta_KORA = as.numeric(beta_kora),
    beta_NTR = as.numeric(beta_ntr)
  )

beta_df <- pivot_longer(
  data = bacon_df,
  cols = c("beta_LLS", "beta_KORA", "beta_NTR"),
  names_to = "cohort",
  names_prefix = "beta_",
  values_to = "beta",
  values_drop_na = TRUE
)

beta_df %>% 
  ggplot(aes(x = cohort, y = beta, fill = cohort)) +
  scale_fill_manual(values = c("#e2455f", "orange", "#66cec8")) +
  geom_boxplot(alpha = 0.8) +
  ggtitle("Effect size distributions by cohort") +
  ylab("\u03b2") + xlab("Cohort") +
  theme( 
    panel.background = element_rect(
      fill="white"),
    panel.border = element_rect(
      color="#1B2021", fill=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "#1B2021"),
    axis.text.x = element_text(size = 12, color = "#1B2021"),
    axis.title = element_text(size=11, color = "#1B2021"),
    plot.title = element_text(size = 16, hjust=0.5,
                              color = "#548687", face = "bold"))

# SEs
bacon_df <- bacon_df %>% 
  mutate(
    SE_LLS = as.numeric(SE_lls),
    SE_KORA = as.numeric(SE_kora),
    SE_NTR = as.numeric(SE_ntr)
  )

SE_df <- pivot_longer(
  data = bacon_df,
  cols = c("SE_LLS", "SE_KORA", "SE_NTR"),
  names_to = "cohort",
  names_prefix = "SE_",
  values_to = "SE",
  values_drop_na = TRUE
)

SE_df %>% 
  ggplot(aes(x = cohort, y = SE, fill = cohort)) +
  scale_fill_manual(values = c("#e2455f", "orange", "#66cec8")) +
  geom_boxplot(alpha = 0.6) +
  ggtitle("SE distributions by cohort") +
  ylab("SE") + xlab("Cohort") +
  theme( 
    panel.background = element_rect(
      fill="white"),
    panel.border = element_rect(
      color="#1B2021", fill=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "#1B2021"),
    axis.text.x = element_text(size = 12, color = "#1B2021"),
    axis.title = element_text(size=11, color = "#1B2021"),
    plot.title = element_text(size = 16, hjust=0.5,
                              color = "#548687", face = "bold"))

##############################################################