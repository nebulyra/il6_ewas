# Interleukin-6 EWAS Meta-analysis
This repository contains all scripts required to reproduce findings in our paper "Interleukin-6 and DNA methylation: links to cytokine signalling and metabolic reprogramming of immune cells".

The scripts are, in order:
1. Perform an EWAS on interleukin-6 levels within the LLS cohort
2. Perform an EWAS on C-reactive protein levels within the LLS cohort
3. Combine summary statistics from the three participating cohorts (LLS, NTR, and KORA)
4. Perform QC on the combined dataset
5. Adjust for bias and inflation in the test statistics using ```bacon```
6. Re-check quality of data following adjustment
7. Perform an EWAS meta-analysis in ```METAL```
8. Combine results from all meta-analysed models
9. Perform sensitivity analyses for extended cell counts and smoking status
10. Identify distinct genomic loci
11. Create a volcano plot of the results
12. Perform a test for enrichment of previous associations in the results
13. Investigate effects of adjusting IL-6 effects for CR-P and vice versa
14. Perform a test for enrichment of chromatin states in the results using a PBMC reference epigenome
15. Perform a test in ```HOMER``` for enrichment of transcription factor binding sites in the results
16. Identify genes within 100kb of the IL-6 associated CpGs
17. Extract all mQTLs within 1000kb of our CpGs from the full GoDMC data regardless of significance
18. Save all mQTLs within 1000kb of our CpGs
19. Perform SuSiE colocalization analyses
20. Save SuSiE colocalization results
21. Perform summary-based Mendelian randomization on the effects of DNAm on IL-6 and vice versa
