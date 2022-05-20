## The implementation of the prediction method presented in the study by Franco et al. 2013

In an analysis of two cohorts, Franco et al. 2013 [58] found 301 transcripts for which the gene expression was significantly correlated with the magnitude of the antibody response to the influenza vaccine. Of the 301 transcripts, 49 of them had max correlation on day 0 with the antibody response. These transcripts were considered a transcriptional signature that had the potential to predict the vaccine response and the predictive power of this gene set was therefore evaluated here. Here there was also no expected directionality of the gene signature, as this was not evaluated in the original study. In the study, the 49 Illumina probe IDs were annotated with 42 unique gene names, 37 of which could be identified with an Ensembl ID in the CMI-PB gene expression datasets. It was also considered to map directly between Illumina probe IDs and Ensembl IDs using the R package biomaRt. In the original study, both Illumina HumanHT-12v3 andHumanHT-12v4 Expression BeadChips were used. However, when mapping directly between the probe and Ensembl IDs, only 28 of the gene symbols were identical to the 42 annotated gene symbols in the study. For this reason, it was determined that using the annotated gene symbols to identify genes in the CMI-PB gene expression datasets was the best option. Of the 37 identified genes, 22 of them were positively correlated with the antibody response and therefore considered upregulated genes in the signature, while the 15 negatively correlated were considered downregulated. The original probe IDs, the identified Ensembl IDs as well as gene symbols are presented in Supplementary Table S9. Preprocessing in the original study included background adjustment and Illumina microarray-specific methods for variance stabilizing transformation [59] and robust spline normalization using the lumi R package [60]. The TPM normalization method was used for the CMI-PB RNA-seq datasets, as the custom gene signature calculation presented in section X.X Calculation of gene signature scores was used.