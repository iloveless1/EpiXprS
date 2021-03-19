# EPIExPRS
DNA methylation has been shown to be in important regulator of gene expression. This relationship is complex, with multiple CpG sites in a region influencing expression rather than by a single CpG site in isolation. The 'EPIExPRS' package provides a means for constructing these multi-CpG models, allowing for the further complexity of context dependent modifiers of the relationship between methylation and gene expression. In particular, EPIExPRS was developed to quantify and account for the context dependent effect of race through the inclusion of CpG-by-race interaction terms. Namely, to construct the original models, paired DNA methylation and gene expression were downloaded from The Cancer Genome Atlas TCGA, using the TCGAbiolinks package, to assess the role of race as a modifier in the multivariate relationship between DNA methylation and gene expression. Eight tumor types were downloaded that had high African American representation to allow us to assess these potential differences. While the subsequent description therefore focuses on race (African American and European American) the context being evaluated, the contribution of any potential dichotomous modifier (e.g. environmental exposure) could be assessed.

# Installation
if (!requireNamespace("BiocManager", quietly = TRUE))  
    install.packages("BiocManager")  
BiocManager::install("EPIExPRS")  
