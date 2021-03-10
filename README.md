# EPIExPRS
Methods for association DNA methylation with gene expression in a multivariate context
---
title: "EPIExPRS: A poly-CpG site approach to gene expression prediction and association"
author: "Ian Loveless"


DNA methylation has been shown to be in important regulator of gene expression. 
This relationship is complex, with multiple CpG sites in a 
region influencing expression rather than by a single CpG site in isolation. The 'EPIExPRS' package provides a means
for constructing these multi-CpG models, allowing for the further complexity of context dependent modifiers of the relationship between methylation and gene expression. In particular, EPIExPRS was developed to quantify and account for the context dependent effect of race through the inclusion of CpG-by-race interaction terms. Namely, to construct the original models, paired DNA methylation and gene expression were downloaded from The Cancer Genome Atlas TCGA, using the TCGbiolinks package, to assess the role of race as a modifier the multivariate relationship between DNA methylation and gene expression. Eight tumor types were downloaded that had high African American representation to allow us to assess these potential differences. While the subsequent description therefore focuses on race (African American and European American) the context being evaluated, the contribution of any potential dichotomous modifier (e.g. environmental exposure) could be assessed. 

The main functionality of this package is to allow one to associate many CpG sites with gene expression, while accounting for possible race-specific effects. This is managed by including multiplicative interation terms between race and CpG sites in a penalized regression model. The function 'EPI.Construct()' 
has three methods, 'Specific' which constructs models that ignore possible race-specific effects, 'Adjusted' where race is included in the model, and 'Interaction' where the aforementioned multiplicative interaction terms are included as well. The function returns 1) the prediction model, 2) the cross-validated predicted expression, 3) the percent of variability explained, 4) & 5) the race specific variability
explained for the 'Adjusted' and 'Interaction' model types, and 6) ensembl annotation for the genes modeled. 

Another useful function, 'model_search()', allows one to search through pre-constructed models to identify gene models in which a given CpG is included.

Similarly, pre-constructed models can be used to predict gene expression from a  450K DNA methylation data set, using the 'EPI.Predict()' function. Prior to use, we advise missing DNA methylation be imputed using the 'Imputation()' function. This function is borrowed from the 'methyLImp' package, toprovide complete predicted expression. For more details, see ?EPI.Predict. 

Finally, predicted expression scores from the EPI.Predict() function can be used to associate predicted expression with a phenotype of interest using the EPI.Assoc fuction. This function currently implements linear regression, logistic regression, and Cox proportional hazards models. For more details, see ?EPI.Assoc
