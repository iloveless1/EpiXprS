#' Function for associating DNA methylation with matched mRNA expression
#'
#' This function allows you to associate DNA methylation with matched mRNA
#' expression. There are three methods to choose from: 'Adjusted' in which
#' one can adjust for race, 'Interaction' in which one can enforce a weak hierarchy
#' and include race interaction terms, and 'Specific' in which one can associate
#' in a race-specific manner.
#'
#' @importFrom utils data
#' @import glmnet
#' @import foreach
#' @import Repitools
#' @import doParallel
#' @import AnnotationHub
#' @import ExperimentHub
#' @import MultiAssayExperiment
#' @import SummarizedExperiment
#' @param Object an object of class Summarized Experiment containing paired 
#' DNA methylation, RNA seq, and Clinical data. 
#' @param parallel whether to run in parallel. Defaults to FALSE
#' @param method method for multivariate associations
#' @param dist the distance from the TSS & TES in Kb. Defaults to 1,000,000Kb
#' @param nfolds the number of folds to use for cross-validation
#' @param impute Whether or not to impute the data. Defaults to TRUE
#' @param beta Whether methylation matrix is beta-values, defaults to TRUE
#' @param Methylation.Annotation A Granges object annotating the methylation array used
#' @param Gene.Annotation A GRanges object containing gene annotation. 
#' @return Large list with up to 5 columns
#'
#' @export Construct
#'
#' @examples
#' data('Testdat',package = 'EpiXprS')
#' Methy <- assays(Testdat)$Methyl450k
#' Methy[1:5,1:5]
#' Counts <- assays(Testdat)$RNASeqGene
#' Counts[1:2,1:5]
#' coldat <- colData(Testdat)
#' coldat$race <- ifelse(coldat$race == 'white',1,0)
#' head(coldat)
#' Construct(x = Methy, y = Counts, clinical = coldat$race,
#' method ='Interaction', dist = NULL, nfolds = NULL, 
#' impute = TRUE, beta = TRUE, parallel = FALSE,
#' Methylation.Annotation = annotation_450k, Gene.Annotation = annotated_RNA)


Construct <- function(Object = Object, method =
                              c('Adjusted', 'Interaction', 'Specific'),
                          dist = NULL, nfolds = NULL, impute = TRUE,
                          beta = TRUE, parallel = FALSE, Methylation.Annotation = NULL,
                      Gene.Annotation = NULL){
    if(class(Object) != "MultiAssayExperiment"){
        stop("Object must be MultiAssayExperiment. Please use 'MakeEpiXprs'")
    }
    
    method = match.arg(method)
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package \"glmnet\" needed for this function to work. Please
             install it.",
             call. = FALSE)
    }

    if(is.null(Methylation.Annotation)) {
        stop("Annotation must be provided for methylation array")
        }
    if(is.null(Gene.Annotation)){
        stop("Annotation must be provided for RNA counts")
    }
    if(is.null(dist)){
        dist <- 1000000
    }

    if(is.null(nfolds)){
        nfolds = 10
    }

    Methy <- as.matrix(assays(Object)$Methyl)
    y <- as.matrix(assays(Object)$RNASeqGene)
    annotated_RNA <- Gene.Annotation[match(rownames(),
                                         Gene.Annotation$hgnc_symbol) ]


    Methylation.Annotation <- Methylation.Annotation[match(rownames(Methy),Methylation.Annotation$name) ]
    

    if(isTRUE(impute)){

    #####Impute missing methylation
    Methy <- imputation(Methy, dist, Methylation.Annotation)
    }
    if(isTRUE(beta)){
    Methy <- log2(Methy) - log2(1 - Methy)
    }
    
    ###Run twice because imputation removes some CpG sites
    Methylation.Annotation <- Methylation.Annotation[match(rownames(Methy),Methylation.Annotation$name) ]
    

    if(isTRUE(parallel)){
    r <- foreach(i=seq_along(length(annotated_RNA)), .combine=rbind, .errorhandling='pass',
                 .packages = c('glmnet')) %dopar% {

                     input <- prepare(y, Methy, Methylation.Annotation,annotated_RNA, idx = i,dist)
                     out = switch(method,
                                  "Adjusted"=Adjust(input$methy_tmp,input$RNA_tmp,0.5,nfolds,
                                                    clin),
                                  "Interaction"=Interact(input$methy_tmp,input$RNA_tmp,0.5,nfolds,
                                                         clin),
                                  "Specific"=Specific(input$methy_tmp,input$RNA_tmp,0.5,nfolds,
                                                      clin)
                     )
        out[['Annotation']] <- annotated_RNA[i,]
                 }
    return(out)
    } else {
        for(i in seq_len(length(annotated_RNA))){
            
            input <- prepare(y, Methy, Methylation.Annotation,annotated_RNA, idx = i, dist)
            out = switch(method,
                         "Adjusted"=Adjust(input$methy_tmp,input$RNA_tmp,0.5,nfolds,
                                           clin),
                         "Interaction"=Interact(input$methy_tmp,input$RNA_tmp,0.5,nfolds,
                                                clin),
                         "Specific"=Specific(input$methy_tmp,input$RNA_tmp,0.5,nfolds,
                                           clin)
            )
            out[['Annotation']] <- annotated_RNA[i,]
        }
        return(out)
    }
}

