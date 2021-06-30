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
#' @param x matrix of DNA methylation beta values. Can include missing values
#' @param y matrix of raw mRNA counts
#' @param clinical Data.frame conataing clinical covariates
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
#' data(BRCA_Methy_Test)
#' data(BRCA_RNA_Test)
#' data(BRCA_Clinical_Test)
#' dim(BRCA_Methy_Test)
#' BRCA_Methy_Test[1:5,1:5]
#' dim(BRCA_Clinical_Test)
#' head(BRCA_Clinical_Test)
#' dim(BRCA_RNA_Test)
#' BRCA_RNA_Test[1,1:5]
#' Construct(x = BRCA_Methy_Test, y = BRCA_RNA_Test,
#' clinical = BRCA_Clinical_Test, method = 'Interaction', beta = FALSE,
#' impute = FALSE)
#'


Construct <- function(x = betas, y = counts, clinical = clin , method =
                              c('Adjusted', 'Interaction', 'Specific'),
                          dist = NULL, nfolds = NULL, impute = TRUE,
                          beta = TRUE, parallel = FALSE, Methylation.Annotation = NULL,
                      Gene.Annotation = NULL){

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

    Methy <- as.matrix(x)

    annotated_RNA <- Gene.Annotation[match(rownames(y),
                                         Gene.Annotation$hgnc_symbol) ]


    Methylation.Annotation <- Methylation.Annotation[match(rownames(Methy),Methylation.Annotation$name) ]
    

    if(isTRUE(impute)){

    #####Impute missing methylation
    Methy <- imputation(Methy, dist, Methylation.Annotation)
    }
    if(isTRUE(beta)){
    Methy <- log2(Methy) - log2(1 - Methy)
    }
    
    

    if(isTRUE(parallel)){
    r <- foreach(i=seq_along(nrow(annotated_RNA)), .combine=rbind, .errorhandling='pass',
                 .packages = c('glmnet')) %dopar% {


        tmp <- annotated_RNA[i ]

        RNA_tmp <- Counts[i, ]

        methy_tmp <- Methy[(seqnames(Methylation.Annotation) == as.character(seqnames(tmp)))  & (end(Methylation.Annotation) <= (end(tmp) + dist)) & (start(Methylation.Annotation) >= (start(tmp) - dist))]

         clin_tmp <- Methylation.Annotation[(seqnames(Methylation.Annotation) == as.character(seqnames(tmp)))  & (end(Methylation.Annotation) <= (end(tmp) + dist)) & (start(Methylation.Annotation) >= (start(tmp) - dist))]
        ###Complete Case Analysis
        data.whole <- cbind(as.matrix(t(methy_tmp)),clinical)
        colnames(data.whole) <- c(colnames(as.matrix(t(methy_tmp))), colnames(clinical))
        data.whole <- data.whole[, colSums(is.na(data.whole)) == 0]
        ####Run Elastic Net
        out = switch(method,
                     "Adjusted"=Adjust(data.whole,RNA_tmp,0.5,nfolds,
                                       ncol(clinical)),
                     "Interaction"=Interact(data.whole,RNA_tmp,0.5,nfolds,
                                            ncol(clinical)),
                     "Specific"=Specific(data.whole,RNA_tmp,0.5,nfolds,
                                       ncol(clinical))
                     )
        out[['Annotation']] <- annotated_RNA[i,]
                 }
    return(out)
    } else {
        for(i in seq_len(nrow(annotated_RNA))){
            tmp <- annotated_RNA[i,]
            RNA_tmp <- Counts[,i]
            methy_tmp <- Methy[(seqnames(Methylation.Annotation) == as.character(seqnames(tmp)))  & (end(Methylation.Annotation) <= (end(tmp) + dist)) & (start(Methylation.Annotation) >= (start(tmp) - dist)),]
            clin_tmp <- methy_clin[(seqnames(Methylation.Annotation) == as.character(seqnames(tmp)))  & (end(Methylation.Annotation) <= (end(tmp) + dist)) & (start(Methylation.Annotation) >= (start(tmp) - dist))]
            ###Complete Case Analysis
            data.whole <- cbind(as.matrix(t(methy_tmp)),clinical)
            colnames(data.whole) <- c(colnames(as.matrix(t(methy_tmp))), colnames(clinical))
            data.whole <- data.whole[, colSums(is.na(data.whole)) == 0]
            out = switch(method,
                         "Adjusted"=Adjust(data.whole,RNA_tmp,0.5,nfolds,
                                           ncol(clinical)),
                         "Interaction"=Interact(data.whole,RNA_tmp,0.5,nfolds,
                                                ncol(clinical)),
                         "Specific"=Specific(data.whole,RNA_tmp,0.5,nfolds,
                                           ncol(clinical))
            )
            out[['Annotation']] <- annotated_RNA[i,]
        }
        return(out)
    }
}

