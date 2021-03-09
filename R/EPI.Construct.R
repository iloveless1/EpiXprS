#' Function for associating DNA methylation with matched mRNA expression
#'
#' This function allows you to associate DNA methylation with matched mRNA
#' expression. There are three methods to choose from: 'Adjusted' in which
#' one can adjust for race, 'Interaction' in which one can enforce a weak hierarchy
#' and include race interaction terms, and 'Specific' in which one can associate
#' in a race-specific manner.
#'
#'
#' @import glmnet
#' @import foreach
#' @import Repitools
#' @import doParallel
#' @param x matrix of DNA methylation beta values. Can include missing values
#' @param y matrix of raw mRNA counts
#' @param clinical Data.frame conataing clinical covariates
#' @param parallel whether to run in parallel. Defaults to FALSE
#' @param method method for multivariate associations
#' @param dist the distance from the TSS & TES in Kb. Defaults to 1,000,000Kb
#' @param nfolds the number of folds to use for cross-validation
#' @param impute Whether or not to impute the data. Defaults to TRUE
#' @param beta Whether methylation matrix is beta-values, defaults to TRUE
#' @param array Methylation array type, 450k or EPIC, defaults to 450K
#' @return Large list with upto 5 columns


EPI.Construct <- function(x = betas, y = counts, clinical = clin , method =
                              c('Adjusted', 'Interaction', 'Specific'),
                          dist = NULL, nfolds = NULL, impute = TRUE,
                          beta = TRUE, parallel = FALSE, array = '450K'){

    method = match.arg(method)
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package \"glmnet\" needed for this function to work. Please
             install it.",
             call. = FALSE)
    }


    if(ncol(x) != ncol(y))
        stop('Number of samples in betas and counts must match')

    if(ncol(x) != nrow(clinical))
        stop('Number of samples in betas and clin must match')

    data("annotated_RNA")
    if(array == "EPIC") {
        methy_clin <- readRDS(url('https://zwdzwd.s3.amazonaws.com/
                                          InfiniumAnnotation/20180909/EPIC/
                                          EPIC.hg38.manifest.rds'))
        } else {
        methy_clin <- readRDS(url('https://zwdzwd.s3.amazonaws.com/
                                          InfiniumAnnotation/20180909/HM450/
                                          HM450.hg38.manifest.rds'))
            }
    methy_clin <- annoGR2DF(methy_clin)
    methy_clin <- methy_clin[rownames(methy_clin) %in% rownames(x),
                             c("chr", "start", "gene_HGNC")]
    if(is.null(dist)){
        dist <- 1000000
    }

    if(is.null(nfolds)){
        nfolds = 10
    }

    Methy <- as.matrix(x)


    counts <- y[rowSums(y==0)/ncol(y) < 0.2, ]
    annotated_RNA <- annotated_RNA[match(rownames(counts),
                                         annotated_RNA$ensembl_id), ]

    methy_clin <- methy_clin[match(rownames(Methy),rownames(methy_clin)), ]


    ####set-up so I can remove this step

    methy_clin$chr <- ifelse(methy_clin$chr == 'X',23,
                             ifelse(methy_clin$chr == 'Y', 24, methy_clin$chr))
    annotated_RNA$chromosome_name <- ifelse(annotated_RNA$chromosome_name == 'X',
            23, ifelse(annotated_RNA$chromosome_name == 'Y', 24,
                       annotated_RNA$chromosome_name))



    if(isTRUE(impute)){

    #####Impute missing methylation
    missing <- rowSums(is.na(methy))

    for(i in seq_len(length(missing))){
        tryCatch({
            if(missing[i]==0) next()

            tmp <- methy_clin[i, ]

            impute <- Methy[(methy_clin$chr == tmp$Chromosome & methy_clin$start
                             <= tmp$Genomic_Coordinate+dist) &
                                (methy_clin$chr == tmp$Chromosome &
                                     methy_clin$start >=
                                     tmp$Genomic_Coordinate-dist), ]

            impute <- impute[rowSums(is.na(impute)) == 0 | rownames(impute)
                             %in% rownames(tmp), ]

            out <- Imputation(t(impute))
            Methy[i,] <- out[ ,match(rownames(tmp),colnames(out))]
        }, error=function(e){})
    }
    Methy <- ifelse(Methy == 1, 0.99, ifelse(Methy == 0, 0.01, Methy))
    }

    if(isTRUE(beta)){
    Methy <- log2(Methy) - log2(1 - Methy)
    }

    if(isTRUE(parallel)){
        Whole <- matrix(list(),nrow(counts),6)

    r <- foreach(i=seq_along(nrow(annotated_RNA)), .combine=rbind, .errorhandling='pass',
                 .packages = c('glmnet')) %dopar% {


        tmp <- annotated_RNA[i, ]

        RNA_tmp <- counts[i, ]

        methy_tmp <- Methy[(methy_clin$chr == tmp$chromosome_name &
                                methy_clin$start <= tmp$end_position+dist) &
                               (methy_clin$chr == tmp$chromosome_name &
                                    methy_clin$start >=
                                    tmp$start_position-dist), ]

         clin_tmp <- methy_clin[(methy_clin$chr == tmp$chromosome_name &
                                     methy_clin$start <= tmp$end_position+dist)
                                & (methy_clin$chr == tmp$chromosome_name &
                                       methy_clin$start
                                   >= tmp$start_position-dist), ]
        ###Complete Case Analysis
        data.whole <- cbind(as.matrix(t(methy_tmp)),clinical)
        colnames(data.whole) <- c(colnames(as.matrix(t(methy_tmp))),
                                  colnames(clinical))
        data.whole <- data.whole[, colSums(is.na(data.whole)) == 0]
        ####Run Elastic Net
        out = switch(method,
                     "Adjusted"=Adjust(data.whole,RNA_tmp,0.5,nfolds,
                                       ncol(clinical)),
                     "Interaction"=Interact(data.whole,RNA_tmp,0.5,nfolds,
                                            ncol(clinical)),
                     "Specific"=Specif(data.whole,RNA_tmp,0.5,nfolds,
                                       ncol(clinical))
                     )
        Whole[[1]] <- ifelse(length(out[1])==0,'EMPTY',out[1])
        Whole[[2]] <- ifelse(length(out[2])==0,'EMPTY',out[2])
        Whole[[3]] <- ifelse(length(out[4])==0,'EMPTY',out[4]) ##### EA D2
        Whole[[4]] <- ifelse(length(out[3])==0,'EMPTY',out[3]) ##### Complete D2
        Whole[[5]] <- ifelse(length(out[5])==0,'EMPTY',out[5]) ##### AA D2
        Whole[[6]] <- annotated_RNA[i,]
        return(out)
    }
    return(r)
    } else {
        Whole <- matrix(list(),nrow(counts),6)
        for(i in seq_len(nrow(annotated_RNA))){
            tmp <- annotated_RNA[i,]
            RNA_tmp <- counts[i,]
            methy_tmp <- Methy[(methy_clin$chr == tmp$chromosome_name
                                & methy_clin$start <= tmp$end_position+dist)
                               & (methy_clin$chr == tmp$chromosome_name
                            & methy_clin$start >= tmp$start_position-dist),]
            clin_tmp <- methy_clin[(methy_clin$chr == tmp$chromosome_name
                                    & methy_clin$start <= tmp$end_position+dist)
                                   & (methy_clin$chr == tmp$chromosome_name
                                & methy_clin$start >= tmp$start_position-dist),]
            ###Complete Case Analysis
            data.whole <- cbind(as.matrix(t(methy_tmp)),clinical)
            colnames(data.whole) <- c(colnames(as.matrix(t(methy_tmp))),
                                      colnames(clinical))
            data.whole <- data.whole[, colSums(is.na(data.whole)) == 0]
            out = switch(method,
                         "Adjusted"=Adjust(data.whole,RNA_tmp,0.5,nfolds,
                                           ncol(clinical)),
                         "Interaction"=Interact(data.whole,RNA_tmp,0.5,nfolds,
                                                ncol(clinical)),
                         "Specific"=Specif(data.whole,RNA_tmp,0.5,nfolds,
                                           ncol(clinical))
            )
            Whole[[i,1]] <- ifelse(length(out[1])==0, 'EMPTY', out[1])
            Whole[[i,2]] <- ifelse(length(out[2])==0, 'EMPTY', out[2])
            Whole[[i,3]] <- ifelse(length(out[4])==0, 'EMPTY', out[4])
            Whole[[i,4]] <- ifelse(length(out[3])==0, 'EMPTY', out[3])
            Whole[[i,5]] <- ifelse(length(out[5])==0, 'EMPTY', out[5])
            Whole[[i,6]] <- annotated_RNA[i, ]
        }
        return(Whole)
    }
}

