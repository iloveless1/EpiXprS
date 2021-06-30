#' This function predicts mRNA expression
#'
#' Using models previously constructed using TCGA reference data, mRNA
#' expression is predicted using DNA methylation
#'
#' @importFrom utils data
#' @param Cancer ExperimentHub object cancer model
#' @param Methy Matrix containing processed methylation beta-values
#' @param clinical Clinical data
#' @return Matrix of imputed expression values
#' @param impute Whether or not to impute the data. Defaults to TRUE
#' @param beta Whether methylation matrix is beta-values, defaults to TRUE
#' @param dist the distance from the TSS & TES in Kb. Defaults to 1,000,000Kb
#' @export Predict
#'
#' @examples
#' data(BRCA_Methy_Test)
#' data(BRCA_Clinical_Test)
#' colnames(BRCA_Clinical_Test) <- c('age','race_list')
#' BRCA_Clinical_Test <- as.data.frame(BRCA_Clinical_Test)
#' Predict(Cancer = 'BRCA', Methy = BRCA_Methy_Test, clinical = BRCA_Clinical_Test,
#' impute = FALSE, beta = FALSE)
#'
#'
Predict <- function(Cancer = object, Methy = methy, clinical = clin,
                                    dist = NULL, impute = TRUE, beta = TRUE){

    if(sum(!c('age','race_list') %in% colnames(clinical))>0)
        stop('age and race must be named "age" and "race_list"')

    EXP <- NULL
    gene <- NULL

    if(is.null(dist)){
        dist <- 1000000
    }

    if(isTRUE(impute)){

        #####Impute missing methylation
        Methy <- imputation(Methy, dist, methy_clin)
    }
    if(isTRUE(beta)){
        Methy <- log2(Methy) - log2(1 - Methy)
    }
    mat <- as.matrix(rbind(Methy,Methy*clinical$race_list,clinical$race_list,clinical$age))
    rownames(mat)[(nrow(Methy)+1):nrow(mat)] <- c(paste0(rownames(Methy),':Race')
                                                  ,'race_list','age')


    for(i in seq_along(nrow(rt))){
        model <- Cancer[[i,1]][[1]]
        predictors <- rbind(rep(1,ncol(mat)),mat[match(model[2:nrow(model),1],rownames(mat)),])
        covariates <- as.numeric(model[,2])
        a <- exp(1)^colSums(predictors*covariates)
        EXP <- rbind(EXP,a)
        genes <- rt[[i,2]]
        gene <- rbind(gene,genes)

    }

    rownames(EXP) <- gene$ensemlb_id
    colnames(EXP) <- colnames(mat)
    return(as.matrix(EXP))
}
