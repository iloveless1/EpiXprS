#' This function predicts mRNA expression
#'
#' Using models previously constructed using TCGA reference data, mRNA
#' expression is predicted using DNA methylation
#'
#' @importFrom utils data
#' @param Cancer The type of cancer model to use for expression prediction
#' @param x Matrix containing processed methylation beta-values
#' @param clinical Clinical data
#' @return Matrix of imputed expression values
#' @param impute Whether or not to impute the data. Defaults to TRUE
#' @param beta Whether methylation matrix is beta-values, defaults to TRUE
#' @param dist the distance from the TSS & TES in Kb. Defaults to 1,000,000Kb
#' @export EPI.Predict
#'
#' @examples
#' EPI.Predict(Cancer = 'PRAD', x = methy, clinical = clin)
#'
EPI_Predict <- function(Cancer = c('PRAD','BRCA','COAD','KIRP','KIRC','HNSC',
                                   'LUAD','UCEC'),x = methy, clinical = clin,
                                    dist = NULL, impute = TRUE, beta = TRUE){

    if(!colnames(clin) %in% c('race','age','ID'))
        stop("Clinical data colnames must be 'ID','race' and 'age'")

    if(is.na(colnames(methy) %in% clin$ID))
        stop('Clinical IDs must match methylation IDs must match')

    if(!colnames(clin) %in% c('age','race_list'))
        stop('age and race must be named "age" and "race_list"')

    EXP <- NULL
    gene <- NULL


    if(Cancer == 'PRAD'){
        load(url("https://github.com/iloveless1/EPIExPRS-Data/raw/main/PRAD.rda"))
        rt = PRAD
    }else if (Cancer == 'BRCA'){
        load(url("https://github.com/iloveless1/EPIExPRS-Data/raw/main/BRCA.rda"))
        rt = BRCA
    }else if (Cancer == 'COAD'){
        load(url("https://github.com/iloveless1/EPIExPRS-Data/raw/main/COAD.rda"))
        rt = COAD
    }else if (Cancer == 'KIRP'){
        load(url("https://github.com/iloveless1/EPIExPRS-Data/raw/main/KIRP.rda"))
        rt = KIRP
    }else if (Cancer == 'KIRC'){
        load(url("https://github.com/iloveless1/EPIExPRS-Data/raw/main/KIRC.rda"))
        rt = KIRC
    }else if (Cancer == 'HNSC'){
        load(url("https://github.com/iloveless1/EPIExPRS-Data/raw/main/HNSC.rda"))
        rt = HNSC
    }else if (Cancer == 'LUAD'){
        load(url("https://github.com/iloveless1/EPIExPRS-Data/raw/main/LUAD.rda"))
        rt = LUAD
    }else if (Cancer == 'UCEC'){
        load(url("https://github.com/iloveless1/EPIExPRS-Data/raw/main/UCEC.rda"))
        rt = UCEC
    } else{stop("Cancer must be one of ('PRAD','BRCA','COAD','KIRP','KIRC','HNSC',
              'LUAD','UCEC')")}


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

    mat <- as.matrix(rbind(Methy,Methy*clin$race_list,clin$race_list,clin$age))
    rownames(mat)[(nrow(Methy)+1):nrow(mat)] <- c(paste0(rownames(Methy),':Race')
                                                  ,'race_list','age')


    for(i in seq_along(nrow(rt))){
        model <- rt[[i,1]][[1]]
        predictors <- rbind(rep(1,ncol(mat)),mat[match(tmp1[2:nrow(model),1],rownames(mat)),])
        covariates <- as.numeric(model[,2])
        a <- exp(1)^colSums(predictors*covariates)
        EXP <- rbind(EXP,a)
        genes <- rt[[i,2]]
        gene <- rbind(gene,genes)

    }

    rownames(EXP) <- gene$ensemlb_id
    colnames(EXP) <- clin$ID
    return(as.matrix(EXP))
}
