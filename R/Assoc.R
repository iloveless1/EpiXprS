#' Function for associating predicted expression with outcome
#'
#' This function takes the predicted expression from EPI_Predict and tests for
#' associations with phenotype of interest.
#'
#'
#' @import stats
#' @import survival
#' @importFrom utils data
#' @param EXP Expression predicted using EPI_Predict
#' @param type Regression method for associating expression with a phenotype
#' @param clinical covariates to be included in the model. Phenotype to test must be in final column (final two for survival)
#' @return Data.frame with the association results for each gene
<<<<<<< HEAD:R/Assoc.R
#' @export Assoc
=======
#' @export methyl_Assoc
>>>>>>> upstream/master:R/methyl_Assoc.R
#' @examples
#' EXP <- matrix(rpois(500,15),nrow = 50)
#' rownames(EXP) <- c(1:50)
#' clin <- as.data.frame(rbinom(50,1,0.5))
<<<<<<< HEAD:R/Assoc.R
#' Assoc(EXP = EXP, type = 'logistic', clinical = clin)
#'
#'
Assoc <- function(EXP = EXP, type = c('logistic','survival','linear'), clinical = clin){
=======
#' methyl_Assoc(EXP = EXP, type = 'logistic', clinical = clin)
#'
#'
methyl_Assoc <- function(EXP = EXP, type = c('logistic','survival','linear'), clinical = clin){
>>>>>>> upstream/master:R/methyl_Assoc.R
    if (!requireNamespace("survival", quietly = TRUE)) {
        stop("Package \"survival\" needed for this function to work. Please install it.",
             call. = FALSE)
    }

        assoc_df <- NULL # Init association dataframe
        # Perform test between each pred_gene_exp column and phenotype
        EXP[is.infinite(EXP)] <- NA
        tryCatch({
        for (i in seq(nrow(EXP))) {
            pred_gene_exp <- EXP[i,]
            merged <- as.data.frame(cbind(pred_gene_exp,clinical))
            if (type == "logistic") {
                if((length(unique(merged[,ncol(merged)])) >2)){
                    stop("For model type logistic, phenotype must be binary")
                }
                results <- summary(stats::glm(as.formula(paste(colnames(merged)[ncol(merged)]," ~ .")), data = merged, family = binomial,na.action = na.omit))$coefficients['pred_gene_exp', ]
            } else if (type == "linear") {
                results <- summary(stats::lm(as.formula(paste(colnames(merged)[ncol(merged)]," ~ .")), data = merged,na.action = na.omit))$coefficients['pred_gene_exp', ]
            } else if (type == "survival") {
            results <- summary(survival::coxph(Surv(merged[,(ncol(merged)-1)],merged[,ncol(merged)]) ~ ., data = merged,na.action = na.omit))$coefficients['pred_gene_exp', ]
            }
            line <- c(rownames(EXP)[i],results)
            assoc_df <- rbind(assoc_df,line)
            colnames(assoc_df)[1] <- 'ensembl ID'
        }
}, error=function(e){})
        return(as.data.frame(assoc_df))
}
