#' Function for associating predicted expression with outcome
#'
#' This function takes the predicted expression from EPI_Predict and tests for
#' associations with phenotype of interest.
#'
#'
#' @import stats
#' @import survival
#' @param EXP Expression predicted using EPI_Predict
#' @param type Regression method for associating expression with a phenotype
#' @param clinical covariates to be included in the model. Phenotype to test must be in final column (final two for survival)
#' @return Data.frame with the association results for each gene
EPI.Assoc <- function(EXP = EXP, type = c('logistic','survival','linear'), clinical = clin){

    if (!requireNamespace("survival", quietly = TRUE)) {
        stop("Package \"survival\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    if(type == 'logistic' & unique(clin[,ncol(clin)]) >2){
        stop("For model type logistic, phenotype must be binary")
    }
        assoc_df <- NULL # Init association dataframe
        # Perform test between each pred_gene_exp column and phenotype
        for (i in seq(nrow(EXP))) {
            pred_gene_exp <- t(EXP[i,])
            merged <- as.data.frame(cbind(clinical, pred_gene_exp))
            if (type == "logistic") {
                results <- summary(stats::glm(colnames(merged)[ncol(merged)] ~. , data = merged, family = binomial))$coefficients['pred_gene_exp', ]
            } else if (type == "linear") {
                results <- summary(stats::lm(colnames(merged)[ncol(merged)] ~., data = merged))$coefficients['pred_gene_exp', ]
            } else if (type == "survival") {
            results <- summary(survival::coxph(Surv(c(merged[(ncol(merged)-1)],merged[ncol(merged)])) ~ ., data = merged))$coefficients['pred_gene_exp', ]
            }
            line <- c(rowname(EXP)[i],results)
            assoc_df <- rbind(assoc_df,line)
        }

        return(as.data.frame(assoc_df))
}
