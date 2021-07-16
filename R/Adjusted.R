#' Function for associating race adjusted DNA methylation with matched mRNA expression
#'
#' This function allows one to associate DNA methylation with matched mRNA expression
#' in a multivariate context while adjusting for race. This function assumes
#' that the variable 'race' is the final column in the matrix.
#'
#'
#' @export Adjust
#' @import stats
#' @importFrom utils data
#' @param x Processed DNA methylation with covariates. Race must be the final term
#' @param y Raw expression count vector
#' @param alpha Elastic net mixing penalty. Defaults to 0.5
#' @param nfolds the number of folds to use for cross-validation
#' @param clin a matrix of clinical covariates
#' @return returns list of five elements for summarizing model construction. 
#' \code{Model} is the predicted model for the gene. \code{Predicted.Expression}
#' is the Cross-validated predicted expression. \code{Overall.D2} The prediction
#'  accuracy for the adjusted model. \code{Class1.D2} is the prediction accuracy
#'  for the reference class. \code{Class2.D2} is the prediction accuracy for the 
#'  non-reference class. 
#' @examples
#' data.whole = cbind(matrix(rnorm(500),nrow=50), rbinom(50,1,0.5))
#' RNA_tmp <- rpois(50,15)
#' Adjust(data.whole, RNA_tmp, alpha = 0.5, nfolds = 5, clin = clinical)
#'



Adjust <- function(x,y,alpha,nfolds,clin){

    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package \"glmnet\" needed for this function to work. Please install it.",
        call. = FALSE)
    }

    tryCatch({

            mat <- model.matrix(y ~ x + clin1 -1)
     
        whole.elastic.fit.cv <- glmnet::cv.glmnet(mat,as.matrix(t(y)), family ='poisson', nfolds = nfolds,alpha=alpha,penalty.factor = c(rep(1, ncol(mat) - ncol(clin)),0,0), parallel = FALSE )
        coef.min = stats::coef(whole.elastic.fit.cv, s = "lambda.min")
        active.min = which(as.numeric(coef.min) != 0)
        Whole.elastic = coef.min[active.min]
        names <- rownames(coef.min)
        a <- data.frame(row.names = names[active.min], Weights = Whole.elastic)
        
        dt <- et <- ct <- rep(NA,nfolds)
        b <- NULL
        tmp <- sample(seq(nfolds),nrow(x),replace=TRUE)

        for (i in seq(nfolds)){
            
            Whole.elastic.fit <- glmnet::cv.glmnet(mat[!tmp %in% i,],t(y[!tmp %in% i]),
            family='poisson', nfolds = nfolds,alpha = alpha,penalty.factor =
            c(rep(1, ncol(mat) - 2),0,0))
            coef.min = stats::coef(whole.elastic.fit.cv, s = "lambda.min")
            out <- predict(Whole.elastic.fit, mat[tmp %in% i,],
            s=Whole.elastic.fit$lambda.min, type = 'response', gamma=0.5)
            rownames(out) <- rownames(mat[tmp %in% i,])
            out <- round(out)
            b <- rbind(b,out)
            y_i <- as.numeric(y[tmp %in% i])
            u_i <- out
            
            ct[i] <- D2(y_i,u_i)
            
            y_i <- as.numeric(y[tmp %in% i])[x[tmp %in% i,][,(ncol(x[tmp %in% i,]))] == 1]
            u_i <- out[x[tmp %in% i,][,(ncol(x[tmp %in% i,]))] == 1]
            
            dt[i] <- D2(y_i,u_i)
            y_i <- as.numeric(y[tmp %in% i])[x[tmp %in% i,][,(ncol(x[tmp %in% i,]))] == 0]
            u_i <- out[x[tmp %in% i,][,(ncol(x[tmp %in% i,]))] == 0]
            
            et[i] <- D2(y_i,u_i)
        }
        c <- ifelse(stats::median(ct, na.rm = TRUE) < 0,0,stats::median(ct, na.rm = TRUE))
        d <- ifelse(stats::median(dt, na.rm = TRUE) < 0,0,stats::median(dt, na.rm = TRUE))
        f <- ifelse(stats::median(et, na.rm = TRUE) < 0,0,stats::median(et, na.rm = TRUE))

    }, error=function(e){})
    output <- list('Model' = a, 'Predicted.Expression' = b, 'Overall.D2' = c, 'Class1.D2' = d, 'Class2.D2' = f)
    return(output)
}

