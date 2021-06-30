#' Function for associating race specific DNA methylation with matched mRNA expression
#'
#' This function allows one to associate DNA methylation with matched mRNA expression
#' in a multivariate context in race specific data.
#' @importFrom utils data
#' @import stats
<<<<<<< HEAD
#' @export Specific
=======
#' @export Specif
>>>>>>> upstream/master
#' @param x Processed DNA methylation with covariates. Race must be the final term
#' @param y Raw expression count vector
#' @param alpha Elastic net mixing penalty. Defaults to 0.5
#' @param nfolds the number of folds to use for cross-validation
#' @param clin.col the number of clinical covariates included
<<<<<<< HEAD
#' @return returns list of three elements for summarizing model construction. 
#' \code{Model} is the predicted model for the gene. \code{Predicted.Expression}
#' is the Cross-validated predicted expression. \code{Overall.D2} The prediction
#'  accuracy for the model. 
=======
#' @return returns list of five elements for summarizing model construction
#' @examples
>>>>>>> upstream/master
#' data.whole = cbind(matrix(rnorm(500),nrow=50), rbinom(50,1,0.5))
#' RNA_tmp <- rpois(50,15)
#' Specif(data.whole, RNA_tmp, alpha = 0.5, nfolds = 5, clin.col = 1)
#'



<<<<<<< HEAD
Specific <- function(x,y,alpha,nfolds,clin.col){
=======
Specif <- function(x,y,alpha,nfolds,clin.col){
>>>>>>> upstream/master
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package \"glmnet\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    tryCatch({
        whole.elastic.fit.cv <- glmnet::cv.glmnet(x,as.matrix(t(y)) , family ='poisson', nfolds = nfolds,alpha=alpha,penalty.factor = c(rep(1, ncol(x) - clin.col),0,0), parallel = FALSE )
        coef.min = stats::coef(whole.elastic.fit.cv, s = "lambda.min")
        active.min = which(as.numeric(coef.min) != 0)
        Whole.elastic = coef.min[active.min]
        names <- rownames(coef.min)
        a <- cbind(names[active.min],Whole.elastic)

        ct <- rep(NA,nfolds)
        b <- NULL
        tmp <- sample(seq(nfolds),nrow(x),replace=TRUE)

        for (i in seq(nfolds)){
<<<<<<< HEAD


            Whole.elastic.fit <- glmnet::cv.glmnet(x[!tmp %in% i,],t(y[!tmp %in% i]), family='poisson', nfolds = nfolds,alpha = alpha,penalty.factor = c(rep(1, ncol(x[!tmp %in% i,]) - clin.col),0,0))
            coef.min = stats::coef(whole.elastic.fit.cv, s = "lambda.min")

            out <- predict(Whole.elastic.fit, x[tmp %in% i,], s=Whole.elastic.fit$lambda.min, type = 'response', gamma=0.5)
            rownames(out) <- rownames(x[tmp %in% i,])
            out <- round(out)
            b <- rbind(b,out)

            y_i <- as.numeric(y[tmp %in% i])
=======
            xtrain <- x[!tmp %in% i,]
            xtest <- x[tmp %in% i,]
            ytrain <- y[!tmp %in% i]
            ytest <- y[tmp %in% i]


            Whole.elastic.fit <- glmnet::cv.glmnet(xtrain,t(ytrain), family='poisson', nfolds = nfolds,alpha = alpha,penalty.factor = c(rep(1, ncol(xtrain) - clin.col),0,0))
            coef.min = stats::coef(whole.elastic.fit.cv, s = "lambda.min")

            out <- predict(Whole.elastic.fit, xtest, s=Whole.elastic.fit$lambda.min, type = 'response', gamma=0.5)
            rownames(out) <- rownames(xtest)
            out <- round(out)
            b <- rbind(b,out)

            y_i <- as.numeric(ytest)
>>>>>>> upstream/master
            u_i <- out

            deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
            D <- 2 * sum(deviance.contribs)
            u_i <- rep(exp(coef.min[1]),length(y_i))
            deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
            Dn <- 2 * sum(deviance.contribs)

            ct[i] <- (1-D/Dn)

        }
<<<<<<< HEAD
        c <- ifelse(stats::median(ct, na.rm = TRUE) < 0,0,stats::median(ct, na.rm = TRUE))
    }, error=function(e){})
    output <- list('Model' = a, 'Predicted.Expression' = b, 'Overall.D2' = c)
=======
        c <- stats::median(ct, na.rm = TRUE)
        c <- ifelse(c < 0,0,c)
    }, error=function(e){})
    output <- list(a,b,c)
>>>>>>> upstream/master
    return(output)
}

