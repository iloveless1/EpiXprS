#' Function for associating race specific DNA methylation with matched mRNA expression
#'
#' This function allows one to associate DNA methylation with matched mRNA expression
#' in a multivariate context in race specific data.
#' @importFrom utils data
#' @import stats
#' @param x Processed DNA methylation with covariates. Race must be the final term
#' @param y Raw expression count vector
#' @param alpha Elastic net mixing penalty. Defaults to 0.5
#' @param nfolds the number of folds to use for cross-validation
#' @param clin.col the number of clinical covariates included
#' @return returns list of five elements for summarizing model construction
#' @examples
#' \dontrun{
#' Specif(methy, counts, alpha = 0.5, nfolds = 10, clin.col = 3)
#' }



Specif <- function(x,y,alpha,nfolds,clin.col){
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package \"glmnet\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    tryCatch({
        whole.elastic.fit.cv <- glmnet::cv.glmnet(x,as.matrix(t(y)) , family ='poisson', nfolds = nfolds,alpha=alpha,penalty.factor = c(rep(1, ncol(xtrain) - clin.col),0,0), parallel = FALSE )
        coef.min = stats::coef(whole.elastic.fit.cv, s = "lambda.min")
        active.min = which(as.numeric(coef.min) != 0)
        Whole.elastic = coef.min[active.min]
        names <- rownames(coef.min)
        a <- cbind(names[active.min],Whole.elastic)

        ct <- rep(NA,nfolds)
        b <- NULL
        tmp <- sample(seq(nfolds),nrow(x),replace=TRUE)

        for (i in seq(nfolds)){
            xtrain <- mat[!tmp %in% i,]
            xtest <- mat[tmp %in% i,]
            ytrain <- y[!tmp %in% i]
            ytest <- y[tmp %in% i]


            Whole.elastic.fit <- glmnet::cv.glmnet(xtrain,t(ytrain), family='poisson', nfolds = nfolds,alpha = alpha,penalty.factor = c(rep(1, ncol(xtrain) - clin.col),0,0))
            coef.min = stats::coef(whole.elastic.fit.cv, s = "lambda.min")

            out <- predict(Whole.elastic.fit, xtest, s=Whole.elastic.fit$lambda.min, type = 'response', gamma=0.5)
            rownames(out) <- rownames(xtest)
            out <- round(out)
            b <- rbind(b,out)

            y_i <- as.numeric(ytest)
            u_i <- out

            deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
            D <- 2 * sum(deviance.contribs)
            u_i <- rep(exp(coef.min[1]),length(y_i))
            deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
            Dn <- 2 * sum(deviance.contribs)

            ct[i] <- (1-D/Dn)

        }
        c <- stats::median(ct, na.rm = TRUE)
    }, error=function(e){})
    output <- list(a,b,c)
    return(output)
}

