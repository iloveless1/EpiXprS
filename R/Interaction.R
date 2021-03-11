#' Function for associating race adjusted DNA methylation with matched mRNA expression
#'
#' This function allows one to associate DNA methylation with matched mRNA expression
#' in a multivariate context while adjusting for race. This function assumes
#' that the variable 'race' is the final column in the matrix.
#'
#' @importFrom utils data
#' @import stats
#' @export Interact
#' @param x Processed DNA methylation with covariates. Race must be the final term
#' @param y Raw expression count vector
#' @param alpha Elastic net mixing penalty. Defaults to 0.5
#' @param nfolds the number of folds to use for cross-validation
#' @param clin.col the number of clinical covariates included
#' @return returns list of five elements for summarizing model construction
#' @examples
#' data.whole = cbind(matrix(rnorm(500),nrow=50), rbinom(50,1,0.5))
#' RNA_tmp <- rpois(50,15)
#' Interact(data.whole, RNA_tmp, alpha = 0.5, nfolds = 10, clin.col = 2)
#'


Interact <- function(x,y,alpha,nfolds,clin.col){
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package \"glmnet\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    tryCatch({
        mat <- as.matrix(cbind(x[,seq(ncol(x)-clin.col)],x[,seq(ncol(x)-clin.col)]*x[,(ncol(x))],x[,(ncol(x)-(clin.col-1)):ncol(x)]))
        colnames(mat)[(ncol(x)-1):(ncol(x)*2-(2*clin.col))] <- paste0(colnames(mat)[(ncol(x)-(clin.col - 1)):(ncol(x)*2-(2*clin.col))],paste0(':',colnames(x)[ncol(x)]))
        whole.elastic.fit.cv <- glmnet::cv.glmnet(mat,as.matrix(t(y)) , family ='poisson', nfolds = nfolds,alpha=alpha,penalty.factor = c(rep(1, ncol(mat) - clin.col),0,0), parallel = FALSE )
        coef.min = stats::coef(whole.elastic.fit.cv, s = "lambda.min")
        active.min = which(as.numeric(coef.min) != 0)
        Whole.elastic = coef.min[active.min]
        names <- rownames(coef.min)
        a <- cbind(names[active.min],Whole.elastic)

        dt <- et <- ct <- rep(NA,nfolds)
        b <- NULL
        tmp <- sample(seq(nfolds),nrow(x),replace=TRUE)

        for (i in seq(nfolds)){
            xtrain <- mat[!tmp %in% i,]
            xtest <- mat[tmp %in% i,]
            ytrain <- y[!tmp %in% i]
            ytest <- y[tmp %in% i]


            Whole.elastic.fit <- glmnet::cv.glmnet(as.matrix(xtrain),as.matrix(t(ytrain)), family='poisson', nfolds = nfolds,alpha = alpha,penalty.factor = c(rep(1, ncol(xtrain) - clin.col),0,0))
            coef.min = stats::coef(whole.elastic.fit.cv, s = "lambda.min")

            out <- predict(Whole.elastic.fit, as.matrix(xtest), s=Whole.elastic.fit$lambda.min, type = 'response', gamma=alpha)
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

            y_i <- as.numeric(ytest)[xtest[,(ncol(xtest))] == 1]
            u_i <- out[xtest[,(ncol(xtest))] == 1]

            deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
            D <- 2 * sum(deviance.contribs)
            u_i <- rep(exp(coef.min[1]),length(y_i))
            deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
            Dn <- 2 * sum(deviance.contribs)

            dt[i] <- (1-D/Dn)

            y_i <- as.numeric(ytest)[xtest[,(ncol(xtest))] == 0]
            u_i <- out[xtest[,(ncol(xtest))] == 0]

            deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
            D <- 2 * sum(deviance.contribs)
            u_i <- rep(exp(coef.min[1]),length(y_i))
            deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
            Dn <- 2 * sum(deviance.contribs)

            et[i] <- (1-D/Dn)


        }
        c <- stats::median(ct, na.rm = TRUE)
        d <- stats::median(dt, na.rm = TRUE)
        f <- stats::median(et, na.rm = TRUE)

    }, error=function(e){})
    output <- list(a,b,c,d,f)
    return(output)
}
