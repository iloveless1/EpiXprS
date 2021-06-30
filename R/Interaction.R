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
<<<<<<< HEAD
#' @return returns list of five elements for summarizing model construction. 
#' \code{Model} is the predicted model for the gene. \code{Predicted.Expression}
#' is the Cross-validated predicted expression. \code{Overall.D2} The prediction
#'  accuracy for the interaction model. \code{Class1.D2} is the prediction accuracy
#'  for the reference class. \code{Class2.D2} is the prediction accuracy for the 
#'  non-reference class. 
#' 
#' @examples
#' data.whole = cbind(matrix(rnorm(500),nrow=50),rbinom(50,1,0.5) ,rbinom(50,1,0.5))
#' colnames(data.whole) <- c(1:10,'sex','race')
#' RNA_idx <- rpois(50,15)
#' Interact(data.whole, RNA_idx, alpha = 0.5, nfolds = 5, clin.col = 2)
=======
#' @return returns list of five elements for summarizing model construction
#' @examples
#' data.whole = cbind(matrix(rnorm(500),nrow=50),rbinom(50,1,0.5) ,rbinom(50,1,0.5))
#' colnames(data.whole) <- c(1:10,'sex','race')
#' RNA_tmp <- rpois(50,15)
#' Interact(data.whole, RNA_tmp, alpha = 0.5, nfolds = 5, clin.col = 2)
>>>>>>> upstream/master
#'


Interact <- function(x,y,alpha,nfolds,clin.col){
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package \"glmnet\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    tryCatch({
<<<<<<< HEAD
        int <- clin[,ncol(clin)]
        clin1 <- clin[,-c(ncol(clin))]
        mat <- model.matrix(y ~ clin1 + x*int)
        
=======
        mat <- as.matrix(cbind(x[,seq(ncol(x)-clin.col)],x[,seq(ncol(x)-clin.col)]*x[,(ncol(x))],x[,(ncol(x)-(clin.col-1)):ncol(x)]))
        colnames(mat)[(ncol(x)-(clin.col-1)):(ncol(x)*2-(2*clin.col))] <- paste0(colnames(mat)[(ncol(x)-(clin.col - 1)):(ncol(x)*2-(2*clin.col))],paste0(':',colnames(x)[ncol(x)]))
        colnames(mat)[c((ncol(mat)-(clin.col-1)):ncol(mat))] <- colnames(x)[c((ncol(x)-(clin.col-1)):ncol(x))]
>>>>>>> upstream/master
        whole.elastic.fit.cv <- glmnet::cv.glmnet(mat,as.matrix(t(y)) , family ='poisson', nfolds = nfolds,alpha=alpha,penalty.factor = c(rep(1, ncol(mat) - clin.col),0,0), parallel = FALSE )
        coef.min = stats::coef(whole.elastic.fit.cv, s = "lambda.min")
        active.min = which(as.numeric(coef.min) != 0)
        Whole.elastic = coef.min[active.min]
        names <- rownames(coef.min)
        a <- cbind(names[active.min],Whole.elastic)

        dt <- et <- ct <- rep(NA,nfolds)
        b <- NULL
<<<<<<< HEAD
        idx <- sample(seq(nfolds),nrow(x),replace=TRUE)

        for (i in seq(nfolds)){

            Whole.elastic.fit <- glmnet::cv.glmnet(as.matrix(mat[!idx %in% i,]),as.matrix(t(y[!idx %in% i])), family='poisson', nfolds = nfolds,alpha = alpha,penalty.factor = c(rep(1, ncol(mat[!idx %in% i,]) - clin.col),0,0))
            coef.min = stats::coef(whole.elastic.fit.cv, s = "lambda.min")

            out <- predict(Whole.elastic.fit, as.matrix(mat[idx %in% i,]), s=Whole.elastic.fit$lambda.min, type = 'response', gamma=alpha)
            rownames(out) <- rownames(mat[idx %in% i,])
            out <- round(out)
            b <- rbind(b,out)

            y_i <- as.numeric(y[idx %in% i])
=======
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
>>>>>>> upstream/master
            u_i <- out

            deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
            D <- 2 * sum(deviance.contribs)
            u_i <- rep(exp(coef.min[1]),length(y_i))
            deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
            Dn <- 2 * sum(deviance.contribs)

            ct[i] <- (1-D/Dn)

<<<<<<< HEAD
            y_i <- as.numeric(y[idx %in% i])[mat[idx %in% i,][,(ncol(mat[idx %in% i,]))] == 1]
            u_i <- out[mat[idx %in% i,][,(ncol(mat[idx %in% i,]))] == 1]
=======
            y_i <- as.numeric(ytest)[xtest[,(ncol(xtest))] == 1]
            u_i <- out[xtest[,(ncol(xtest))] == 1]
>>>>>>> upstream/master

            deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
            D <- 2 * sum(deviance.contribs)
            u_i <- rep(exp(coef.min[1]),length(y_i))
            deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
            Dn <- 2 * sum(deviance.contribs)

            dt[i] <- (1-D/Dn)

<<<<<<< HEAD
            y_i <- as.numeric(y[idx %in% i])[mat[idx %in% i,][,(ncol(mat[idx %in% i,]))] == 0]
            u_i <- out[mat[idx %in% i,][,(ncol(mat[idx %in% i,]))] == 0]
=======
            y_i <- as.numeric(ytest)[xtest[,(ncol(xtest))] == 0]
            u_i <- out[xtest[,(ncol(xtest))] == 0]
>>>>>>> upstream/master

            deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
            D <- 2 * sum(deviance.contribs)
            u_i <- rep(exp(coef.min[1]),length(y_i))
            deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
            Dn <- 2 * sum(deviance.contribs)

            et[i] <- (1-D/Dn)


        }
<<<<<<< HEAD
        c <- ifelse(stats::median(ct, na.rm = TRUE) < 0,0,stats::median(ct, na.rm = TRUE))
        d <- ifelse(stats::median(dt, na.rm = TRUE) < 0,0,stats::median(dt, na.rm = TRUE))
        f <- ifelse(stats::median(ft, na.rm = TRUE) < 0,0,stats::median(ft, na.rm = TRUE))

    }, error=function(e){})
    output <- list('Model' = a, 'Predicted.Expression' = b, 'Overall.D2' = c, 'EA.D2' = d, 'AA.D2' = f)
=======
        c <- stats::median(ct, na.rm = TRUE)
        d <- stats::median(dt, na.rm = TRUE)
        f <- stats::median(et, na.rm = TRUE)
        c <- ifelse(c < 0,0,c)
        d <- ifelse(d < 0,0,d)
        f <- ifelse(f < 0,0,f)

    }, error=function(e){})
    output <- list(a,b,c,d,f)
>>>>>>> upstream/master
    return(output)
}
