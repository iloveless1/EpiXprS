#' Prepare CpG sites and expression for modeling
#' 
#' @param Counts A matrix of RNA counts
#' @param Methy A matrix of Methylation M-values
#' @param Methylation.Annotation A GRanges object containing CpG site locations.
#' must have same order as \code{Methy}
#' @param annotated_RNA A GRanges object containing RNA annotation. Must have 
#' same order as \code{Counts}
#' @return A list containing the prepared methylation and RNA
#' 
prepare <- function(Counts = Counts, Methy = Methy, Methylation.Annotation = Methylation.Annotation, annotated_RNA = annotated_RNA, idx = idx, dist){
    tmp <- annotated_RNA[idx,]
    RNA_tmp <- Counts[idx,]
    methy_tmp <- Methy[(as.character(seqnames(Methylation.Annotation)) == as.character(seqnames(tmp)))  & (GenomicRanges::end(Methylation.Annotation) <= (as.numeric(GenomicRanges::end(tmp)) + as.numeric(dist))) & (GenomicRanges::start(Methylation.Annotation) >= (GenomicRanges::start(tmp) - dist)),]
    ###Complete Case Analysis
    methy_tmp <- t(methy_tmp)
    methy_tmp <- methy_tmp[, colSums(is.na(methy_tmp)) == 0]
    
    return(list(RNA_tmp = RNA_tmp,methy_tmp = methy_tmp))
}


#' Calculate Accuracy
#' 
#' @param y_i Observed counts
#' @param u_i Predicted counts
#' @return the McFadden's Pseudo R-squared
#' 
#' @examples 
#' y_i <- rpois(50,12)
#' u_i <- rpois(50,11.5)
#' D2(y_i,u_i)

D2 <- function(y_i,u_i){
deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
D <- 2 * sum(deviance.contribs)
u_i <- rep(exp(coef.min[1]),length(y_i))
deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
Dn <- 2 * sum(deviance.contribs)

return(1-D/Dn)
}







