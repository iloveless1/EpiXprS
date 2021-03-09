#' annotated_RNA.
#'
#' Annotation for hg38 aligned mRNA.
#'
#' @format A data frame with five variables:
#' \describe{
#' \item{\code{ensembl_id}}{Ensebl Gene ID},
#' \item{\code{chromosome_name}}{chromosome},
#' \item{\code{start_position}}{Transcription Start Site},
#' \item{\code{end_position}}{Transcription Stop Site},
#' \item{\code{hgnc_symbol}}{Gene ID}
#'}
#'
"annotated_RNA"

#' BRCA_Methy_Test.
#'
#' Test data for vignetta from TCGA-BRCA Methylation.
#'
#' @format A matrix with 468 sample and 258 CpG site M-values:
#' \describe{
#' \item{\code{BRCA_Methy_Test}}{Ensebl Gene ID},
#'}
#'
"BRCA_Methy_Test"

#' BRCA_Clinical_Test.
#'
#' Clinical Covariates for BRCA Samples
#'
#' @format A data frame with two variables:
#' \describe{
#' \item{\code{age}}{Age for samples},
#' \item{\code{race}}{Race of samples},
#'}
#'
"BRCA_Clinical_Test"

#' BRCA_RNA_Test.
#'
#' A data.frame with 468 samples a expression counts for TSPAN6:
#'
#' @format A data frame with one variables:
#' \describe{
#' \item{\code{ENSG00000000003}}{Expression counts for ENSG00000000003},
#'}
#'
"BRCA_RNA_Test"
