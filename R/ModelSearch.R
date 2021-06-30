#' Allows user to search provided models for CpG sites that are associated with
#' expression of genes in any of the cancer models. 
#' @importFrom utils data
#' @import stringr
#' @param features CpG sites to be identified in the model of choice
#' @param Cancer ExperimentHub object cancer model
#' @param type Feature type to search for0
#' @return data.frame containing all instance of CpG sites in models.
#' \code{CpG.site} The CpG site included in the model and if it has race specific effects
#' \code{Weight} The weight of the CpG site on the expression of the gene
#' \code{Ensembl.ID} The Ensembl ID
#' \code{Chromosome} The chromosome of the gene
#' \code{TSS} The transcription start site of the gene
#' \code{TES} The transcription end site of the gene
#' \code{Gene} The gene name
#' \code{Overall.D2} The prediction accuracy for the race interaction model
#' \code{EA.D2} The European American specific prediction accuracy
#' \code{AA.D2} The African American specific prediction accuract
#' @export ModelSearch
#' @examples
#' library(ExperimentHub)
#' eh = ExperimentHub()
#' x <- query(eh, 'EpiXprS')
#' Cancer <- x[["EH6022"]]
#' ModelSearch('cg21837192', Cancer = Cancer)
#'

ModelSearch <- function(features, Cancer=object, type = c('CpG','Gene')){

   
all_list <- list()
for(i in seq_len(length(Cancer))){
    if(length(which(substr(rownames(Cancer[[i]][[1]]),1,str_length(features[1])) %in% features)) >0){
        k  <- length(which(substr(rownames(Cancer[[i]][[1]]),1,str_length(features[1])) %in% features))
        out <- data.frame('CpG Site' = rownames(Cancer[[i]][[1]])[substr(rownames(Cancer[[i]][[1]]),1, 
                                        str_length(features[1])) %in% features],
                          'Weight' = Cancer[[i]][[1]][substr(rownames(Cancer[[i]][[1]]),1, 
                                                str_length(features[1])) %in% features,],
                          'Ensembl ID' = rep(unlist(Cancer[[i]][[2]][2]),k),
                          'Chromosome' = rep(unlist(Cancer[[i]][[2]][3]),k),
                          'TSS' = rep(unlist(Cancer[[i]][[2]][4]),k),
                          'TES' = rep(unlist(Cancer[[i]][[2]][5]),k),
                          'Gene' = rep(unlist(Cancer[[i]][[2]][6]),k),
                          'Overall D2' = rep(unlist(Cancer[[i]][[2]][7]),k),
                          'EA D2' = rep(unlist(Cancer[[i]][[2]][8]),k),
                          'AA D2' = rep(unlist(Cancer[[i]][[2]][9]),k))
                          
        all_list[[i]] <- out
        
    } 
}
all <- do.call(rbind,all_list)

return(data.frame(all, row.names = c(1:nrow(all))))
}
