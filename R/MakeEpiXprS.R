#' Creates a MultiAssay object containing RNA counts, DNA methylation, and 
#' clinical data. 
#' 
#' @import MultiAssayExperiment
#' @param RNASeqGene A matrix of RNA counts
#' @param Methyl A matrix of methylation values
#' @param Clinical A data.frame of clinical covariates
#' @param samples A vector of sample names ordered the same as RNASeqGene and 
#' Methyl450k
#' @return An object of class MultiAssayExperiment containing paired RNA counts 
#' and DNA methylation
#' 
#' @examples 
#' data('Testdat', package = 'EpiXprS')
#' Methyl <- assays(Testdat)$Methyl
#' Counts <- assays(Testdat)$RNASeqGene
#' Clinical <- colData(Testdat)
#' samples <- sampleMap(Testdat)$primary[1:100]
#' Data <- MakeEpiXprS(RNASeqGene = Counts, Methyl = Methyl, 
#' Clinical = Clinical, samples = samples)

MakeEpiXprS <- function(RNASeqGene = Counts, Methyl = Methyl, Clinical = Clinical, samples = samples){
    if(ncol(RNASeqGene) != ncol(Methyl450k)){
        stop('Must include same samples')
    }
    if(ncol(RNASeq) != nrow(Clinical)){
        stop('Must include same samples')
    }
    
    
    
    names(RNASeqGene) <- samples
    names(Methyl) <- samples
    rownames(Clinical) <- samples
    sampmap <- data.frame(assay = c(rep('RNASeqGene',ncol(RNASeqGene)),rep('Methyl',ncol(Methyl))), 
                          primary = rep(samples,2),
                          colname = c(colnames(RNASeqGene),colnames(Methyl)))
    
    
    assayList = list(RNASeqGene = RNASeqGene, Methyl = Methyl)
    
    Epi.Object <- MultiAssayExperiment(experiments=ExperimentList(assayList),
                                    colData=Clinical,
                                    sampleMap=sampmap)
    return(Epi.Object)
}
