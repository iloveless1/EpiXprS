#' Allows user to search provided models for CpG sites that are associated with
#' expression of genes in any of the eight cancers
#' @importFrom utils data
#' @param features CpG sites to be identified in the model of choice
#' @param Cancer The model to be searched for associated CpG sites
#' @return data.frame containing all instance of CpG sites in models
#' @export model_search
#' @examples
#' \dontrun{
#' model_search('cg21837192', Cancer = 'PRAD')
#' }

model_search <- function(features, Cancer=c('PRAD','BRCA','COAD','KIRP','KIRC',
                                            'HNSC','LUAD','UCEC')){

    if(Cancer == 'PRAD'){
        load(url("https://github.com/iloveless1/EPIExPRS-Data/raw/main/PRAD.rda"))
        rt = PRAD
    }else if (Cancer == 'BRCA'){
        load(url("https://github.com/iloveless1/EPIExPRS-Data/raw/main/BRCA.rda"))
        rt = BRCA
    }else if (Cancer == 'COAD'){
        load(url("https://github.com/iloveless1/EPIExPRS-Data/raw/main/COAD.rda"))
        rt = COAD
    }else if (Cancer == 'KIRP'){
        load(url("https://github.com/iloveless1/EPIExPRS-Data/raw/main/KIRP.rda"))
        rt = KIRP
    }else if (Cancer == 'KIRC'){
        load(url("https://github.com/iloveless1/EPIExPRS-Data/raw/main/KIRC.rda"))
        rt = KIRC
    }else if (Cancer == 'HNSC'){
        load(url("https://github.com/iloveless1/EPIExPRS-Data/raw/main/HNSC.rda"))
        rt = HNSC
    }else if (Cancer == 'LUAD'){
        load(url("https://github.com/iloveless1/EPIExPRS-Data/raw/main/LUAD.rda"))
        rt = LUAD
    }else if (Cancer == 'UCEC'){
        load(url("https://github.com/iloveless1/EPIExPRS-Data/raw/main/UCEC.rda"))
        rt = UCEC
    } else{stop("Cancer must be one of ('PRAD','BRCA','COAD','KIRP','KIRC','HNSC',
              'LUAD','UCEC')")}


features <- c(features,paste0(features,':Race'))
all <- NULL
genes <- NULL
for(i in seq_len(nrow(rt))){
    tmp <- as.data.frame(rt[[i,1]])

    if(length(which(features %in% tmp$V1) > 0)){
        genes <- cbind(tmp[which(tmp$V1 %in% features),],as.data.frame(rt[[i,2]]))
    }
    all <- rbind(all,genes)
}
tmp <- unique(all[c('V1','ensembl_gene')])
all <- all[rownames(all) %in% rownames(tmp),]
all <- all[,-c(3)]
colnames(all) <- c('CpG site', 'Effect on Expression', 'Ensembl Gene ID',
                   'Chromosome', 'TSS', 'TES', 'HGNC Symbol')
return(as.data.frame(all))
}
