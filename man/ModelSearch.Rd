% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModelSearch.R
\name{ModelSearch}
\alias{ModelSearch}
\title{Allows user to search provided models for CpG sites that are associated with
expression of genes in any of the cancer models.}
\usage{
ModelSearch(features, Cancer = object)
}
\arguments{
\item{features}{CpG sites to be identified in the model of choice}

\item{Cancer}{ExperimentHub object cancer model}
}
\value{
data.frame containing all instance of CpG sites in models.
\code{CpG.site} The CpG site included in the model and if it has race specific effects
\code{Weight} The weight of the CpG site on the expression of the gene
\code{Ensembl.ID} The Ensembl ID
\code{Chromosome} The chromosome of the gene
\code{TSS} The transcription start site of the gene
\code{TES} The transcription end site of the gene
\code{Gene} The gene name
\code{Overall.D2} The prediction accuracy for the race interaction model
\code{EA.D2} The European American specific prediction accuracy
\code{AA.D2} The African American specific prediction accuract
}
\description{
Allows user to search provided models for CpG sites that are associated with
expression of genes in any of the cancer models.
}
\examples{
library(ExperimentHub)
eh = ExperimentHub()
x <- query(eh, 'EpiXprS')
Cancer <- x[["EH6022"]]
ModelSearch('cg21837192', Cancer = Cancer)

}
