% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Helpers.R
\name{D2}
\alias{D2}
\title{Calculate Accuracy}
\usage{
D2(y_i, u_i)
}
\arguments{
\item{y_i}{Observed counts}

\item{u_i}{Predicted counts}
}
\value{
the McFadden's Pseudo R-squared
}
\description{
Calculate Accuracy
}
\examples{
y_i <- rpois(50,12)
u_i <- rpois(50,11.5)
D2(y_i,u_i)
}
