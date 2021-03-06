% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathwayFisher.R
\name{pathwayFisher}
\alias{pathwayFisher}
\title{pathway analysis}
\usage{
pathwayFisher(significant, whole, fdr = 1.1, database = NULL,
  pathSizeMin = 15, pathSizeMax = 200)
}
\arguments{
\item{significant}{a vector of hit genes}

\item{whole}{a vector of background genes}

\item{fdr}{desired FDR cutoff. Only pathways passed the FDR cutoff will survive}

\item{database}{a list of pathway database, each list element contains a vector of genes}
}
\value{
pathway enrichment result table
}
\description{
pathway enrichment analysis using Fisher exaxt test
}
\details{
pathway enrichment analysis using Fisher exaxt test
}
\examples{
whole <- letters
significant <- letters[1:5]
database <- list(d1=letters[2:20],d2=letters[5:26])
pathwayFisher(significant, whole, database=database)

}
\author{
Caleb
}
