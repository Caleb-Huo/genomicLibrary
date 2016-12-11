prepareFisherTable <- function(one_database,significant,whole)
{
  len.gene=length(whole)
  len.DEgene=length(significant)
  len.NonDEgene=len.gene-len.DEgene

  DE_IN=length(intersect(significant,one_database))
  DE_OUT=len.DEgene-DE_IN
  NonDE_IN=length(one_database)-DE_IN
  NonDE_OUT=len.NonDEgene-NonDE_IN
  FisherTable=matrix(c(DE_IN,DE_OUT,NonDE_IN,NonDE_OUT),ncol=2)
  return(FisherTable)
}#####this function is to prepare fot the fisher's exact test

