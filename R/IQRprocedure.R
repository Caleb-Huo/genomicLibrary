##' @export
IQRprocedure <- function(wholeGenes, datamatrix){
	stopifnot(length(wholeGenes) == nrow(datamatrix))

	geneSplit = split(1:length(wholeGenes),wholeGenes)

	IQRMatrix <- function(Matr)
	{
	  if(!is.matrix(Matr)) return (Matr)
	  iqrSconre=apply(Matr,1,IQR)
	  return (Matr[which.max(iqrSconre),])
	}

	data_genes <- t(sapply(geneSplit,function(x) IQRMatrix(datamatrix[x,])))
	data_genes
}

## examples
## GSE4922_genes <- IQRprocedure(wholeGenes, data_GSE4922)




