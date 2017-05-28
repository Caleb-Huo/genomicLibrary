##' @export
probeAnnotation <- function(probes, platformSYMBOL){
	x <- platformSYMBOL
	# Convert to a list
	xx <- base::as.list(x)
	
	matchProbeIndex = match(probes,names(xx))
	naive <- function(x){
		if(is.null(x)){
			return(NA)
		}
		return(x)
	}
	geneNames = sapply(xx[matchProbeIndex],naive)
	geneNames
}

## examples
## allProbes <- rownames(data_GSE4922)
## wholeGenes = probeAnnotation(allProbes)

