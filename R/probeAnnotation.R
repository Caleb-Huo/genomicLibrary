##' @export
probeAnnotation <- function(probes, platformSYMBOL){
	x <- platformSYMBOL
	# Get the probe identifiers that are mapped to a gene symbol
	mapped_probes <- mappedkeys(x)
	# Convert to a list
	xx <- as.list(x[mapped_probes])
	
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

