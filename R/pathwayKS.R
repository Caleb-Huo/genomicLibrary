##' pathway enrichment analysis using KS test
##'
##' pathway enrichment analysis using KS exaxt test
##' @title pathway analysis
##' @param wholePvalues a vector of p values for all genes
##' @param wholeName a vector of all genes symbols
##' @param fdr desired FDR cutoff. Only pathways passed the FDR cutoff will survive
##' @param database a list of pathway database, each list element contains a vector of genes
##' @return pathway enrichment result table
##' @author Caleb
##' @export
##' @examples
##' wholePvalues <- rnorm(26)
##' wholeName <- letters[1:26]
##' database <- list(d1=letters[2:20],d2=letters[5:26])
##' pathwayKS(wholePvalues, wholeName, database=database, fdr=1)
##'
pathwayKS <- function(wholePvalues,wholeName,fdr=1.1,database=NULL,pathSizeMin=15,pathSizeMax=200){
    ####wholePvalues: p values for all genes
    ####wholeName: all genes symbols
  if(is.null(fdr)) fdr=0.05
  if(is.null(database)){
    return ("Please specify database")
  }
	
	if(F){
	  whole <- toupper(wholeName)
	  database <- lapply(database,function(x) toupper(x))		
	}

  ######################Update the gene sets by dropping genes that do not appear in whole
  gene.overlap2=lapply(database,function(x) intersect(x,whole))
  gene.overlap=levels(factor(unlist(gene.overlap2)))

  #####################Filter out gene sets that contain less than 5 genes or more than 200 genes
  genesets <- lapply(database,function(x) intersect(x,gene.overlap))
  index=which(sapply(genesets,function(x) length(x)<=pathSizeMax & length(x)>=pathSizeMin))
  genesets_pro <- genesets[index]
  names(genesets_pro) <-names(genesets)[index]
  

  path_pval <- sapply(genesets_pro,function(x){
	  sigMatchIndex <- match(x, whole)
	  ks.test(wholePvalues[sigMatchIndex], wholePvalues[-sigMatchIndex],alternative="greater")$p.value
  } )
  
  path_qval <- p.adjust(path_pval, method = "BH")
  sortIndex = order(path_pval)
  gene_set <- names(genesets_pro)

  Pathway = gene_set
  Path_pval = path_pval
  Path_qval = path_qval
  NumOfGenesInPath = sapply(genesets_pro,length)
  res = data.frame(Pathway,Path_pval,Path_qval,NumOfGenesInPath)
  
  resOrder <- res[sortIndex,]
  
  sigIndex <- resOrder$Path_qval<=fdr
  
  if(sum(sigIndex)>0){
	  return(resOrder[resOrder$Path_qval<=fdr, ])  	
  } else {
	  return(0)
  }
}

