##' pathway enrichment analysis using Fisher exaxt test
##'
##' pathway enrichment analysis using Fisher exaxt test
##' @title pathway analysis
##' @param significant a vector of hit genes
##' @param whole a vector of background genes
##' @param fdr desired FDR cutoff. Only pathways passed the FDR cutoff will survive
##' @param database a list of pathway database, each list element contains a vector of genes
##' @return pathway enrichment result table
##' @author Caleb
##' @export
##' @examples
##' whole <- letters
##' significant <- letters[1:5]
##' whole0 <- runif(length(whole),0.05,1)
##' whole0[1:length(significant)] <- runif(length(significant),0,0.05)
##' names(whole0) <- whole
##' database <- list(d1=letters[2:20],d2=letters[5:26])
##' pathwayFisherAUC(significant, whole0, database=database)
##'
pathwayFisherAUC <- function(significant,whole0,fdr=1.1,database=NULL,pathSizeMin=15,pathSizeMax=200){
  ####significant: just genenames
  ####whole: p values with genenames
  if(is.null(fdr)) fdr=0.05
  if(is.null(database)){
    return ("Please specify database=Mm.gmtl.c2 or database=Mm.gmtl.c5")
  }
  
  whole <- names(whole0)
  significant <- toupper(significant)
  whole <- toupper(whole)
  database <- lapply(database,function(x) toupper(x))

  ######################Update the gene sets by dropping genes that do not appear in whole
  gene.overlap2=lapply(database,function(x) intersect(x,whole))
  gene.overlap=levels(factor(unlist(gene.overlap2)))
  #####################Filter out gene sets that contain less than 5 genes or more than 200 genes
  genesets <- lapply(database,function(x) intersect(x,gene.overlap))
  index=which(sapply(genesets,function(x) length(x)<=pathSizeMax & length(x)>=pathSizeMin))
  genesets_pro <- genesets[index]
  names(genesets_pro) <-names(genesets)[index]

  auc <- sapply(genesets_pro, function(x) roc(whole %in% x, 1 - whole0, direction='<')$auc)

  path_pval <- sapply(genesets_pro,function(x) fisher.test(prepareFisherTable(x,significant,whole),alternative="greater")$p.value)
  path_qval <- p.adjust(path_pval, method = "BH")
  sortIndex = order(path_pval)
  gene_set <- names(genesets_pro)

  Pathway = gene_set
  Path_pval = path_pval
  Path_qval = path_qval
  TotalNumOfDEGenes = length(significant)
  NumOfGenesInPath = sapply(genesets_pro,length)
  NumOfDEGenesInPath = sapply(genesets_pro,function(x) length(intersect(x,significant)))
  DEGenesInPath = sapply(genesets_pro,function(x) paste(intersect(x,significant),collapse='/'))
  res = data.frame(Pathway,Path_pval,Path_qval,auc,TotalNumOfDEGenes,NumOfGenesInPath,NumOfDEGenesInPath,DEGenesInPath)
  
  resOrder <- res[sortIndex,]
  
  sigIndex <- resOrder$Path_qval<=fdr
  
  if(sum(sigIndex)>0){
	  return(resOrder[resOrder$Path_qval<fdr, ])  	
  } else {
	  return(0)
  }
}

