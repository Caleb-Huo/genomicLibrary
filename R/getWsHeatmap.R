##' @export
getWsHeatmap <- function(astudy,aCs,ws=NULL, Rowv=NA,geneOrder=NULL,maxIndex=3,...){
   if(is.null(ws))
     ws=rep(1,nrow(astudy))
   astudy = astudy[ws!=0,]
   coll <- NULL
   resHeatmap <- NULL
   label <- NULL
   uniCs = unique(aCs)
   uniorderCs = sort(uniCs)

   for(i in uniorderCs){  
     if(is.null(resHeatmap)){
       resHeatmap = astudy[,i==aCs]
       label = rep(i,sum(i==aCs))
     } 
     else {
       resHeatmap = cbind(resHeatmap,astudy[,i==aCs])
       label = c(label,rep(i,sum(i==aCs)))
     }
   }
   for(alabel in unique(label))
     coll[which(label==alabel)]=palette()[alabel]    
   finalRes = t(scale(t(resHeatmap)))
   if(!is.null(geneOrder)){
     finalRes = finalRes[geneOrder,]
   }
   outIndex = abs(finalRes)>maxIndex 
   finalRes[outIndex] = maxIndex*sign(finalRes)[outIndex]
   B=16
   #return(gplots::heatmap.2(finalRes, col="greenred" ,trace="none",Rowv=Rowv,ColSideColors=coll,
   #                         Colv=NA,keysize=1.3,...)  )
   return(heatmap(finalRes,Rowv=Rowv,ColSideColors=coll,
 	  col= rgb(c(rep(0, B), (0:B)/B), c((B:0)/16, rep(0, B)), rep(0, 2*B+1))
   ,scale='none',Colv=NA,...)  )
 }
