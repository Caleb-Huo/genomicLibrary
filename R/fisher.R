##' @export
fisher <- function(aP){
  fisherDf = 2*length(aP)
  S = -2*sum(log(aP))
  return(pchisq(q=S,df=fisherDf,lower.tail=FALSE))
}
