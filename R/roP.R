##' @export
roP <- function(aP,rth){
  return(pbeta(q=sort(aP)[rth],shape1=rth,shape2=length(aP)-rth+1))
}
