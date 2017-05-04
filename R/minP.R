##' @export
minP <- function(aP){
  return(pbeta(q=min(aP),shape1=1,shape2=length(aP)))
}
