##' @export
maxP <- function(aP){
  return(pbeta(q=max(aP),shape1=length(aP),shape2=1))
}
