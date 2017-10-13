#' helper function to compute variance weighted eucledian distance
#' \code{varWeightedEucledian} sums logs without loss of precision
#'
#' @param vexp a numeric vector of variance explained for a set of principle components (e.g summary(pc)$importance[2,])
#' @param a a vector of principle component loadings
#' @param b a vector of principle component loadings
#' @return a scalar

varWeightedEucledian<-function(vexp,a,b){
  if(length(a)!=length(b))
    stop("Loading vectors are not equal")
  if(length(a)!=vexp)
    stop("Loading and variance vectors are not equal")
  sqrt((a-b)^2 %*% vexp)
}

#' \code{pairwise_euc_matrix} compute a pairwise matrix of variance weighted distances across PC loadings
#'
#' @param x numeric matrix of loadings (e.g pc$x)
#' @param vexp a numeric vector of variance explained for a set of principle components (e.g summary(pc)$importance[2,])
#' @param pcidx (optional) parameter sepcifying a single principle component to use when computing distance (note no longer variance weighted)
#' @return a numeric matrix
#' @export

pairwise_euc_matrix<-function(x,vexp,pcidx){
  rnames<-rownames(x)
  M<-matrix(data=NA,nrow=length(rnames),ncol=length(rnames),dimnames=list(rnames,rnames))
  diag(M)<-0
  for(i in seq_along(rnames)){
    for(j in seq_along(rnames)){
      if(i!=j){
        if(missing(pcidx)){
          M[i,j] <- varWeightedEucledian(vexp,x[i,],x[j,])
        }else{
          M[i,j] <- varWeightedEucledian(1,x[i,pcidx],x[j,pcidx])
        }
      }
    }
  }
  return(M)
}
