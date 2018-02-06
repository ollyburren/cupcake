#' compute counts for alternative (minor) allele within a population
#' \code{n_alt} compute expected count of alt alleles for a given sample
#'
#' @param N sample size
#' @param m minor allele freq
#' @return a vector of allele counts

n_alt<-function(N,m) N * m^2


#' compute counts for heterozygous allele within a population
#' \code{n_het} compute expected count of heterozygous\ alleles for a given sample
#'
#' @param N sample size
#' @param m minor allele freq
#' @return a vector of allele counts

n_het<-function(N,m) 2 * N * m * (1-m)

#' compute counts for different reference (major) allele within a population
#' \code{n_ref} compute expected count of ref alleles for a given sample
#'
#' @param N sample size
#' @param m minor allele freq
#' @return a vector of allele counts

n_ref<-function(N,m) N * (1-m)^2

#' compute expected value of outcome variable Y
#' \code{y_bar} compute expected count of ref alleles for a given sample
#'
#' @param b - beta estimate
#' @param N sample size
#' @param m minor allele freq
#' @return a vector of expected outcome values for Y


y_bar<-function(b,N,m) (b * (n_het(N,m) + (2 * n_alt(N,m))))/N

#' estimate allelic sum of squares
#' \code{ssx} estimate allelic sum of squares
#'
#' @param N sample size
#' @param m minor allele freq
#' @return a vector of ssx

ssx<-function(N,m) 2 * m * (N-1) * (1-m)

# standard error of epsilon
# N - Sample size
# m - Minor allele frequency
# seb - Standard error for beta


#' estimate standard error of epsilon
#' \code{sigma_e} estimate standard error of epsilon
#'
#' @param N sample size
#' @param seb - standard error of beta
#' @param m - minor allele freq
#' @return a vector of standard error of epsilon

sigma_e<-function(N,seb,m) seb * sqrt(ssx(N,m))

#' estimate P(Y>thresh|X=x)
#' \code{threshP} estimate the probability that Y exceeds a threshold value given an alleleic configuration
#'
#' @param thresh - threshold value to use.
#' @param X (allelic configuration) (0,1,2) assuming additive model
#' @param seb - standard error of beta
#' @param m - minor allele freq
#' @param N sample size
#' @return a vector of probabilities

threshP<-function(thresh,X,N,b,seb,m){
  SE.epsilon <- sigma_e(N,seb,m)
  pnorm(thresh,mean=X*b,sd=SE.epsilon,lower.tail = FALSE)
}


#' estimate P(Y>thresh|X=x) on odds scale
#' \code{threshOdds} estimate the probability that Y exceeds a threshold value given an alleleic configuration on the odds scale
#'
#' @param thresh - threshold value to use.
#' @param X (allelic configuration) (0,1,2) assuming additive model
#' @param seb - standard error of beta
#' @param m - minor allele freq
#' @param N sample size
#' @return a vector of odds


threshOdds<-function(thresh,X,N,b,seb,m){
  num<-threshP(thresh,X,N,b,seb,m)
  num/(1-num)
}


#' compute OR associated with a given linear coefficient and threshold value
#' \code{threshBeta}
#'
#' @param thresh - threshold value to use.
#' @param seb - standard error of beta
#' @param m - minor allele freq
#' @param N sample size
#' @return a vector of log(OR)


threshBeta<-function(thresh,N,b,seb,m){
  lyx1<-log(threshOdds(thresh,1,N,b,seb,m))
  lyx0<-log(threshOdds(thresh,0,N,b,seb,m))
  lyx1-lyx0
}


## compute the stderr of thresh beta
#' compute standard error associated with Beta
#' \code{threshSigmaBeta}
#'
#' @param thresh - threshold value to use.
#' @param m - minor allele freq
#' @param n1 - number of cases
#' @param n0 - number of controls
#' @param tb - number of controls
#' @return a vector standard errors


threshSigmaBeta<-function(thresh,n1,n0,m,tb){
  ## number above thresh
  N<-n0+n1
  ploidy<-2 # given current framework this is constant
  a<-(n0 * (1-m))/N
  b<-(n0 * m)/N
  c<-((n1*a)/(a+(b*exp(tb))))/N
  d<-((n1*b)/(a+(b*exp(tb))))/N
  var_maf<-sqrt(rowSums(do.call('cbind',lapply(list(a,b,c,d),function(fi) 1/fi))))
  var_ss<-sqrt(1/N)
  var_ploidy<-sqrt(1/ploidy)
  var_maf * var_ss * var_ploidy
}

## overall conversion routine

##  convert a linear beta to an OR using summary statistics
#' \code{convertBetaToOR}
#'
#' @param thresh - threshold value to use (if missing use mean of Y).
#' @param m - minor allele freq
#' @param b - beta coefficient from linear model fit
#' @param N - Number of samples
#' @param seb - standard error of beta coefficient
#' @return a vector of lists containing OR,betas and other intermediates
#' @export


convertBetaToOR<-function(thresh,N,b,seb,m){
  if (missing(thresh))
    thresh <- y_bar(b,N,m)
  tb<-threshBeta(thresh,N,b,seb,m)
  x_0<-threshP(thresh,0,N,b,seb,m) * n_ref(N,m)
  x_1<-threshP(thresh,1,N,b,seb,m) * n_het(N,m)
  x_2<-threshP(thresh,2,N,b,seb,m) * n_alt(N,m)
  n1<-round(x_0 + x_1 + x_2)
  n0<-N-n1
  tsb<-threshSigmaBeta(thresh,n1,n0,m,tb)
  Z<-tb/tsb
  P<-2*(pnorm(abs(Z),lower.tail = FALSE))
  list(OR=exp(tb),beta=tb,se.beta=tsb,P=P,Z=Z,thresh=thresh,n0=n0,n1=n1)
}
