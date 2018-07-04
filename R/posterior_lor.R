#' This function computes shape parameters for beta distribution matching a binomial proportion B(n,f)/n
#' \code{control_prior_shape} shape parameters for beta distribution matching a binomial proportion B(n,f)/n
#'
#' @param f a scalar - allele frequency in the controls
#' @param n a scalar - number if indivudals (not observations) used to estimate f
#' @return list of shape parameters for a beta distributon
#' \enumerate{
#' \item a0 - alpha shape parameter in controls
#' \item b0 - beta shape parameter in controls
#' }
#' @export
control_prior_shape <- function(f,n){
    v=f*(1-f)/(2*n)
    a0=-f*(f^2-f+v)/v
    b0=(f^2-f+v)*(f-1)/v
    list(a0=a0,b0=b0)
}

#' This function computes expected log odds ratio for a set of shape parameters
#' \code{e_lor} computes expected log odds ratio for a set of shape parameters
#'
#' @param b1 a scalar - beta distribution beta shape parameter for prior on allele frequency in cases
#' @param a0 a scalar - beta distribution alpha shape parameter for prior on allele frequency in controls
#' @param a0 a scalar - beta distribution beta shape parameter for prior on allele frequency in controls
#' @param a1 a scalar - beta distribution alpha shape parameter for prior on allele frequency in cases
#' @return scalar - expected log(OR)

e_lor <- function(b1,a0,b0,a1){
    abs(digamma(a0) - digamma(b0) - digamma(a1) + digamma(b1))
}

#' This function computes a probability for a given configuration of beta distribution
#' shape parameters estimate the P(OR > target.or).
#' \code{lor_constraint} estimate the P(OR > target.or) given shape parameters for case and control allele frequency
#'
#' @param a1 a scalar - beta distribution alpha shape parameter for prior on allele frequency in cases
#' @param a0 a scalar - beta distribution alpha shape parameter for prior on allele frequency in controls
#' @param b0 a scalar - beta distribution beta shape parameter for prior on allele frequency in controls
#' @param p0 a vector - a set of samples from beta(a0,b0)
#' @param target.or a scalar - an odds ratio threshold to compute P(sim.or > target.or)
#' @return scalar - P(sim.or>target.or)

lor_constraint <- function(a1,a0,b0,p0,target.or,nsim){
    b1 <- optimise(e_lor,c(1,b0),a1=a1,a0=a0,b0=b0)$minimum
    p1 <- rbeta(nsim,shape1=a1,shape2=b1)
    lor <- log(p0) - log(1-p0) + log(1-p1) - log(p1)
    mean(abs(lor) > log(target.or))
}

#' This function estimates shape parameters for a prior distribution of control allele frequencies
#' \code{est_a1b1} estimate shape parameters a1,b1
#' @param a0 a scalar - beta distribution alpha shape parameter for prior on allele frequency in controls
#' @param b0 a scalar - beta distribution beta shape parameter for prior on allele frequency in controls
#' @param target.or a scalar - an odds ratio threshold to compute P(sim.or > target.or)
#' @param target.prob a scalar - a probability that a sampled variant will exceed target.or
#' @return list - shape parameters for beta distribution satisfying target.or,target.prob and control prior distribution constraints.

est_a1b1 <- function(a0,b0,target.or,target.prob,nsim){
    # get an independendent sample of the prior of f_0
    p0 <- rbeta(nsim,shape1 = a0 ,shape2 = b0)
    ## estimate compatible a1 shape parameter for f1 given target OR and probability
    a1 <- 1
    while((pr <- lor_constraint(a1,a0,b0,p0,target.or,nsim)) > target.prob){
        if(pr==1){
            message(sprintf("Never satisfies %.2f %.2f",target.or,target.prob))
            return(list(a1=NULL,b1=NULL))
        }
        a1 <- a1 + 1
    }
    b1 <- optimise(e_lor,interval=c(1,b0),a1=a1,a0=a0,b0=b0)$minimum
    list(a1=a1,b1=b1)
}


#' This function samples from the posterior distribution of log(OR) based on the observation of a case genotype
#' \code{post_lor} samples from the posterior distribution of log(OR) given genotype.
#' @param gt a scalar - takes the values 0 - reference hom, 1 - het, 2 - alternative hom.
#' @param a1 a scalar - beta distribution alpha shape parameter for prior on allele frequency in cases
#' @param b1 a scalar - beta distribution beta shape parameter for prior on allele frequency in cases
#' @param p0 a vector - a set of samples from beta(a0,b0) i.e. prior on control allele frequency
#' @param nsim a scalar - number of log odds ratios to simulate from the posterior distribution
#' @return vector of simulated log(OR) from the posterior.

post_lor <- function(gt=c(0,1,2),a1,b1,p0,nsim){
    posta <- 2-gt + a1
    postb <- gt + b1
    p1 <- rbeta(nsim,shape1=posta,shape2=postb)
    list(lor=log(p0*(1-p1)/(p1*(1-p0))),p1=p1)
}

#' This function samples from the posterior distribution of a given control allele frequency, for different genotype configurations
#' \code{lor_f} sample from log(or) posterior distribution given an allele frequency.
#' @param f0 a scalar - allele frequency in cases to simulate
#' @param n a scalar - number of individuals sampled
#' @param nsim a scalar - number of simulations to perform to define the posterior and select parameters
#' @param target.or a scalar - an odds ratio threshold to compute P(sim.or > target.or)
#' @param target.prob a scalar - a probability that a sampled variant will exceed target.or
#' @return a list object
#' \enumerate{
#' \item 00 - mean lor for reference hom
#' \item 01 - mean lor for het
#' \item 11 - mean lor for alternative hom
#' \item quant00 - quantiles for posterior distribution of allele frequency in cases for ref hom
#' \item quant01 - quantiles for posterior distribution of allele frequency in cases for het
#' \item quant11 - quantiles for posterior distribution of allele frequency in cases for alt hom
#' }
#' @export

lor_f <- function(f0,n,nsim,target.or,target.prob,n.steps){
  p0.shape <- control_prior_shape(f0,n.sample)
  a0 <- p0.shape$a0
  b0 <- p0.shape$b0
  a1b1 <- opt_a1b1(a0,b0,target.or,target.prob,n.steps)
  a1 <- a1b1$a1
  b1 <- a1b1$b1
  c("00"=digamma(a0) - digamma(b0) - digamma(2+a1) + digamma(b1),
    "01"=digamma(a0) - digamma(b0) - digamma(1+a1) + digamma(1+b1),
    "11"=digamma(a0) - digamma(b0) - digamma(a1) + digamma(2+b1))
}

#' This function estimates shape parameters for a prior distribution of control allele frequencies
#' \code{opt_a1b1} estimate shape parameters a1,b1
#' @param a0 a scalar - beta distribution alpha shape parameter for prior on allele frequency in controls
#' @param b0 a scalar - beta distribution beta shape parameter for prior on allele frequency in controls
#' @param target.or a scalar - an odds ratio threshold to compute P(sim.or > target.or)
#' @param target.prob a scalar - a probability that a sampled variant will exceed target.or
#' @param n.steps an integer - length of search grid to employ - larger gives a more accurate integral estimate
#' @return list - shape parameters for beta distribution satisfying target.or,target.prob and control prior distribution constraints.

opt_a1b1 <- function(a0,b0,target.or,target.prob,n.steps){
    ## estimate compatible a1 shape parameter for f1 given target OR and probability
    p <- seq(0,1,length.out=n.steps)
    p <- p[-c(1,length(p))]
    a1 <- optimise(fopt,interval=c(1,a0),target.prob=target.prob,p=p,a0=a0,b0=b0,target.or=target.or)$minimum
    b1 <- optimise(e_lor,interval=c(1,b0),a1=a1,a0=a0,b0=b0)$minimum
    list(a1=a1,b1=b1)
}


#' This function computes the optimal shape parameter a1 for a given set of constraints
#' shape parameters estimate the P(OR > target.or).
#' \code{fopt} estimate the P(OR > target.or) given shape parameters for case and control allele frequency
#'
#' @param a1 a scalar - beta distribution alpha shape parameter for prior on allele frequency in cases
#' @param target.prob a scalar - a probability that a sampled variant will exceed target.or
#' @param target.or a scalar - an odds ratio threshold to compute P(sim.or > target.or)
#' @param p - a vector of values to optimise over
#' @param a0 a scalar - beta distribution alpha shape parameter for prior on allele frequency in controls
#' @param b0 a scalar - beta distribution beta shape parameter for prior on allele frequency in controls
#' @return scalar - P(sim.or>target.or)


fopt <- function(a1,target.prob,target.or,p,a0,b0) {
    b1 <- optimise(e_lor,interval=c(1,b0),a1=a1,a0=a0,b0=b0)$minimum
    denom <- dbeta(p,shape1=a0,shape2=b0)
    py <- sapply(p,function(x){
      l <- x/(target.or * (1-x) + x)
      u <- x/((1-x)/target.or + x)
      pbl <- pbeta(l,shape1=a1,shape2=b1,lower.tail=TRUE)
      pbu <- pbeta(u,shape1=a1,shape2=b1,lower.tail=FALSE)
      (pbl + pbu)
    })
    abs(sum(py*denom)/sum(denom) - target.prob)
}
