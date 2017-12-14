library(data.table)
library(GenomicRanges)
library(magrittr)



#' helper function to sum logs without loss of precision
#' \code{logsum} sums logs without loss of precision
#'
#' @param x a vector of logs to sum
#' @return a scalar

logsum <- function(x) {
    my.max <- max(x) ##take out the maximum value in log form)
    my.res <- my.max + log(sum(exp(x - my.max )))
    return(my.res)
}

#' compute posterior probabilities using Wakefield's approximate Bayes Factors
#' \code{wakefield_pp} computes posterior probabilities for a given SNP to be causal for a given SNP under the assumption of a single causal variant.
#'
#' @param p a vector of univariate pvalues from a GWAS
#' @param f a vector of minor allele frequencies taken from some reference population.
#' @param N a scalar or vector for total sample size of GWAS
#' @param s a scalar representing the proportion of cases (n.cases/N)
#' @param pi a scalar representing the prior probability (DEFAULT \eqn{1 \times 10^{-4}})
#' @param sd.prior a scalar representing our prior expectation of \eqn{\beta} (DEFAULT 0.2).
#' The method assumes a normal prior on the population log relative risk centred at 0 and the DEFAULT
#' value sets the variance of this distribution to 0.04, equivalent to a 95\%  belief that the true relative risk
#' is in the range of 0.66-1.5 at any causal variant.
#' @return a vector of posterior probabilities.
#' @export

wakefield_pp <- function(p,f, N, s,pi=1e-4,sd.prior=0.2) {
    if(length(p) != length(f))
      stop("p and f must be vectors of the same size")
    # compute V
    V <- 1 / (2 * N * f * (1 - f) * s * (1 - s))
    # convert p vals to z
    z <- qnorm(0.5 * p, lower.tail = FALSE)
    ## Shrinkage factor: ratio of the prior variance to the total variance
    r <- sd.prior^2 / (sd.prior^2 + V)
    ## Approximate BF
    lABF = 0.5 * (log(1-r) + (r * z^2))
    ## tABF - to add one we create another element at the end of zero for which pi_i is 1
    tABF <- c(lABF,0)
    vpi_i<-c(rep(pi_i,length(lABF)),1)
    sBF <- logsum(tABF + log(vpi_i))
    exp(lABF+log(pi_i)-sBF)
}

#' compute reciprocal posterior probabilities using Wakefield's approximate Bayes Factors
#' \code{wakefield_null_pp} computes posterior probabilities for a given SNP to be NOT be causal for a given SNP under the assumption of a single causal variant.
#'
#' @param p a vector of univariate pvalues from a GWAS
#' @param f a vector of minor allele frequencies taken from some reference population.
#' @param N a scalar or vector for total sample size of GWAS
#' @param s a scalar representing the proportion of cases (n.cases/N)
#' @param pi_i a scalar representing the prior probability (DEFAULT \eqn{1 \times 10^{-4}})
#' @param sd.prior a scalar representing our prior expectation of \eqn{\beta} (DEFAULT 0.2).
#' The method assumes a normal prior on the population log relative risk centred at 0 and the DEFAULT
#' value sets the variance of this distribution to 0.04, equivalent to a 95\%  belief that the true relative risk
#' is in the range of 0.66-1.5 at any causal variant.
#' @return a vector of posterior probabilities.
#' @export

wakefield_null_pp <- function(p,f, N, s,pi_i=1e-4,sd.prior=0.2) {
    if(length(p) != length(f))
      stop("p and f must be vectors of the same size")
    # compute V
    V <- 1 / (2 * N * f * (1 - f) * s * (1 - s))
    # convert p vals to z
    z <- qnorm(0.5 * p, lower.tail = FALSE)
    ## Shrinkage factor: ratio of the prior variance to the total variance
    r <- sd.prior^2 / (sd.prior^2 + V)
    ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
    lABF = 0.5 * (log(1-r) + (r * z^2))
    po = exp(lABF + log(pi_i) - log(1-pi_i))
    pp = po/(1+po)
}

#' This function computes the posterior prob that a SNP is causal in a set of traits
#' \code{basis_pp} computes the posterior probability that a SNP is causal across a set of traits
#'
#' @param bf a vector of approximate Bayes Factors using Wakefield's method.
#' @param a scalar or vector of posterior probabilites
#' @return pi a scalar empirical prior

basis_pp<-function(BF,emp_pi){
  lABF<-log(BF)
  tABF <- c(lABF,0)
  vpi_i<-c(rep(emp_pi,length(lABF)),1)
  sBF <- logsum(tABF + log(vpi_i))
  exp(lABF+log(emp_pi)-sBF)
}

#' This function computes an approximate Bayes Factor for whether a given SNP is causal
#' in a set of basis traits assuming a single causal variant.
#' \code{bayesian_shrinkage} computes an approximate Bayes Factor for whether a given SNP is causal
#' in a set of basis traits assuming a single causal variant. Using an emprically derived
#' prior computes a posterior probability.
#'
#' @param DT basis data.table object
#' @param tquant a scalar representing quantile on which to truncate infinite Bayes Factors (DEFAULT 0.9999)
#' @return data.table object

bayesian_shrinkage<-function(DT,tquant=0.9999){
  tmp<-DT[,list(pid=pid,lp0=log(1-wakefield_null_pp(p.val,maf,n,n1/n))),by=c('trait','ld.block')]
  # if pvalue is 1 (as OR is 1) we get numerical errors / NA. Removing is fine as the they should be 0 therefore 1
  tmp<-tmp[,list(q_i=1-exp(sum(lp0,na.rm=TRUE))),by=c('ld.block','pid')]
  ## compute an empirical prior
  emp<-mean(tmp$q_i)
  ## prior odds
  po<-emp/(1-emp)
  ## note that emp is for h1 that beta != 0 therefore we need to take reciprocal as equation assumes pi_0 - see notes
  po<-1/po
  tmp[,uABF:=po*(q_i/(1-q_i))]
  ## set an upper limit to BF (here its upper 0.0001 percentile)
  BIG<-quantile(tmp[is.finite(tmp$uABF),]$uABF,prob=0.9999)
  tmp[is.infinite(uABF) | uABF > BIG, uABF:=BIG]
  tmp[,bshrink:=basis_pp(uABF,emp),by=ld.block][,.(bshrink),by=pid]
}

#' This function computes an estimate of allele count of unexposed controls
#' \code{ca} estimate of allele count of unexposed controls
#'
#' @param n0 a vector or scalar of number of control samples
#' @param f a vector of reference allele frequencies
#' @return a numeric vector

ca<-function(n0,f){
    n0*(1-f)
}

#' This function computes an estimate of allele count of  exposed controls
#' \code{cb} estimate of allele count of exposed controls
#'
#' @param n0 a vector or scalar of number of control samples
#' @param f a vector of reference allele frequencies
#' @return a numeric vector

cb<-function(n0,f){
    n0*f
}

#' This function computes an estimate of allele count of  unexposed cases
#' \code{cc} estimate of allele count of unexposed cases
#'
#' @param n1 a vector or scalar of number of control samples
#' @param a a vector of estimates for allele count of  unexposed controls
#' @param b a vector of estimates for allele count of  exposed controls
#' @return a numeric vector
#' see also \code{\link{ca}} and \code{\link{cb}}

cc<-function(n1,a,b,theta){
    (n1*a)/(a+(b*theta))
}

#' This function computes an estimate of allele count of  exposed cases
#' \code{cc} estimate of allele count of exposed cases
#'
#' @param n1 a vector or scalar of number of control samples
#' @param a a vector of estimates for allele count of  unexposed controls
#' @param b a vector of estimates for allele count of  exposed controls
#' @return a numeric vector
#' see also \code{\link{ca}} and \code{\link{cb}}

cd<-function(n1,a,b,theta){
    (n1*b)/(a+(b*theta))
}


#' This function standard error due to minor allele frequency empirically
#' \code{maf_se_empirical} emprically estimate the standard error that is due to minor allele frequency
#' \eqn{\sqrt{\frac{1}{a} +  \frac{1}{b} + \frac{1}{b} + \frac{1}{b} }}
#'
#' @param n0 a vector or scalar of number of control samples
#' @param n1 a vector or scalar of number of case samples
#' @param f a vector of reference allele frequencies
#' @param a vector of Odds Ratios
#' @return a numeric vector
#' see also \code{\link{ca}}, \code{\link{cb}}, \code{\link{cc}} and \code{\link{cd}}

maf_se_empirical<-function(n0,n1,f,theta){
    n<-n0+n1
    a<-ca(n0,f)/n
    b<-cb(n0,f)/n
    c<-cc(n1,a,b,theta)/n
    d<-cd(n1,a,b,theta)/n
    recip.sm<-do.call('cbind',lapply(list(a,b,c,d),function(fi) 1/fi))
    sqrt(rowSums(recip.sm))
    #return(log(theta)/res)
}

#' Compute minor allele frequency shrinkage
#' \code{maf_se_estimate} computes a shrinkage metric for a given list of minor allele frequencies'
#'
#' @param f a vector of minor allele frequencies taken from some reference population.
#' @return a vector of shrinkage metrics

maf_se_estimate <- function(f){
  #1/sqrt(f * (1-f))
  sqrt(1/f + 1/(1-f))
}


#' convert p value  to a signed Z score
#' \code{p2z} p value to a signed Z score
#'
#' @param p a vector of p values
#' @param lor a vector of log odds ratios
#' @return a vector of signed Z scores

p2z <- function(p,lor){
  z <- qnorm(0.5 * p.val, lower.tail = FALSE)
  if(missing(lor))
    return(z)
  return(z * sign(lor))
}

#' convert z to p value
#' \code{p2z} z to p value
#'
#' @param z a vector of Z scores
#' @return a vector of p values

z2p <- function(z){
  2* pnorm(abs(z), lower.tail = FALSE)
}

# this function reads in maf data

add_ref_maf <- function(snp_support_file,DT){
  if(!file.exists(snp_support_file))
    stop(sprintf("Cannot find file %s",snp_support_file))
  ss<-fread(snp_support_file)
  ## use data table to merge the two files
  ss[,pid:=paste(chr,position,sep=':')]
  ss<-ss[,.(pid,maf)]
  setkey(ss,pid)
  tmp<-DT[ss]
  if(nrow(tmp)!=nrow(DT))
    stop("Something went wrong perhaps there are duplicates (by position) in your snp support file or in GWAS input")
  return(tmp)
}

add_ld_block <- function(ld_support_file,DT){
  if(!file.exists(ld_support_file))
    stop(sprintf("Cannot find file %s",ld_support_file))
  ld<-fread(ld_support_file)
  ld.gr<-with(ld,GRanges(seqnames=Rle(chr),ranges=IRanges(start=as.numeric(start),end=as.numeric(end))))
  snps.gr<-with(DT,GRanges(seqnames=Rle(chr),ranges=IRanges(start=position,width=1)))
  ol<-as.matrix(findOverlaps(snps.gr,ld.gr))
  DT[ol[,1],ld.block := ol[,2]]
}

# this function gets GWAS data using a manifest file. If a trait list is supplied
# then gets just those traits, if trait list is missing assumes that you want just
# basis
#' \code{get_gwas_data} integrate GWAS summary data with support files
#' @param manifest_file character vector file path to GWAS manifest file
#' @param snp_manifest_file character vector file path to snp manifest
#' @param ld_support_file character vector file path to ld regions
#' @param data_dir character vector file path to location of GWAS summary stats
#' @param trait_list character vector of specific traits in manifest file to include
#' @return data.table object
#' @export

get_gwas_data <- function(manifest_file,snp_manifest_file,ld_support_file,data_dir,trait_list){
  if(missing(trait_list)){
    man<-fread(manifest_file)[basis_trait==1,]
  }else{
    man<-fread(manifest_file)[trait %in% trait_list,]
  }
  if(nrow(man)==0)
    stop(sprintf("Cannot find any traits in manifest %s for %s",manifest_file,paste(trait_list,collapse=',')))
  man[,file:=file.path(data_dir,file)]
  ret<-rbindlist(lapply(1:nrow(man),function(i){
    message(sprintf("Processing %s",man[i,]$trait))
    tDT<-fread(man[i,]$file)
    tDT[,pid:=paste(chr,position,sep=':')]
    tDT[,c('trait','n','n1') := man[i,.(trait,cases+controls,cases)]]
  }))
  setkey(ret,pid)
  ## next add minor allele frequencies
  message("Adding reference minor allele freq.")
  ret<-add_ref_maf(snp_manifest_file,ret)
  message("Assigning LD Blocks")
  ret<-add_ld_block(ld_support_file,ret)
  ret
}


#' This function computes various shrinkage metrics
#' \code{compute_shrinkage_metrics} computes various shrinkage metrics
#'
#' @param data.table object for basis traits as returned by \code{\link{get_gwas_data}}
#' @return a data.table object.
#' \enumerate{
#' \item pid - unique id using chr and position (useful for merging back)
#' \item bshrink - Bayesian shrinkage based on association across all basis traits
#' \item emp_maf_se - empirically derived standard error for MAF
#' \item est_maf_se - estimated standard error for MAF
#' \item emp_shrinkage - overall shrinkage using emp_maf_se
#' \item est_shrinkage - overall shrinkage using est_maf_se
#' }
#' see also \code{\link{maf_se_empirical}}, \code{\link{maf_se_estimate}} and \code{\link{bayesian_shrinkage}}.
#' @export

compute_shrinkage_metrics<-function(DT){
  message("Computing maf_se_empirical")
  emp_maf_se.DT<-DT[,list(pid=pid,emp_maf_se=maf_se_empirical(n-n1,n1,maf,or))][,list(emp_maf_se=mean(emp_maf_se)),by=pid]
  setkey(emp_maf_se.DT,pid)
  ## second way to do it is to compute based on function fitting.
  message("Computing maf_se_estimated")
  est_maf_se.DT<-unique(DT[,list(est_maf_se=maf_se_estimate(maf)),by=pid])
  setkey(est_maf_se.DT,pid)
  maf_se.DT<-emp_maf_se.DT[est_maf_se.DT]
  ## next compute basis shrinkage vector
  message("Computing Bayesian shrinkage")
  bs.DT<-bayesian_shrinkage(DT)
  setkey(bs.DT,pid)
  shrinkage.DT<-bs.DT[maf_se.DT]
  shrinkage.DT[,c('emp_shrinkage','est_shrinkage'):=list(bshrink/emp_maf_se,bshrink/est_maf_se),by=pid]
  setkey(shrinkage.DT,pid)
  return(shrinkage.DT)
}

#' This function creates a trait snp matrix
#' \code{create_ts_matrix} creates a trait snp matrix that is suitable for basis creation and projection
#'
#' @param bDT data.table object for basis traits as returned by \code{\link{get_gwas_data}}
#' @param sDT data.table object of matching shrinkage estimates returned by \code{\link{compute_shrinkage_metrics}}
#' @param method scalar vector (either emp or est), emp uses empirically generated MAF SE, est uses and estimate.
#' @return a matrix.
#' @export

create_ds_matrix <- function(bDT,sDT,method=c('emp','est')){
  if(missing(method)){
    method='emp'
  }
  message(sprintf("Using %s",method))
  vmethod = sprintf("%s_shrinkage",method)
  stmp<-sDT[,c('pid',vmethod),with=FALSE]
  tmp<-bDT[stmp]
  tmp$metric <- tmp[[vmethod]] * log(tmp$or)
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  tmp.mat <- as.matrix(B[,-1]) %>% t()
  colnames(tmp.mat) <- snames
  return(tmp.mat)
}
