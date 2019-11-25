#' @import data.table


library(data.table)

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
#' @param pi_i a scalar representing the prior probability (DEFAULT \eqn{1 \times 10^{-4}})
#' @param sd.prior a scalar representing our prior expectation of \eqn{\beta} (DEFAULT 0.2).
#' The method assumes a normal prior on the population log relative risk centred at 0 and the DEFAULT
#' value sets the variance of this distribution to 0.04, equivalent to a 95\%  belief that the true relative risk
#' is in the range of 0.66-1.5 at any causal variant.
#' @return a vector of posterior probabilities.
#' @export

wakefield_pp <- function(p,f, N, s,pi_i=1e-4,sd.prior=0.2) {
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



#' This function computes an alternative to the Bayesian shrinkage method which can be too agressive.
#' \code{ws_shrinkage} computes a shrinkage based on a weighted sum (ws) of posteriors for each disease
#' this is then normalised by the total posterior for a given LD block
#'
#' @param DT basis data.table object
#' @param pi_i - The prior probability that a variant is causal
#' @return data.table object

ws_shrinkage <- function(DT,pi_i=1e-4){
  tmp <- DT[,list(pid=pid,ppi=wakefield_pp(p.value,maf,n,n1/n,pi_i)),by=c('trait','ld.block')]
  tmp[,wj:=sum(ppi),by=c('trait','ld.block')]
  tmp[,list(ws_ppi=sum(ppi * wj)/sum(wj)),by=c('pid','ld.block')]
}



#' This function computes standard error under the null
#' \code{se_null} analytically compute standard error of \eqn{\beta} under \eqn{\mathbb{E}(\beta) = 0}
#' \eqn{\sqrt{\frac{1}{2fN_0} +  \frac{1}{2N_{0}(1-f)} + \frac{1}{2fN_1} +  \frac{1}{2N_{1}(1-f)} }}
#'
#' @param N a vector or scalar of total number od samples
#' @param n1 a vector or scalar of number of case samples
#' @param f a vector of reference allele frequencies
#' @return a numeric vector
#' @export

se_null<-function(N,n1,f){
  n0<-N-n1
  a<-1/(2*f*n0)
  b<-1/(2*(1-f)*n0)
  c<-1/(2*f*n1)
  d<-1/(2*(1-f)*n1)
  sqrt(rowSums(cbind(a,b,c,d)))
}


#' Compute minor allele frequency shrinkage
#' \code{maf_se_estimate} computes a shrinkage metric for a given list of minor allele frequencies'
#'
#' @param f a vector of minor allele frequencies taken from some reference population.
#' @return a vector of shrinkage metrics

maf_se_estimate <- function(f){
  #1/sqrt(f * (1-f))
  sqrt(1/f + 1/(1-f)) * 2
}

#' Compute minor allele frequency shrinkage using sample size
#' \code{maf_se_estimate_sample_size} computes component of standard error of beta due to minor allele frequency
#'
#' @param N a vector or scalar of total number od samples
#' @param p a vector or scalar of p values
#' @param theta a vector or scalar of odds ratios
#' @param f a vector of reference allele frequencies
#' @return a numeric vector

maf_se_estimate_sample_size <- function(N,p,theta,f){
  Z <- qnorm(p/2,lower.tail=FALSE)
  se.beta <- log(theta)/Z
  se_maf_ss <- sqrt(2 * N) * se.beta
  ## can get numeric errors if theta = 1 or such like in this case compute using maf_se estimate under null
  idx <- which(is.infinite(se_maf_ss) | is.nan(se_maf_ss) | se_maf_ss>100)
  se_maf_ss[idx] <- maf_se_estimate(f[idx])
  se_maf_ss
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

# this function gets adds reference data from a snp support data.table to GWAS summ stats
#' \code{add_ref_annotations} integrate GWAS summary data with support file
#' @param ss data.table object for snp manifest
#' @param DT data.table containing GWAS summary stats
#' @return data.table object

add_ref_annotations <- function(ss,DT){
  #if(!file.exists(snp_support_file))
  #  stop(sprintf("Cannot find file %s",snp_support_file))
  #ss<-fread(snp_support_file)
  ## use data table to merge the two files
  #ss[,pid:=paste(chr,position,sep=':')]
  ss[,maf:=ifelse(ref_a1.af>0.5,1-ref_a1.af,ref_a1.af)]
  ss<-ss[,.(pid,maf,ld.block)]
  setkey(ss,pid)
  # we filter here as this allows us to use this to
  # knock traits where we have surplus SNPs into the correct format
  tmp<-DT[ss][!is.na(or),]
  if(nrow(tmp)!=nrow(DT))
    stop("Something went wrong perhaps there are duplicates (by position) in your snp support file or in GWAS input")
  return(tmp)
}

# this function gets GWAS data using a manifest file. If a trait list is supplied
# then gets just those traits, if trait list is missing assumes that you want just
# basis
#' \code{get_gwas_data} integrate GWAS summary data with support files
#' @param trait_manifest_file character vector file path to GWAS manifest file
#' @param snp_manifest_file character vector file path to snp manifest
#' @param data_dir character vector file path to location of GWAS summary stats
#' @param filter_snps_by_manifest boolean - whether to prefilter the by snp manifest.
#' This should be true if you wish to take a subset of SNPs
#' @return data.table object
#' @export

get_gwas_data <- function(trait_manifest_file,snp_manifest_file,data_dir,filter_snps_by_manifest=FALSE){
  #if(missing(trait_list)){
  #  man<-fread(manifest_file)[basis_trait==1 & include=='Y',]
  #}else{
  #  man<-fread(manifest_file)[trait %in% trait_list & include=='Y',]
  #}
  if(!file.exists(trait_manifest_file))
    stop(sprintf("Cannot find trait manifest file %s",trait_manifest_file))
  man.DT <- fread(trait_manifest_file)
  if(nrow(man.DT)==0)
    stop(sprintf("Cannot find any traits in manifest %s for %s",manifest_file,paste(trait_list,collapse=',')))
  man.DT[,file:=file.path(data_dir,file)]
  ret<-rbindlist(lapply(1:nrow(man.DT),function(i){
    message(sprintf("Processing %s",man.DT[i,]$trait))
    tDT<-fread(man.DT[i,]$file)
    tDT[,c('trait','n','n1') := man.DT[i,.(trait,cases+controls,cases)]]
  }))
  setkey(ret,pid)
  ss <- readRDS(snp_manifest_file)
  if(filter_snps_by_manifest){
    ret <- ret[pid %in% ss$pid,]
  }
  ## next add minor allele frequencies
  message("Adding reference snp manifest annotations")
  ret<-add_ref_annotations(ss,ret)
  ret
}


#' This function computes various shrinkage metrics
#' \code{compute_shrinkage_metrics} computes various shrinkage metrics
#'
#' @param DT data.table object for basis traits as returned by \code{\link{get_gwas_data}}
#' @param pi_i numeric the prior probability that a variant is causal  see \code{\link{ws_shrinkage}}
#' @return a data.table object.
#' \enumerate{
#' \item pid - unique id using chr and position (useful for merging back)
#' \item bshrink - Bayesian shrinkage based on association across all basis traits
#' \item emp_maf_se - empirically derived standard error for MAF
#' \item est_maf_se - estimated standard error for MAF
#' \item emp_shrinkage - overall shrinkage using emp_maf_se
#' \item est_shrinkage - overall shrinkage using est_maf_se
#' }
#' see also \code{\link{maf_se_estimate_sample_size}}.
#' @export

compute_shrinkage_metrics<-function(DT,pi_i=1e-4){
  message("Computing maf_se_estimated using or, sample size and p.value ")
  ss_est_maf_se.DT<-DT[,list(pid=pid,ss_emp_maf_se=maf_se_estimate_sample_size(n,p.value,or,maf)),by=pid][,list(ss_emp_maf_se=mean(abs(ss_emp_maf_se))),by=pid]
  ss_est_maf_se.DT[,recip.ss_emp_maf_se:=1/ss_emp_maf_se]
  setkey(ss_est_maf_se.DT,pid)
  message("Computing weighted pp shrinkage")
  ws.DT <- ws_shrinkage(DT,pi_i)
  setkey(ws.DT,pid)
  shrinkage.DT <- ss_est_maf_se.DT[ws.DT]
  shrinkage.DT[,shrinkage:=ws_ppi/ss_emp_maf_se,by=pid]
  setkey(shrinkage.DT,pid)
  return(shrinkage.DT[,.(pid,shrinkage)])
}

#' This function creates a trait snp matrix
#' \code{create_ts_matrix} creates a trait snp matrix that is suitable for basis creation and projection
#'
#' @param bDT data.table object for basis traits as returned by \code{\link{get_gwas_data}}
#' @param sDT data.table object of matching shrinkage estimates returned by \code{\link{compute_shrinkage_metrics}}
#' @param method scalar vector (either emp or est), emp uses empirically generated MAF SE, est uses and estimate.
#' @return a matrix.
#' @export

create_ds_matrix <- function(bDT,sDT,method){
  if(missing(method)){
    method='shrinkage'
  }
  message(sprintf("Using %s",method))
  if(method=='none')
    sDT[,none:=1]
  #vmethod = sprintf("%s_shrinkage",method)
  stmp<-sDT[,c('pid',method),with=FALSE]
  tmp<-bDT[stmp]
  tmp$metric <- tmp[[method]] * log(tmp$or)
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  tmp.mat <- as.matrix(B[,-1]) %>% t()
  colnames(tmp.mat) <- snames
  return(tmp.mat)
}

#' This function creates a basis
#' \code{create_basis} creates a trait snp matrix that is suitable for basis creation and projection
#'
#' @param gwas.DT data.table object for basis traits as returned by \code{\link{get_gwas_data}}
#' @param shrink.DT data.table object of matching shrinkage estimates returned by \code{\link{compute_shrinkage_metrics}}
#' @param apply.shrinkage boolean on whether to apply shrinkage or not [Default TRUE]
#' @return a prcomp object representing the basis
#' @export

create_basis <- function(gwas.DT,shrink.DT,apply.shrinkage=TRUE){
  if(apply.shrinkage){
    basis.mat.emp <- create_ds_matrix(gwas.DT,shrink.DT,'shrinkage')
  }else{
    message("Warning: No shrinkage will be applied!")
    basis.mat.emp <- create_ds_matrix(gwas.DT,shrink.DT,'none')
  }
  ## need to add control where beta is zero
  basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
  prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
}

#' This function projects an aligned trait onto the basis
#' \code{project_basis} projects an external trait into basis space
#'
#' @param gwas.DT data.table object representing summary statistics from external GWAS with columns as defined in \code{\link{get_gwas_data}}
#' @param shrink.DT data.table object of matching shrinkage estimates returned by \code{\link{compute_shrinkage_metrics}}
#' @param pc prcomp object returned by \code{\link{create_basis}}
#' @param traitname character label for this trait defaults to 'test_trait'
#' @param apply.shrinkage - boolean value on whether to apply shrinkage prior to projection.
#' @return a matrix of PC scores for the projection.
#' @export

project_basis <- function(gwas.DT,shrink.DT,pc,traitname='test_trait',apply.shrinkage=TRUE){
  ## if there is already a shrinkage column this will cause an issue
  tmp <- copy(gwas.DT)
  if(any(names(gwas.DT)=='shrinkage'))
    stop("gwas.DT already contains a column named 'shrinkage' please remove or rename")
  tmp <- merge(gwas.DT,shrink.DT,by='pid',all.y=TRUE)
  ##check how many SNPs missing
  #if((!is.na(tmp$shrinkage) %>% sum) < 0.95 * nrow(shrink.DT))
  #  warning("more than 5% variants are missing")
  ## check to see we need to create beta
  if(!any(names(tmp)=='beta')){
     if(any(names(tmp)=='or')){
    	tmp[,beta:=log(or)]
     }else{
        stop("GWAS input must have either an or or beta column")
     }
  }
  ## check beta
  missing<-tmp[is.na(beta) | !is.finite(beta),]$pid
  lmiss <- length(missing)
  if(lmiss){
    sprintf("%d (%.1f%%) SNPs have missing or infinite betas setting these to 0",lmiss,(lmiss/nrow(shrink.DT))*100) %>% message()
    ## where snp is missing or infinite make it zero
    tmp[missing,beta:=0]
  }
  if(apply.shrinkage){
    tmp[,metric:=shrinkage * beta]
  }else{
    warning("no shrinkage will be applied")
    tmp[,metric:=beta]
  }
  tmp[,trait:= traitname]
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  if(!identical(colnames(mat.emp),rownames(pc$rotation)))
    stop("Something wrong basis and projection matrix don't match")
  all.proj <- predict(pc,newdata=mat.emp)
  if(any(names(tmp)=='seb')){
    tmp <- tmp[,.(pid,beta,seb,p.value,shrinkage,trait)]
  }else{
    tmp <- tmp[,.(pid,beta,p.value,shrinkage,trait)]
  }
  if(lmiss>0)
    tmp <- tmp[!pid %in% missing,]
  list(proj=all.proj,data=tmp,missing=missing)
}

#' A streamlined function to project a trait onto a sparse basis
#' \code{project_sparse}
#' @param beta a vector of beta estimates
#' @param seb a vector of standard error of the beta estimates
#' @param pids a vector of primary identifiers for SNPs with the same order as beta and seb
#' @section Notes:
#' This function assumes that the following objects are defined in the current environment
#' \itemize{
#'   \item rot.pca - Matrix of rotations, usually obtained from PCA via prcomp.
#'   \item beta.centers - Vector of basis SNP beta centres, labelled by pid.
#'   \item shrinkage - Vector of basis SNP shrinkage values, labelled by pid.
#'   \item LD - Matrix of covariance between basis SNPs
#' }
#' This function assumes that the order snps in arguments is the same. Whilst missing SNPs
#' are allowed this will degrade the projection a warining is issued when more than 5% of SNPs are missing
#' @return  a data.table with the following columns#' \itemize{
#'   \item PC - principal component label
#'   \item var.proj - Variance of the projection score.
#'   \item delta - The difference between projection score and pseudo control score.
#'   \item p.overall - The p value over all components for the projection score.
#'   \item z - z score for projection score
#'   \item p - p value for projection score.
#' }
#' @export

project_sparse <- function(beta,seb,pids){
### assumes basis-sparse-13-0.999.RData has been loaded and that LD,
### rot.pca, beta.centers, shrinkage are defined in current environment
### beta = new beta, seb=se(beta), pids=snp ids in order of beta
    if(length(beta)!=length(seb) || length(beta)!=length(pids) || !length(beta))
        stop("arguments must be equal length vectors > 0")
    if(!all(pids %in% SNP.manifest$pid))
        stop("all pids must be members of sparse basis (SNP.manifest$pid)")
    if(length(pids) < 0.95 * nrow(rot.pca))
        warning("more than 5% sparse basis snps missing")
    b <- beta * shrinkage[pids] - beta.centers[pids] # * shrinkage[pids]
    proj <- b %*% rot.pca[pids,]
    v <- seb * shrinkage[pids] * rot.pca[pids,]
    var.proj  <- t(v) %*% LD[pids,pids] %*% v
    ctl <-  (-beta.centers[pids])  %*% rot.pca[pids,]
    delta <- (proj-ctl)[1,]
    chi2 <- (t(delta) %*% solve(var.proj) %*% delta)[1,1]
    ret <- data.table::data.table(PC=colnames(proj),
                      proj=proj[1,],
                      var.proj=Matrix::diag(var.proj),
                      delta=delta,
                      p.overall=pchisq(chi2,df=13,lower.tail=FALSE))
    ret$z=ret$delta/sqrt(ret$var.proj)
    ret$p=pnorm(abs(ret$z),lower.tail=FALSE) * 2
    copy(ret)
}
