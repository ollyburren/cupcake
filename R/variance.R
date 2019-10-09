library(snpStats)

#' analytically compute the variance of a projection given a standard error of beta for a given trait
#' \code{compute_proj_var}
#'
#' @param gwas.DT a data.table - variants and summary statistics
#' @param shrink.DT a data.table - as returned by \code{\link{compute_shrinkage_metrics}}
#' @param man.DT a data.table - snp manifest file that maps snps to ld.block
#' @param w.DT a data.table - a data table of weights/rotations from a  basis first column is 'pid' and subsequent columns for each principal component.
#' @param ref_gt_dir scalar - path to a dir of R objects named CHR_1kg.RData containing reference GT in snpMatrix format
#' @param method scalar - shrinkage method to use (default shrinkage) but also valid is none
#' @param quiet a scalar - boolean whether to show progress messages
#' @return a scalar of variances for each principal component
#' @export

## for testing

compute_seb_proj_var <- function(gwas.DT,shrink.DT,man.DT,w.DT,ref_gt_dir,method='shrinkage',quiet=FALSE){
  if(method=='none')
    message("Warning: Method 'none' selected no weighting applied")
  gwas.DT <- merge(gwas.DT,shrink.DT[,.(pid,shrinkage=get(`method`))],by='pid')
  ## where possible compute se_beta
  ## sometimes our input will be for a quantitative trait
  if(any(names(gwas.DT)=='or')){
    gwas.DT[,seb:=abs(log(or)/qnorm(p.value/2,lower.tail=FALSE))]
  }else if(any(names(gwas.DT)=='beta')){
    gwas.DT[,seb:=abs(beta/qnorm(p.value/2,lower.tail=FALSE))]
  }else{
    message("Cannot find required summary statistics aborting")
    return()
  }
  gwas.DT[!is.finite(seb),seb:=0]
  gwas.DT[,c('chr','pos'):=tstrsplit(pid,':')]
  ## next add ld.block information
  gwas.DT <- merge(gwas.DT,man.DT[,.(pid,ld.block)],by='pid')
  ## finally add rotations
  gwas.DT <- merge(gwas.DT,w.DT,by='pid')
  s.DT <- split(gwas.DT,gwas.DT$chr)
  all.chr <- lapply(names(s.DT),function(chr){
    if(!quiet)
      message(sprintf("Processing %s",chr))
    ss.file<-file.path(REF_GT_DIR,sprintf("%s.RDS",chr))
    if(!quiet)
      message(ss.file)
    sm <- readRDS(ss.file)
    ## there are sometimes duplicates that we need to remove
    pids <- colnames(sm)
    dup.idx<-which(duplicated(pids))
    if(length(dup.idx)>0){
      if(!quiet)
        message(sprintf("Warning removing %d duplicated SNPs",length(dup.idx)))
      sm <- sm[,-dup.idx]
      pids <- pids[-dup.idx]
    }
    # by ld block
    by.ld <- split(s.DT[[chr]],s.DT[[chr]]$ld.block)
    chr.var <- lapply(by.ld,function(block){
      if(!quiet)
        message(sprintf("Processing %s",block$ld.block %>% unique))
      sm.map <- match(block$pid,pids)
      if(any(is.na(sm.map))){
        message("SNPs in manifest that don't have genotypes")
      }
      r <- ld(sm[,sm.map],sm[,sm.map],stats="R")
      if(any(is.na(r)))
        if(!quiet)
          message(sprintf("Found %s where R^2 is NA",sum(is.na(r))))
      r[is.na(r)]<-0
      # compute closest pos-def covariance matrix
      #Sigma <- as.matrix(mvs_sigma(Matrix(r)))
      pc.cols <- which(grepl("^PC[0-9]+$",names(block)))
      sapply(names(block)[pc.cols],function(pc) (tcrossprod(block[,.(tot=get(`pc`) * shrinkage * seb)]$tot) * r) %>% sum)
    })
    colSums(do.call('rbind',chr.var))
  })
  colSums(do.call('rbind',all.chr))
}

#' analytically compute the variance of a projection given a reference set of genotypes
#' \code{compute_proj_var}
#'
#' @param man.DT a data.table - manifest of variants included in the basis
#' @param w.DT a data.table - a data table of weights/rotations from a  basis first column is 'pid' and subsequent columns for each principal component.
#' @param shrink.DT a data.table - as returned by \code{\link{compute_shrinkage_metrics}}
#' @param ref_gt_dir scalar - path to a dir of R objects named CHR_1kg.RData containing reference GT in snpMatrix format
#' @param method scalar - shrinkage method to use (default ws_emp)
#' @param quiet a scalar - boolean whether to show progress messages
#' @return a scalar of variances for each principal component
#' @export

compute_proj_var <- function(man.DT,w.DT,shrink.DT,ref_gt_dir,method='shrinkage',quiet=TRUE){
  M <- merge(man.DT,w.DT,by='pid')
  M <- merge(M,shrink.DT[,.(pid,shrink=get(`method`))],by='pid')
  M <- M[,maf:=ifelse(ref_a1.af>0.5,1-ref_a1.af,ref_a1.af)]
  ## create a set of analytical se_beta_maf
  beta_se_maf <- function(f) sqrt(1/f + 1/(1-f)) * sqrt(1/2)
  M <- M[,beta_se_maf:=beta_se_maf(maf)]
  M <- M[,c('chr','pos'):=tstrsplit(pid,':')]
  s.DT <- split(M,M$chr)
  all.chr <- lapply(names(s.DT),function(chr){
    if(!quiet)
      message(sprintf("Processing %s",chr))
    ss.file<-file.path(ref_gt_dir,sprintf("%s.RDS",chr))
    message(ss.file)
    sm <- readRDS(ss.file)
    ## there are sometimes duplicates that we need to remove
    pids <- colnames(sm)
    dup.idx<-which(duplicated(pids))
    if(length(dup.idx)>0){
      if(!quiet)
        message(sprintf("Warning removing %d duplicated SNPs",length(dup.idx)))
      sm <- sm[,-dup.idx]
      pids <- pids[-dup.idx]
    }
    # by ld block
    by.ld <- split(s.DT[[chr]],s.DT[[chr]]$ld.block)
    chr.var <- lapply(by.ld,function(block){
      if(!quiet)
        message(sprintf("Processing %s",block$ld.block %>% unique))
      sm.map <- match(block$pid,pids)
      if(any(is.na(sm.map))){
        message("SNPs in manifest that don't have genotypes")
      }
      r <- ld(sm[,sm.map],sm[,sm.map],stats="R")
      if(any(is.na(r)))
        if(!quiet)
          message(sprintf("Found %s where R^2 is NA",sum(is.na(r))))
      r[is.na(r)]<-0
      # compute closest pos-def covariance matrix
      #Sigma <- as.matrix(mvs_sigma(Matrix(r)))
      pc.cols <- which(grepl("^PC[0-9]+$",names(M)))
      sapply(names(M)[pc.cols],function(pc) (tcrossprod(block[,.(tot=get(`pc`) * shrink * beta_se_maf)]$tot) * Sigma) %>% sum)
    })
    colSums(do.call('rbind',chr.var))
  })
  colSums(do.call('rbind',all.chr))
}
