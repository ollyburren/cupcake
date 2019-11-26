library(snpStats)

#' An internal helper function to get obtain standard error of beta estimates from a GWAS.DT data.frame object
#' \code{get_seb}
#' @param gwas.DT a data.table - variants and summary statistics
#' @return a data.table with seb column added if not already present

get_seb <- function(gwas.DT){
  seb <- or <- p.value <- NULL
  if(!any(names(gwas.DT)=='seb')){
    if(any(names(gwas.DT)=='or')){
      gwas.DT[,seb:=abs(log(or)/stats::qnorm(p.value/2,lower.tail=FALSE))]
    }else if(any(names(gwas.DT)=='beta')){
      gwas.DT[,seb:=abs(beta/stats::qnorm(p.value/2,lower.tail=FALSE))]
    }else{
      message("Cannot find required summary statistics aborting")
      return()
    }
  }
  gwas.DT[!is.finite(seb),seb:=0]
}

#' analytically compute the variance of a projection given a standard error of beta for a given trait for sparse basis
#' \code{compute_seb_proj_var_sparse}
#'
#' @param gwas.DT a data.table - variants and summary statistics
#' @param shrink.DT a data.table - as returned by \code{\link{compute_shrinkage_metrics}}
#' @param w.DT a data.table - a data table of weights/rotations from a  basis first column is 'pid' and subsequent columns for each principal component.
#' @param sm snpMatrix object - a snp matrix object of reference genotypes for
#' @param method scalar - shrinkage method to use (default shrinkage) but also valid is none
#' @param quiet a scalar - boolean whether to show progress messages
#' @return a scalar of variances for each principal component
#' @export

compute_seb_proj_var_sparse <- function(gwas.DT,shrink.DT,w.DT,sm,method='shrinkage',quiet=FALSE){
  shrinkage <- seb <- pid <- ld <- NULL
  if(method=='none')
    message("Warning: Method 'none' selected no weighting applied")
  gwas.DT <- merge(gwas.DT,shrink.DT[,list(pid,shrinkage=get(`method`))],by='pid')
  gwas.DT <- get_seb(gwas.DT)
  gwas.DT[,c('chr','pos'):=tstrsplit(pid,':')]
  gwas.DT <- merge(gwas.DT,w.DT,by='pid')
  pids <- colnames(sm)
  sm.map <- match(gwas.DT$pid,pids)
  if(any(is.na(sm.map))){
    message("SNPs in manifest that don't have genotypes")
  }
  r <- snpStats::ld(sm[,sm.map],sm[,sm.map],stats="R")
  if(any(is.na(r)))
    if(!quiet)
      message(sprintf("Found %s where R^2 is NA",sum(is.na(r))))
  r[is.na(r)]<-0
  pc.cols <- which(grepl("^PC[0-9]+$",names(gwas.DT)))
  sapply(names(gwas.DT)[pc.cols],function(pc) (tcrossprod(gwas.DT[,list(tot=get(`pc`) * shrinkage * seb)]$tot) * r) %>% sum)
}



#' analytically compute the variance of a projection given a standard error of beta for a given trait
#' \code{compute_seb_proj_var}
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


compute_seb_proj_var <- function(gwas.DT,shrink.DT,man.DT,w.DT,ref_gt_dir,method='shrinkage',quiet=FALSE){
  shrinkage <- pid <- ld.block <- seb <- NULL
  if(method=='none')
    message("Warning: Method 'none' selected no weighting applied")
  gwas.DT <- merge(gwas.DT,shrink.DT[,list(pid,shrinkage=get(`method`))],by='pid')
  gwas.DT <- get_seb(gwas.DT)
  gwas.DT[,c('chr','pos'):=tstrsplit(pid,':')]
  ## next add ld.block information
  gwas.DT <- merge(gwas.DT,man.DT[,list(pid,ld.block)],by='pid')
  ## finally add rotations
  gwas.DT <- merge(gwas.DT,w.DT,by='pid')
  s.DT <- split(gwas.DT,gwas.DT$chr)
  all.chr <- lapply(names(s.DT),function(chr){
    if(!quiet)
      message(sprintf("Processing %s",chr))
    ss.file<-file.path(ref_gt_dir,sprintf("%s.RDS",chr))
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
      r <- snpStats::ld(sm[,sm.map],sm[,sm.map],stats="R")
      if(any(is.na(r)))
        if(!quiet)
          message(sprintf("Found %s where R^2 is NA",sum(is.na(r))))
      r[is.na(r)]<-0
      # compute closest pos-def covariance matrix
      #Sigma <- as.matrix(mvs_sigma(Matrix(r)))
      pc.cols <- which(grepl("^PC[0-9]+$",names(block)))
      sapply(names(block)[pc.cols],function(pc) (tcrossprod(block[,list(tot=get(`pc`) * shrinkage * seb)]$tot) * r) %>% sum)
    })
    colSums(do.call('rbind',chr.var))
  })
  colSums(do.call('rbind',all.chr))
}
