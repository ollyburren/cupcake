library(snpStats)

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
#' are allowed this will degrade and projection when more than 5% of SNPs are missing a warnings
#' will be generated
#' @returns  a data.table with the following columns:
#' \itemize{
#'   \item PC - principal component label
#'   \item var.proj - Variance of the projection score.
#'   \item delta - The difference between projection score and pseudo control score.
#'   \item p.overall - The p value over all components for the projection score.
#'   \item z - z score for projection score
#'   \item z - p value for projection score.
#' }


project_sparse <- function(beta,seb,pids) {
###' assumes basis-sparse-13-0.999.RData has been loaded and that LD,
###' rot.pca, beta.centers, shrinkage are defined in current environment
###' beta = new beta, seb=se(beta), pids=snp ids in order of beta
    require(Matrix)
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
    ret <- data.table(PC=colnames(proj),
                      proj=proj[1,],
                      var.proj=diag(var.proj),
                      delta=delta,
                      p.overall=pchisq(chi2,df=13,lower.tail=FALSE))
    ret$z=ret$delta/sqrt(ret$var.proj)
    ret$p=pnorm(abs(ret$z),lower.tail=FALSE) * 2
    copy(ret)
}


#' An internal helper function to get obtain standard error of beta estimates from a GWAS.DT data.frame object
#' \code{get_seb}
#' @param gwas.DT a data.table - variants and summary statistics
#' @return a data.table with seb column added if not already present

get_seb <- function(gwas.DT){
  if(!any(names(gwas.DT)=='seb')){
    if(any(names(gwas.DT)=='or')){
      gwas.DT[,seb:=abs(log(or)/qnorm(p.value/2,lower.tail=FALSE))]
    }else if(any(names(gwas.DT)=='beta')){
      gwas.DT[,seb:=abs(beta/qnorm(p.value/2,lower.tail=FALSE))]
    }else{
      message("Cannot find required summary statistics aborting")
      return()
    }
  }
  gwas.DT[!is.finite(seb),seb:=0]
}

#' analytically compute the variance of a projection given a standard error of beta for a given trait for sparse basis
#' \code{compute_proj_var}
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
  if(method=='none')
    message("Warning: Method 'none' selected no weighting applied")
  gwas.DT <- merge(gwas.DT,shrink.DT[,.(pid,shrinkage=get(`method`))],by='pid')
  gwas.DT <- get_seb(gwas.DT)
  gwas.DT[,c('chr','pos'):=tstrsplit(pid,':')]
  gwas.DT <- merge(gwas.DT,w.DT,by='pid')
  pids <- colnames(sm)
  sm.map <- match(gwas.DT$pid,pids)
  if(any(is.na(sm.map))){
    message("SNPs in manifest that don't have genotypes")
  }
  r <- ld(sm[,sm.map],sm[,sm.map],stats="R")
  if(any(is.na(r)))
    if(!quiet)
      message(sprintf("Found %s where R^2 is NA",sum(is.na(r))))
  r[is.na(r)]<-0
  pc.cols <- which(grepl("^PC[0-9]+$",names(gwas.DT)))
  sapply(names(gwas.DT)[pc.cols],function(pc) (tcrossprod(gwas.DT[,.(tot=get(`pc`) * shrinkage * seb)]$tot) * r) %>% sum)
}



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
  # if(!any(names(gwas.DT)=='seb')){
  #   if(any(names(gwas.DT)=='or')){
  #     gwas.DT[,seb:=abs(log(or)/qnorm(p.value/2,lower.tail=FALSE))]
  #   }else if(any(names(gwas.DT)=='beta')){
  #     gwas.DT[,seb:=abs(beta/qnorm(p.value/2,lower.tail=FALSE))]
  #   }else{
  #     message("Cannot find required summary statistics aborting")
  #     return()
  #   }
  # }
  # gwas.DT[!is.finite(seb),seb:=0]
  gwas.DT <- get_seb(gwas.DT)
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
