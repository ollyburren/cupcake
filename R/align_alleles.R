## this contains code to align alleles between studies

#' \code{comp} get complement for a set of variants
#'
#' @param cv a vector of variants
#' @return a vector of complement variants

comp<-function(cv){
  tmp<-cv
  a<-c('A','G','C','T')
  b<-c('T','C','G','A')
  idx<-split(1:length(cv),cv)
  for(i in seq_along(a)){
    cv[idx[[a[i]]]]<-b[i]
  }
  cv
}

#' \code{align_alleles} compute variance weighted euclidean distance between two PC loading vectors a and b
#'
#' @param gwas.DT a data.table of GWAS data see \code{\link{get_gwas_data}}
#' @param ref.DT a data.table of reference alleles (see details)
#' @param check boolean scalar whether to add column of transformations for checking
#' @return gwas.DT with allele aligned DT
#' @export

#' @details this routine lines up on reference `a1` this is arbitrary but as long as consistency is maintained
#' then the basis can be correctly created


align_alleles<-function(gwas.DT,ref.DT,check=FALSE){
  ## first make it so that a1 and a2 on gwas.DT are converted to risk.allele and other.allele
  gwas.DT[,c('risk.allele','other.allele'):=list(a1,a2)]
  gwas.DT[or<1,c('risk.allele','other.allele','or'):=list(a2,a1,signif(1/or,digits=3))]
  ## next merge in ref information
  gwas.DT <- gwas.DT[,.(id,chr,position,p.val,or,risk.allele,other.allele)]
  ref.DT[,pid:=paste(chr,position,sep=':')]
  gwas.DT[,pid:=paste(chr,position,sep=':')]
  setkey(ref.DT,pid)
  setkey(gwas.DT,pid)
  gwas.DT<-gwas.DT[ref.DT]
  flip.idx<-with(gwas.DT,which(risk.allele == a2 & other.allele ==a1))
  gwas.DT[flip.idx,c('risk.allele','other.allele','or','check'):=list(a1,a2,signif(1/or,digits=3),'flip')]
  non.match.idx<-which(gwas.DT$risk.allele != gwas.DT$a1)
	flip.comp.idx<-with(gwas.DT,which(comp(risk.allele) == a2 & comp(other.allele) ==a1))
  actual.flip.comp.idx<-intersect(non.match.idx,flip.comp.idx)
  gwas.DT[actual.flip.comp.idx,c('risk.allele','other.allele','or','check'):=list(a1,a2,signif(1/or,digits=3),'comp.flip')]
  rev.comp.idx<-with(gwas.DT,which(comp(risk.allele) == a1 & comp(other.allele) ==a2))
  if(!check){
    gwas.DT[rev.comp.idx,c('risk.allele','other.allele','check'):=list(a1,a2,'comp')][,.(id,chr,position,p.val,or)]
  }else{
    gwas.DT[is.na(check),check:='none']
    gwas.DT[rev.comp.idx,c('risk.allele','other.allele','check'):=list(a1,a2,'comp')][,.(id,chr,position,p.val,or,check)]
  }
  return(gwas.DT)
}
