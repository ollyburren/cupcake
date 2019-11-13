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

## this contains code to align alleles between studies

#' \code{ambiguous} get ambiguous variants
#'
#' @param a1 a vector of alleles AGCT
#' @param a2 a vector of alleles AGCT
#' @return a boolean vector

ambiguous<-function(a1,a2){
  tmp<-data.table(a1=a1,a2=a2,amb=FALSE)
  tmp[(a1=='A' & a2=='T') | (a1=='T' & a2=='A'),amb:=TRUE]
  tmp[(a1=='G' & a2=='C') | (a1=='C' & a2=='G'),amb:=TRUE]
  return(tmp$amb)
}

#' \code{align_alleles} compute variance weighted euclidean distance between two PC loading vectors a and b
#'
#' @param gwas.DT a data.table of GWAS data see \code{\link{get_gwas_data}}
#' @param ref.DT a data.table of reference alleles (see details)
#' @param check boolean scalar whether to add column of transformations for checking
#' @return gwas.DT with allele aligned DT

#' @details this routine lines up on reference `a1` this is arbitrary but as long as consistency is maintained
#' then the basis can be correctly created


align_alleles<-function(gwas.DT,ref.DT,check=FALSE){
  ## first make it so that a1 and a2 on gwas.DT are converted to risk.allele and other.allele
  #gwas.DT[,c('risk.allele','other.allele'):=list(a1,a2)]
  gwas.DT[,c('risk.allele','other.allele','ambig'):=list(a2,a1,ambiguous(a1,a2))]
  #gwas.DT[or>1,c('risk.allele','other.allele','or'):=list(a2,a1,1/or)]
  gwas.DT[or<1,c('risk.allele','other.allele','or'):=list(a1,a2,1/or)]
  ## next merge in ref information
  gwas.DT <- gwas.DT[,.(id,chr,position,p.val,or,risk.allele,other.allele,ambig)]
  ref.DT[,pid:=paste(chr,position,sep=':')]
  gwas.DT[,pid:=paste(chr,position,sep=':')]
  setkey(ref.DT,pid)
  setkey(gwas.DT,pid)
  gwas.DT<-gwas.DT[ref.DT]
  flip.idx<-with(gwas.DT,which(risk.allele == a1 & other.allele ==a2))
  gwas.DT[flip.idx,c('risk.allele','other.allele','or','check'):=list(a1,a2,1/or,'flip')]
  non.match.idx<-which(gwas.DT$risk.allele != gwas.DT$a1)
	flip.comp.idx<-with(gwas.DT,which(comp(risk.allele) == a1 & comp(other.allele) ==a2))
  actual.flip.comp.idx<-intersect(non.match.idx,flip.comp.idx)
  gwas.DT[actual.flip.comp.idx,c('risk.allele','other.allele','or','check'):=list(a1,a2,1/or,'comp.flip')]
  rev.comp.idx<-with(gwas.DT,which(comp(risk.allele) == a2 & comp(other.allele) ==a1))
  gwas.DT <- gwas.DT[rev.comp.idx,c('risk.allele','other.allele','check'):=list(a1,a2,'comp')]
  if(!check){
    return(gwas.DT[,.(id,chr,position,p.val,or)])
  }else{
    gwas.DT[is.na(check),check:='none']
    return(gwas.DT[,.(id,chr,position,p.val,or,check,risk.allele,other.allele,a1,a2,ambig)])
  }
}

## align OR
#' \code{align_allele} given a data.table of pid,a1,a2 get a list of SNPs where alleles are flipped wrt to basis reference
#' @param DT - a data.table object - columns chr,position,a1 and a2
#' @param ref.DT - a data.table of reference alleles (see details)
#' @return an index wrt to DT of OR that require flipping
#' @export

flip_allele<-function(DT,ref.DT){
   ref.DT[,pid:=paste(chr,position,sep=':')]
   DT[,pid:=paste(chr,position,sep=':')]
   setkey(ref.DT,pid)
   setkey(DT,pid)
   # a1 is from dataset i.a1 is from reference
   gwas.DT<-DT[ref.DT][,.(pid,a1,a2,i.a1,i.a2)]
   no.change.pid<-gwas.DT[a1==i.a1 & a2==i.a2,]$pid
   flip.pid <- gwas.DT[a2==i.a1 & a1==i.a2,]$pid
   rev.no.change.pid <- gwas.DT[comp(a1)==i.a1 & comp(a2)==i.a2,]$pid
   rev.flip.pid <- gwas.DT[comp(a2)==i.a1 & comp(a1)==i.a2,]$pid
   nc.pid <- union(no.change.pid,rev.no.change.pid)
   c.pid <- union(flip.pid,rev.flip.pid)
   if(length(intersect(nc.pid,c.pid))>0)
    stop("Can't work out allele flipping")
   #check to see if there are any issues
   if(any(is.na(gwas.DT$pid)))
    stop("Missing some SNPs when aligning with basis reference")
   return(which(DT$pid %in% c.pid))
 }
