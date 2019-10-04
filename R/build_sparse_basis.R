## create sparse basis that will be an output from the paper

## Some data has been shared with us exclusively and cannot be redistributed.
library(cupcake)
library(magrittr)
SPARSE_BASIS_EXTDATA <- '../inst/extdata/sparse_imd_basis'
SNP_MANIFEST_FILE <- file.path(SPARSE_BASIS_EXTDATA,'support/sparse_snp_manifest.RDS')
TRAIT_MANIFEST_FILE <- file.path(SPARSE_BASIS_EXTDATA,'support/13_trait_sparse_trait_manifest.tab')
SHRINKAGE_FILE <- file.path(SPARSE_BASIS_EXTDATA,'support/13_trait_sparse_trait_manifest.tab')
GWAS_DATA_DIR <- file.path(SPARSE_BASIS_EXTDATA,'/gwas_data/13_trait_sparse_shrinkage.RDS')

gwas.DT<-get_gwas_data(TRAIT_MANIFEST_FILE,SNP_MANIFEST_FILE,GWAS_DATA_DIR)
