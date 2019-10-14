# cupcake

This package contains routines to perform principal component analysis (PCA) using genome-wide association study summary statistics (GWAS). The resultant **basis** is a lower-dimensional representation of input traits, that summarises their genetic relationships. The basis' main utility is the projection of other traits, especially those for which sample sizes are modest and convential GWAS analyses are challenging. These projected traits can be used for the biological characterisation of exisiting components or as a vehicle for investigating the genetic architecture of clinically related traits.

# Installation

```
library(devtools)
install_github('ollyburren/cupcake',build_vignettes=TRUE)
vignette('create-basis',package='cupcake')
````
