
<!-- README.md is generated from README.Rmd. Please edit that file -->

# typGumbo

<!-- badges: start -->
<!-- badges: end -->

typGumbo is an R package used internally within the Young-Pearse
research lab at Brigham and Women’s Hospital and Harvard Medical School.
All scripts were all written by Richard V. Pearse II PhD. Most, maybe
all, are highly specialized pipelines for -omics level data analysis. As
such they may not be extensible to your specific use case. If that
changes in the future I will alter this message. I named it typGumbo
because it’s a spicy stew. If you have questions feel free to [email
me](mailto:richard.pearse@gmail.com)

## Installation

typGumbo is not on CRAN or BioConductor, and won’t be unless I make it
more generalizable, so it can only be loaded via git download.

The development version of typGumbo can be obtained from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("genejockey33000/typGumbo")
#> Skipping install of 'typGumbo' from a github remote, the SHA1 (e213e7fe) has not changed since last install.
#>   Use `force = TRUE` to force installation
```

## Example

This is a basic example which shows you how use a few of the functions:

``` r
library(typGumbo)
## Generate a sleuth Object from kallisto output
so <- run.sleuth(in.dir = "results", meta = "meta", level = "gene")

## Generate a "typ" class R object for use with other Gumbo functions
#from a sleuth object (meta data and level are pulled automatically)
obj <- make.typ(so, regress = c("RNAseqLibBatch", "diffRound"))

#from a sleuth object and regressing for sequencing round and differentiation batch
obj <- make.typ(so)

#from a matrix (measurements in rows)
obj <- make.typ(x = {inputmatrix}, level = "gene", meta = meta, regress = c("RNAseqLibBatch", "diffRound"))
```
