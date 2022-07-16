---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# CausNet

<!-- badges: start -->
<!-- badges: end -->

The goal of CausNet is to find globally optimal Bayesian networks via dynamic programming with parent set constraints.

## Installation

You can install the development version from GitHub with:

```{r installation, eval=FALSE}
require("devtools")
install_github("https://github.com/nand1155/CausNet")
```

## Example

```{r}
library(CausNet)

# simulate data
set.seed(1234)
mydata = simdat(300,5,1)
# run Causnet

links.s = sfun(mydata,  surdata=NULL, scoreFn = "bic", pheno = FALSE, fdr = FALSE, alpha = 0.6, alpha1 = NULL, alpha2 = NULL, pp = NULL, multBNs = TRUE)
netplot_jm(links.s[[1]]) # if multBNs = TRUE
netplot_jm(links.s) # if multBNs = False
```



