
<!-- README.md is generated from README.Rmd. Please edit that file -->

# causnet

<!-- badges: start -->

<!-- badges: end -->

The goal of causNet is to …

## Installation

You can install the development version from GitHub with:

``` r
require("devtools")
install_github("USCbiostats/causnet")
```

\~You can install the released version of causnet from
[CRAN](https://CRAN.R-project.org) with:\~

``` r
install.packages("causnet")
```

## Example

``` r
library(causnet)

# simulate data
set.seed(1234)
mydata = simdat(n.var = 5)

# causnet results
links.s = sfun(mydata, alpha = 1)

links.s


netplot_jm(links.s)

```

