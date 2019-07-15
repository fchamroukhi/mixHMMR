
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

<!-- badges: end -->

## Overview

R code for the **clustering** and **segmentation** of time series
(including with regime changes) by mixture of gaussian Hidden Markov
Models Regression (MixHMMR) and the EM algorithm, i.e functional data
clustering and segmentation.

## Installation

You can install the development version of mixHMMR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fchamroukhi/mixHMMR")
```

To build *vignettes* for examples of usage, type the command below
instead:

``` r
# install.packages("devtools")
devtools::install_github("fchamroukhi/mixHMMR", 
                         build_opts = c("--no-resave-data", "--no-manual"), 
                         build_vignettes = TRUE)
```

Use the following command to display vignettes:

``` r
browseVignettes("mixHMMR")
```

## Usage

``` r
library(mixHMMR)

data("toydataset")

K <- 3 # Number of clusters
R <- 3 # Number of regimes/states
p <- 2 # Degree of the polynomial regression
variance_type <- "heteroskedastic" # "heteroskedastic" or "homoskedastic" model

ordered_states <- TRUE
n_tries <- 1
max_iter <- 1000
init_kmeans <- TRUE
threshold <- 1e-6
verbose <- TRUE

mixhmmr <- emMixHMMR(toydataset$x, t(toydataset[,2:ncol(toydataset)]), K, R, p,
                     variance_type, ordered_states, init_kmeans, n_tries, max_iter,
                     threshold, verbose)
#> EM: Iteration : 1 || log-likelihood : -18777.2031829953
#> EM: Iteration : 2 || log-likelihood : -15126.4977686182
#> EM: Iteration : 3 || log-likelihood : -15108.0257218485
#> EM: Iteration : 4 || log-likelihood : -15099.3886708306
#> EM: Iteration : 5 || log-likelihood : -15096.2232155737
#> EM: Iteration : 6 || log-likelihood : -15094.6653710273
#> EM: Iteration : 7 || log-likelihood : -15093.146042446
#> EM: Iteration : 8 || log-likelihood : -15092.5333943885
#> EM: Iteration : 9 || log-likelihood : -15092.3574425426
#> EM: Iteration : 10 || log-likelihood : -15092.2739302088
#> EM: Iteration : 11 || log-likelihood : -15092.2246217209
#> EM: Iteration : 12 || log-likelihood : -15092.1936060862
#> EM: Iteration : 13 || log-likelihood : -15092.1737135961
#> EM: Iteration : 14 || log-likelihood : -15092.1608712011

mixhmmr$summary()
#> ------------------------
#> Fitted mixHMMR model
#> ------------------------
#> 
#> MixHMMR model with K = 3 clusters and R = 3 regimes:
#> 
#>  log-likelihood nu       AIC      BIC      ICL
#>       -15092.16 62 -15154.16 -15197.6 -15197.6
#> 
#> Clustering table (Number of curves in each clusters):
#> 
#>  1  2  3 
#> 10 10 10 
#> 
#> Mixing probabilities (cluster weights):
#>          1         2         3
#>  0.3333333 0.3333333 0.3333333
#> 
#> 
#> --------------------
#> Cluster 1 (K = 1):
#> 
#> Regression coefficients:
#> 
#>     Beta(R = 1) Beta(R = 2) Beta(R = 3)
#> 1      6.002044    13.62910    4.934443
#> X^1    9.248449   -41.99225    4.615747
#> X^2  -41.739255    47.63765   -2.884340
#> 
#> Variances:
#> 
#>  Sigma2(R = 1) Sigma2(R = 2) Sigma2(R = 3)
#>      0.9961152      1.008064      1.006024
#> 
#> --------------------
#> Cluster 2 (K = 2):
#> 
#> Regression coefficients:
#> 
#>     Beta(R = 1) Beta(R = 2) Beta(R = 3)
#> 1      5.139985     2.17167    4.380359
#> X^1   -3.125108    21.90556    1.459464
#> X^2    9.744229   -24.65557   -0.852382
#> 
#> Variances:
#> 
#>  Sigma2(R = 1) Sigma2(R = 2) Sigma2(R = 3)
#>      0.9509166      1.052753     0.9495063
#> 
#> --------------------
#> Cluster 3 (K = 3):
#> 
#> Regression coefficients:
#> 
#>     Beta(R = 1) Beta(R = 2) Beta(R = 3)
#> 1     6.8643627    64.43719    6.218654
#> X^1   1.2876344  -462.55911   -3.949075
#> X^2  -0.1413059   893.93553    1.560385
#> 
#> Variances:
#> 
#>  Sigma2(R = 1) Sigma2(R = 2) Sigma2(R = 3)
#>      0.9770486     0.8116849      1.029812

mixhmmr$plot()
```

<img src="man/figures/README-unnamed-chunk-5-1.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-2.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-3.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-4.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-5.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-6.png" style="display: block; margin: auto;" />
