
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Overview

<!-- badges: start -->

<!-- badges: end -->

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

data("simulatedtimeseries")

K <- 3 # Number of clusters
R <- 3 # Number of regimes/states
p <- 2 # Degree of the polynomial regression
variance_type <- "heteroskedastic" # "heteroskedastic" or "homoskedastic" model

ordered_states = TRUE
n_tries = 1
max_iter = 1000
init_kmeans = TRUE
threshold = 1e-6
verbose = TRUE

mixhmmr <- emMixHMMR(simulatedtimeseries$X, t(simulatedtimeseries[, 2:ncol(simulatedtimeseries)]), K, R, p, variance_type, ordered_states, n_tries, max_iter, init_kmeans, threshold, verbose)
#> EM: Iteration : 1 || log-likelihood : -27603.7693242885
#> EM: Iteration : 2 || log-likelihood : -22464.3472396589
#> EM: Iteration : 3 || log-likelihood : -22229.1696060833
#> EM: Iteration : 4 || log-likelihood : -22212.3560575002
#> EM: Iteration : 5 || log-likelihood : -22204.7382239976
#> EM: Iteration : 6 || log-likelihood : -22192.3908663086
#> EM: Iteration : 7 || log-likelihood : -22179.0984264533
#> EM: Iteration : 8 || log-likelihood : -22161.3266480818
#> EM: Iteration : 9 || log-likelihood : -22140.2069698656
#> EM: Iteration : 10 || log-likelihood : -22118.2274362038
#> EM: Iteration : 11 || log-likelihood : -22091.4388392804
#> EM: Iteration : 12 || log-likelihood : -21996.798418692
#> EM: Iteration : 13 || log-likelihood : -21866.4806541681
#> EM: Iteration : 14 || log-likelihood : -21826.8183365824
#> EM: Iteration : 15 || log-likelihood : -21804.0499569094
#> EM: Iteration : 16 || log-likelihood : -21802.7361620057
#> EM: Iteration : 17 || log-likelihood : -21802.6973948926
#> EM: Iteration : 18 || log-likelihood : -21802.6958496405

mixhmmr$summary()
#> ------------------------
#> Fitted mixHMMR model
#> ------------------------
#> 
#> MixHMMR model with K = 3 clusters and R = 3 regimes:
#> 
#>  log-likelihood nu      AIC       BIC       ICL
#>        -21802.7 62 -21864.7 -21923.97 -21923.97
#> 
#> Clustering table:
#>  1  2  3 
#> 20 15 15 
#> 
#> Mixing probabilities (cluster weights):
#>    1   2   3
#>  0.4 0.3 0.3
#> 
#> 
#> --------------------
#> Cluster 1 (K = 1):
#> 
#> Regression coefficients:
#> 
#>       Beta(R = 1)   Beta(R = 2)   Beta(R = 3)
#> 1   10.0608592226  7.154392e+00  1.099150e+01
#> X^1 -0.0041688587 -2.054854e-03 -1.539832e-02
#> X^2  0.0000585695  6.573220e-06  2.938932e-05
#> 
#> Variances:
#> 
#>  Sigma2(R = 1) Sigma2(R = 2) Sigma2(R = 3)
#>      0.9272584     0.9732739      1.075292
#> 
#> --------------------
#> Cluster 2 (K = 2):
#> 
#> Regression coefficients:
#> 
#>       Beta(R = 1)   Beta(R = 2)   Beta(R = 3)
#> 1    8.078276e+00  1.041280e+01 -1.0174961243
#> X^1 -3.665485e-03  8.028359e-03  0.0614216445
#> X^2  4.327332e-05 -2.581304e-05 -0.0001168651
#> 
#> Variances:
#> 
#>  Sigma2(R = 1) Sigma2(R = 2) Sigma2(R = 3)
#>      0.9593165       0.97573      1.031369
#> 
#> --------------------
#> Cluster 3 (K = 3):
#> 
#> Regression coefficients:
#> 
#>       Beta(R = 1)   Beta(R = 2)   Beta(R = 3)
#> 1    8.099465e+00  8.0594251376 21.7881531429
#> X^1 -4.282857e-03  0.0432660127 -0.0922833598
#> X^2  4.141703e-05 -0.0001555985  0.0001790529
#> 
#> Variances:
#> 
#>  Sigma2(R = 1) Sigma2(R = 2) Sigma2(R = 3)
#>      0.9787244      1.011285      1.089561

mixhmmr$plot()
```

<img src="man/figures/README-unnamed-chunk-5-1.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-2.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-3.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-4.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-5.png" style="display: block; margin: auto;" />
