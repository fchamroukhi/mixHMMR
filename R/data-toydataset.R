#' A dataset composed of simulated time series with regime changes.
#'
#' A dataset composed of 30 simulated time series with regime changes.
#'
#' @format A data frame with 350 rows and 31 variables:
#' \describe{
#'   \item{x}{The covariate variable which is the time in that case.}
#'   \item{y1 to y10}{10 times series having the same shape for which a
#'   normally distributed random noise has been added.}
#'  \item{y11 to y20}{10 time series generated as follows:
#'     \itemize{
#'       \item  First regime: 120 values of Normally distributed random numbers
#'       with mean 5.
#'       \item Second regime: 70 values of Normally distributed random numbers
#'       with mean 7.
#'       \item Third regime: 160 values of Normally distributed random numbers
#'       with mean 5.
#'     }
#'   }
#'  \item{y21 to y30}{10 time series generated as follows:
#'     \itemize{
#'       \item  First regime: 80 values of Normally distributed random numbers
#'       with mean 7.
#'       \item Second regime: 130 values of Normally distributed random numbers
#'       with mean 5.
#'       \item Third regime: 140 values of Normally distributed random numbers
#'       with mean 4.
#'     }
#'   }
#' }
#'
"toydataset"
