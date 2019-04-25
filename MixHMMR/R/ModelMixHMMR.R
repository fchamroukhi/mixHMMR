source("R/FData.R")
source("R/enums.R")

ModelMixHMMR <- setRefClass(
  "ModelMixHMMR",
  contains = "FData",
  # Define the fields
  fields = list(
    K = "numeric", # Number of clusters
    R = "numeric", # Number of regimes (HMM states)
    p = "numeric", # dimension of beta (order of polynomial regression)
    variance_type = "numeric",
    nu = "numeric" # degree of freedom
  )
)

ModelMixHMMR <- function(fData, K, R, p, variance_type) {
  if (variance_type == variance_types$homoskedastic) {
    nu <<- (K - 1) + K * ((R - 1) + R * (R - 1) + R * (p + 1) + 1)
  }
  else{
    nu <<- (K - 1) + K * ((R - 1) + R * (R - 1) + R * (p + 1) + R)
  }

  new(
    "ModelMixHMMR",
    Y = fData$Y,
    X = fData$X,
    m = fData$m,
    n = fData$n,
    vecY = fData$vecY,
    K = K,
    R = R,
    p = p,
    variance_type = variance_type,
    nu = nu
  )
}
