#' @export
ModelMixHMMR <- setRefClass(
  "ModelMixHMMR",
  fields = list(
    paramMixHMMR = "ParamMixHMMR",
    statMixHMMR = "StatMixHMMR"
  ),
  methods = list(
    plot = function() {

      # yaxislim <- c(min(paramMixHMMR$fData$Y) - 2 * mean(sqrt(apply(paramMixHMMR$fData$Y, 1, var))), max(paramMixHMMR$fData$Y) + 2 * mean(sqrt(apply(paramMixHMMR$fData$Y, 1, var))))

      matplot(t(paramMixHMMR$fData$Y), type = "l", lty = "solid", col = "black", xlab = "x", ylab = "y(t)")
      title(main = "Original time series")


      colorsvec <- rainbow(paramMixHMMR$K)
      matplot(t(paramMixHMMR$fData$Y), type = "l", lty = "dotted", col = colorsvec[statMixHMMR$klas], xlab = "x", ylab = "y(t)")
      title(main = "Clustered time series")

      for (k in 1:paramMixHMMR$K) {
        matplot(t(paramMixHMMR$fData$Y[statMixHMMR$klas == k, ]), type = "l", lty = "dotted", col = colorsvec[k], xlab = "x", ylab = "y(t)")
        title(main = sprintf("Cluster %1.1i", k))
        lines(statMixHMMR$smoothed[, k], lwd = 1.5)
      }
    }
  )
)
