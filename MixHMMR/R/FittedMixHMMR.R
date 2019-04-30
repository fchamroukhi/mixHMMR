FittedMixHMMR <- setRefClass(
  "FittedMixHMMR",
  fields = list(
    modelMixHMMR = "ModelMixHMMR",
    paramMixHMMR = "ParamMixHMMR",
    statMixHMMR = "StatMixHMMR"
  ),
  methods = list(
    plot = function() {

      # yaxislim <- c(min(modelMixHMMR$Y) - 2 * mean(sqrt(apply(modelMixHMMR$Y, 1, var))), max(modelMixHMMR$Y) + 2 * mean(sqrt(apply(modelMixHMMR$Y, 1, var))))

      matplot(t(modelMixHMMR$Y), type = "l", lty = "solid", col = "black", xlab = "Time", ylab = "y(t)")
      title(main = "Original time series")


      colorsvec <- rainbow(modelMixHMMR$K)
      matplot(t(modelMixHMMR$Y), type = "l", lty = "dotted", col = colorsvec[statMixHMMR$klas], xlab = "Time", ylab = "y(t)")
      title(main = "Clustered time series")

      for (k in 1:modelMixHMMR$K) {
        matplot(t(modelMixHMMR$Y[statMixHMMR$klas == k, ]), type = "l", lty = "dotted", col = colorsvec[k], xlab = "Time", ylab = "y(t)")
        title(main = sprintf("Cluster %1.1i", k))
        lines(statMixHMMR$smoothed[, k], lwd = 1.5)
      }
    }
  )
)

FittedMixHMMR <- function(modelMixHMMR, paramMixHMMR, statMixHMMR) {
  new("FittedMixHMMR", modelMixHMMR = modelMixHMMR, paramMixHMMR = paramMixHMMR, statMixHMMR = statMixHMMR)
}
