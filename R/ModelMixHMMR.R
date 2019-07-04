#' A Reference Class which represents a fitted MixHMMR model.
#'
#' ModelMixHMMR represents an estimated MixHMMR model.
#'
#' @field paramMixHMMR A [ParamMixHMMR][ParamMixHMMR] object. It contains the
#'   estimated values of the parameters.
#' @field statMixHMMR A [StatMixHMMR][StatMixHMMR] object. It contains all the
#'   statistics associated to the MixHMMR model.
#' @seealso [ParamMixHMMR], [StatMixHMMR]
#' @export
ModelMixHMMR <- setRefClass(
  "ModelMixHMMR",
  fields = list(
    paramMixHMMR = "ParamMixHMMR",
    statMixHMMR = "StatMixHMMR"
  ),
  methods = list(
    plot = function() {
      "Plot method."

      # yaxislim <- c(min(paramMixHMMR$fData$Y) - 2 * mean(sqrt(apply(paramMixHMMR$fData$Y, 1, var))), max(paramMixHMMR$fData$Y) + 2 * mean(sqrt(apply(paramMixHMMR$fData$Y, 1, var))))

      matplot(paramMixHMMR$fData$X, t(paramMixHMMR$fData$Y), type = "l", lty = "solid", col = "black", xlab = "x", ylab = "y(t)")
      title(main = "Original time series")


      colorsvec <- rainbow(paramMixHMMR$K)
      matplot(paramMixHMMR$fData$X, t(paramMixHMMR$fData$Y), type = "l", lty = "dotted", col = colorsvec[statMixHMMR$klas], xlab = "x", ylab = "y(t)")
      title(main = "Clustered time series")

      for (k in 1:paramMixHMMR$K) {
        matplot(paramMixHMMR$fData$X, t(paramMixHMMR$fData$Y[statMixHMMR$klas == k,]), type = "l", lty = "dotted", col = colorsvec[k], xlab = "x", ylab = "y(t)")
        title(main = sprintf("Cluster %1.1i", k))
        lines(paramMixHMMR$fData$X, statMixHMMR$smoothed[, k], lwd = 1.5)
      }
    },

    summary = function() {
      "Summary method."

      digits = getOption("digits")

      title <- paste("Fitted mixHMMR model")
      txt <- paste(rep("-", min(nchar(title) + 4, getOption("width"))), collapse = "")

      # Title
      cat(txt)
      cat("\n")
      cat(title)
      cat("\n")
      cat(txt)

      cat("\n")
      cat("\n")
      cat(
        paste0("MixHMMR model with K = ", paramMixHMMR$K, ifelse(paramMixHMMR$K > 1, " clusters", " cluster"), " and R = ", paramMixHMMR$R, ifelse(paramMixHMMR$R > 1, " regimes", " regime"), ":"))
      cat("\n")
      cat("\n")

      tab <- data.frame("log-likelihood" = statMixHMMR$loglik, "nu" = paramMixHMMR$nu,
          "AIC" = statMixHMMR$AIC, "BIC" = statMixHMMR$BIC, "ICL" = statMixHMMR$ICL1,
          row.names = "", check.names = FALSE)
      print(tab, digits = digits)

      cat("\nClustering table:")
      print(table(statMixHMMR$klas))

      cat("\nMixing probabilities (cluster weights):\n")
      pro <- data.frame(t(paramMixHMMR$alpha))
      colnames(pro) <- 1:paramMixHMMR$K
      print(pro, digits = digits, row.names = FALSE)

      cat("\n\n")

      txt <- paste(rep("-", min(nchar(title), getOption("width"))), collapse = "")

      for (k in 1:paramMixHMMR$K) {
        cat(txt)
        cat("\nCluster ", k, " (K = ", k, "):\n", sep = "")

        cat("\nRegression coefficients:\n\n")
        if (paramMixHMMR$p > 0) {
          row.names = c("1", sapply(1:paramMixHMMR$p, function(x) paste0("X^", x)))
        } else {
          row.names = "1"
        }

        if (paramMixHMMR$p > 0) {
          row.names = c("1", sapply(1:paramMixHMMR$p, function(x) paste0("X^", x)))
          betas <- data.frame(paramMixHMMR$beta[, , k], row.names = row.names)
        } else {
          row.names = "1"
          betas <- data.frame(t(paramMixHMMR$beta[, , k]), row.names = row.names)
        }
        colnames(betas) <- sapply(1:paramMixHMMR$R, function(x) paste0("Beta(R = ", x, ")"))
        print(betas, digits = digits)

        cat(paste0(ifelse(paramMixHMMR$variance_type == "homoskedastic", "\n", "\nVariances:\n\n")))
        sigma2 <- data.frame(t(paramMixHMMR$sigma2[, k]))
        if (paramMixHMMR$variance_type == "homoskedastic") {
          colnames(sigma2) <- "Sigma2"
          print(sigma2, digits = digits, row.names = FALSE)
        } else {
          colnames(sigma2) = sapply(1:paramMixHMMR$R, function(x) paste0("Sigma2(R = ", x, ")"))
          print(sigma2, digits = digits, row.names = FALSE)
        }
        cat("\n")
      }

    }
  )
)
