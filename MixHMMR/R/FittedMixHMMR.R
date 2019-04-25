FittedMixHMMR <- setRefClass(
  "FittedMixHMMR",
  fields = list(
    modelMixHMMR = "ModelMixHMMR",
    paramMixHMMR = "ParamMixHMMR",
    statMixHMMR = "StatMixHMMR"
  ),
  methods = list(
    plot = function() {

      matplot(t(modelMixHMMR$Y), col = "black", type = "l", lty = "solid", xlab = "Time", ylab = "y(t)")
      title(main = "Original time series")


      couleurs = rainbow(modelMixHMMR$K)
      matplot(t(modelMixHMMR$Y), col = couleurs[statMixHMMR$klas], type = "l", lty = "dotted", xlab = "Time", ylab = "y(t)")
      title(main = "Clustered time series")


    #   yaxislim <- c(mean(modelHMMR$Y) - 2 * sd(modelHMMR$Y), mean(modelHMMR$Y) + 2 * sd(modelHMMR$Y))
    #
    #   par(mfrow = c(2, 1))
    #
    #   # Predicted time series and predicted regime probabilities
    #   plot.default(modelHMMR$Y, type = "l", ylab = "y", xlab = "", ylim = yaxislim)
    #   lines(statHMMR$predicted, type = "l", lwd = 2, col = "red")
    #   title(main = "Original and predicted HMMR time series")
    #
    #   # Prediction probabilities of the hidden process (segmentation)
    #   colors <- rainbow(modelHMMR$K)
    #   plot.default(statHMMR$predict_prob[, 1], col = colors[1], type = "l", lwd = 1.5, main = "Prediction probabilities", ylab = "Prob")
    #   for (k in 2:modelHMMR$K) {
    #     lines(statHMMR$predict_prob[, k], col = colors[k])
    #   }
    #
    #   # Filtered time series and filtering regime probabilities
    #   par(mfrow = c(2, 1))
    #   plot.default(modelHMMR$Y, type = "l", ylab = "y", xlab = "", ylim = yaxislim)#black
    #   title(main = "Original and filtered HMMR time series")
    #   lines(statHMMR$filtered, col = "red", lwd = 2)
    #
    #   # Filtering probabilities of the hidden process (segmentation)
    #   plot.default(statHMMR$filter_prob[, 1], col = colors[1], type = "l", lwd = 1.5, main = "Filtering probabilities", xlab = "t", ylab = "Prob")
    #   for (k in 2:modelHMMR$K) {
    #     lines(statHMMR$filter_prob[, k], col = colors[k], type = "l", lwd = 1.5) #Post Probs: Pr(Z_{t}=k|y_1,\ldots,y_t)
    #   }
    #
    #   # Data, regressors, and segmentation
    #   par(mfrow = c(2, 1))
    #   plot.default(modelHMMR$Y, type = "l", ylab = "y", xlab = "", ylim = yaxislim)
    #   title(main = "Time series, HMMR regimes, and smoothing probabilites")
    #   for (k in 1:modelHMMR$K) {
    #     model_k <- statHMMR$regressors[, k]
    #     #prob_model_k = HMMR$param$piik[,k]
    #
    #     index <- statHMMR$klas == k
    #     active_model_k <- model_k[index]#prob_model_k >= prob);
    #     active_period_model_k <- seq(1:modelHMMR$m)[index]#prob_model_k >= prob);
    #
    #     if (length(active_model_k) != 0) {
    #       lines(model_k, col = colors[k], lty = "dotted", lwd = 1)
    #       lines(active_period_model_k, active_model_k, col = colors[k], type = "l", lwd = 3)
    #     }
    #   }
    #
    #   # Probablities of the hidden process (segmentation)
    #   plot.default(statHMMR$tau_tk[, 1], main = "Smoothing probabilities", xlab = "", ylab = "Prob", col = colors[1], type = "l", lwd = 1.5)
    #   if (modelHMMR$K > 1) {
    #     for (k in 2:modelHMMR$K) {
    #       lines(statHMMR$tau_tk[, k], col = colors[k], type = "l", lwd = 1.5)
    #     }
    #   }
    #   #Post Probs: Pr(Z_{t}=k|y_1,\ldots,y_n)
    #
    #   ## data, regression model, and segmentation
    #   par(mfrow = c(2, 1))
    #   plot.default(modelHMMR$Y, type = "l", ylab = "y", xlab = "", ylim = yaxislim)
    #   title(main = "Original and smoothed HMMR time series, and segmentation")
    #   lines(statHMMR$smoothed, col = "red", lwd = 2)
    #
    #   # Transition time points
    #   tk <- which(diff(statHMMR$klas) != 0)
    #   for (i in 1:length(tk)) {
    #     abline(v = tk[i], lty = "dotted", lwd = 2, col = "red")
    #   }
    #
    #   # Probablities of the hidden process (segmentation)
    #   plot.default(statHMMR$klas, type = "l", lwd = 2, col = "red", xlab = "", ylab = "Estimated class labels")
    }
  )
)

FittedMixHMMR <- function(modelMixHMMR, paramMixHMMR, statMixHMMR) {
  new("FittedMixHMMR", modelMixHMMR = modelMixHMMR, paramMixHMMR = paramMixHMMR, statMixHMMR = statMixHMMR)
}
