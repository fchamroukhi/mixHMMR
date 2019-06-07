#' @export
StatMixHMMR <- setRefClass(
  "StatMixHMMR",
  fields = list(
    tau_ik = "matrix",
    gamma_ikjr = "array",
    log_w_k_fyi = "matrix",
    exp_num_trans = "array",
    exp_num_trans_from_l = "array",
    loglik = "numeric",
    stored_loglik = "list",
    cputime = "numeric",
    klas = "matrix", # klas: [nx1 double]
    z_ik = "matrix", # z_ik: [nxK]
    smoothed = "matrix",
    mean_curves = "array",
    BIC = "numeric",
    AIC = "numeric",
    ICL1 = "numeric"
  ),
  methods = list(

    initialize = function(paramMixHMMR = ParamMixHMMR(fData = FData(numeric(1), matrix(1)), K = 2, R = 1, p = 2, variance_type = 1)) {

      tau_ik <<- matrix(NA, paramMixHMMR$fData$n, paramMixHMMR$K)
      gamma_ikjr <<- array(NA, dim = c(paramMixHMMR$fData$n * paramMixHMMR$fData$m, paramMixHMMR$R, paramMixHMMR$K))
      log_w_k_fyi <<- matrix(NA, paramMixHMMR$fData$n, paramMixHMMR$K)
      exp_num_trans <<- array(NA, dim = c(paramMixHMMR$R, paramMixHMMR$R, paramMixHMMR$fData$n, paramMixHMMR$K))
      exp_num_trans_from_l <<- array(NA, dim = c(paramMixHMMR$R, paramMixHMMR$fData$n, paramMixHMMR$K))
      loglik <<- -Inf
      stored_loglik <<- list()
      cputime <<- Inf
      klas <<- matrix(NA, paramMixHMMR$fData$n, 1) # klas: [nx1 double]
      z_ik <<- matrix(NA, paramMixHMMR$fData$n, paramMixHMMR$K) # z_ik: [nxK]
      smoothed <<- matrix(NA, paramMixHMMR$fData$m, paramMixHMMR$K)
      mean_curves <<- array(NA, dim = c(paramMixHMMR$fData$m, paramMixHMMR$R, paramMixHMMR$K))
      BIC <<- -Inf
      AIC <<- -Inf
      ICL1 <<- -Inf

    },

    MAP = function() {
      N <- nrow(tau_ik)
      K <- ncol(tau_ik)
      ikmax <- max.col(tau_ik)
      ikmax <- matrix(ikmax, ncol = 1)
      z_ik <<- ikmax %*% ones(1, K) == ones(N, 1) %*% (1:K) # partition_MAP
      klas <<- ones(N, 1)
      for (k in 1:K) {
        klas[z_ik[, k] == 1] <<- k
      }
    },

    computeStats = function(paramMixHMMR, cputime_total) {

      cputime <<- mean(cputime_total)

      for (k in 1:paramMixHMMR$K) {

        betakr <- paramMixHMMR$beta_kr[, , k]
        weighted_segments <- apply(gamma_ikjr[, , k] * (repmat(paramMixHMMR$phi, paramMixHMMR$fData$n, 1) %*% betakr), 1, sum)
        dim(weighted_segments) <- c(paramMixHMMR$fData$m, paramMixHMMR$fData$n)
        weighted_clusters <- (matrix(1, paramMixHMMR$fData$m, 1) %*% t(tau_ik[, k])) * weighted_segments
        smoothed[, k] <<- apply(weighted_clusters, 1, sum) / sum(tau_ik[, k])
      }

      # BIC AIC et ICL*
      BIC <<- loglik - (paramMixHMMR$nu * log(paramMixHMMR$fData$n) / 2)
      AIC <<- loglik - paramMixHMMR$nu
      # ICL*
      # Compute the comp-log-lik
      cik_log_w_k_fyi <- z_ik * log_w_k_fyi
      comp_loglik <- sum(cik_log_w_k_fyi)
      ICL1 <<- comp_loglik - paramMixHMMR$nu * log(paramMixHMMR$fData$n) / 2 #n*m/2!

    },

    EStep = function(paramMixHMMR) {

      exp_num_trans_ck  <- array(0, dim = c(paramMixHMMR$R, paramMixHMMR$R, paramMixHMMR$fData$n))
      exp_num_trans_from_l_ck <- matrix(0, paramMixHMMR$R, paramMixHMMR$fData$n)

      for (k in 1:paramMixHMMR$K) {
        # Run a hmm for each sequence
        log_fkr_yij <- matrix(0, paramMixHMMR$R, paramMixHMMR$fData$m)
        fkr_yij <- matrix(0, paramMixHMMR$R, paramMixHMMR$fData$m)

        Li <- matrix(0, paramMixHMMR$fData$n, 1) # To store the loglik for each example (curve)

        beta_kr <- paramMixHMMR$beta_kr[, , k]
        num_log_post_prob <- matrix(0, paramMixHMMR$fData$n, paramMixHMMR$K)

        for (i in 1:paramMixHMMR$fData$n) {
          y_i <- paramMixHMMR$fData$Y[i, ]

          for (r in 1:paramMixHMMR$R) {
            betakr <- beta_kr[, r]

            if (paramMixHMMR$variance_type == variance_types$homoskedastic) {
              sigma_kr <- paramMixHMMR$sigma_kr[, k]
              sk <- sigma_kr
            }
            else{
              sigma_kr <- paramMixHMMR$sigma_kr[, k]
              sk <- sigma_kr[r]
            }
            z <- ((y_i - t(paramMixHMMR$phi %*% betakr)) ^ 2) / sk

            log_fkr_yij[r, ] <- -0.5 * matrix(1, 1, paramMixHMMR$fData$m) * (log(2 * pi) + log(sk)) - 0.5 * z # Log pdf yij | c_i = k et z_i = r

          }

          log_fkr_yij <- pmin(log_fkr_yij, log(.Machine$double.xmax))
          log_fkr_yij <- pmax(log_fkr_yij, log(.Machine$double.xmin))
          fkr_yij <- exp(log_fkr_yij)


          # Forwards backwards ( calcul de logProb(Yi)...)
          fb <- forwardsBackwards(paramMixHMMR$pi_k[, k], paramMixHMMR$A_k[, , k], fkr_yij)

          gamma_ik <- fb$tau_tk
          xi_ik <- fb$xi_tk
          fwd_ik <- fb$alpha_tk
          backw_ik <- fb$beta_tk
          loglik_i <- fb$loglik

          Li[i] <- loglik_i # Loglik of the ith curve (logProb(Yi))

          gamma_ikjr[(((i - 1) * paramMixHMMR$fData$m + 1):(i * paramMixHMMR$fData$m)), , k] <<- t(gamma_ik) # [n*m R K] : "segments" post prob for each cluster k

          exp_num_trans_ck[, , i] <- apply(xi_ik, MARGIN = c(1, 2), sum) # [R R n]
          exp_num_trans_from_l_ck[, i] <- gamma_ik[, 1] # [R x n]

        }

        exp_num_trans_from_l[, , k] <<- exp_num_trans_from_l_ck # [R n K]
        exp_num_trans[, , , k] <<- exp_num_trans_ck # [R R n K]

        # For computing the global loglik
        # w_k_fyi[, k] <- paramMixHMMR$w_k[k] * exp(Li)#[nx1]
        log_w_k_fyi[, k] <<- log(paramMixHMMR$w_k[k]) + Li
      }

      log_w_k_fyi <<- pmin(log_w_k_fyi, log(.Machine$double.xmax))
      log_w_k_fyi <<- pmax(log_w_k_fyi, log(.Machine$double.xmin))

      tau_ik <<- exp(log_w_k_fyi) / (apply(exp(log_w_k_fyi), 1, sum) %*% matrix(1, 1, paramMixHMMR$K)) # Cluster post prob

      # Log-likelihood for the n curves
      loglik <<- sum(log(apply(exp(log_w_k_fyi), 1, sum)))

    }
  )
)
