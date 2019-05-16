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
    # tau_tk = "matrix", # tau_tk: smoothing probs: [nxK], tau_tk(t,k) = Pr(z_i=k | y1...yn)
    # alpha_tk = "matrix", # alpha_tk: [nxK], forwards probs: Pr(y1...yt,zt=k)
    # beta_tk = "matrix", # beta_tk: [nxK], backwards probs: Pr(yt+1...yn|zt=k)
    # xi_tkl = "array", # xi_tkl: [(n-1)xKxK], joint post probs : xi_tk\elll(t,k,\ell)  = Pr(z_t=k, z_{t-1}=\ell | Y) t =2,..,n
    # f_tk = "matrix", # f_tk: [nxK] f(yt|zt=k)
    # log_f_tk = "matrix", # log_f_tk: [nxK] log(f(yt|zt=k))
    # loglik = "numeric", # loglik: log-likelihood at convergence
    # stored_loglik = "list", # stored_loglik: stored log-likelihood values during EM
    # cputime = "numeric", # cputime: for the best run
    # klas = "matrix", # klas: [nx1 double]
    # z_ik = "matrix", # z_ik: [nxK]
    # state_probs = "matrix", # state_probs: [nxK]
    # BIC = "numeric", # BIC
    # AIC = "numeric", # AIC
    # regressors = "matrix", # regressors: [nxK]
    # predict_prob = "matrix", # predict_prob: [nxK]: Pr(zt=k|y1...y_{t-1})
    # predicted = "matrix", # predicted: [nx1]
    # filter_prob = "matrix", # filter_prob: [nxK]: Pr(zt=k|y1...y_t)
    # filtered = "matrix", # filtered: [nx1]
    # smoothed_regressors = "matrix", # smoothed_regressors: [nxK]
    # smoothed = "matrix" # smoothed: [nx1]
    # #           X: [nx(p+1)] regression design matrix
    # #           nu: model complexity
    # #           parameter_vector
  ),
  methods = list(
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
    #######
    # compute the final solution stats
    #######
    computeStats = function(modelMixHMMR, paramMixHMMR, phi, cputime_total) {

      cputime <<- mean(cputime_total)

      for (k in 1:modelMixHMMR$K) {

        betakr <- paramMixHMMR$beta_kr[, , k]
        weighted_segments <- apply(gamma_ikjr[, , k] * (repmat(phi, modelMixHMMR$n, 1) %*% betakr), 1, sum)
        dim(weighted_segments) <- c(modelMixHMMR$m, modelMixHMMR$n)
        weighted_clusters <- (matrix(1, modelMixHMMR$m, 1) %*% t(tau_ik[, k])) * weighted_segments
        smoothed[, k] <<- apply(weighted_clusters, 1, sum) / sum(tau_ik[, k])
      }

      # BIC AIC et ICL*
      BIC <<- loglik - (modelMixHMMR$nu * log(modelMixHMMR$n) / 2)
      AIC <<- loglik - modelMixHMMR$nu
      # ICL*
      # Compute the comp-log-lik
      cik_log_w_k_fyi <- z_ik * log_w_k_fyi
      comp_loglik <- sum(cik_log_w_k_fyi)
      ICL1 <<- comp_loglik - modelMixHMMR$nu * log(modelMixHMMR$n) / 2 #n*m/2!

    },
    #######
    # EStep
    #######
    EStep = function(modelMixHMMR, paramMixHMMR, phi) {

      exp_num_trans_ck  <- array(0, dim = c(modelMixHMMR$R, modelMixHMMR$R, modelMixHMMR$n))
      exp_num_trans_from_l_ck <- matrix(0, modelMixHMMR$R, modelMixHMMR$n)

      # w_k_fyi <- matrix(0, modelMixHMMR$n, modelMixHMMR$K)

      for (k in 1:modelMixHMMR$K) {
        # run a hmm for each sequence
        log_fkr_yij <- matrix(0, modelMixHMMR$R, modelMixHMMR$m)
        fkr_yij <- matrix(0, modelMixHMMR$R, modelMixHMMR$m)
        #
        Li <- matrix(0, modelMixHMMR$n, 1)# to store the loglik for each example (curve)
        #
        beta_kr <- paramMixHMMR$beta_kr[, , k]
        num_log_post_prob <- matrix(0, modelMixHMMR$n, modelMixHMMR$K)

        for (i in 1:modelMixHMMR$n) {
          y_i <- modelMixHMMR$Y[i, ]

          for (r in 1:modelMixHMMR$R) {
            betakr <- beta_kr[, r]

            if (modelMixHMMR$variance_type == variance_types$homoskedastic) {
              sigma_kr <- paramMixHMMR$sigma_kr[, k]
              sk <- sigma_kr
            }
            else{
              sigma_kr <- paramMixHMMR$sigma_kr[, k]
              sk <- sigma_kr[r]
            }
            z <- ((y_i - t(phi %*% betakr)) ^ 2) / sk

            log_fkr_yij[r, ] <- -0.5 * matrix(1, 1, modelMixHMMR$m) * (log(2 * pi) + log(sk)) - 0.5 * z # log pdf yij | c_i = k et z_i = r
            # fkr_yij[r, ] <- dnorm(y_i, t(phi %*% beta_kr), sqrt(sk))

          }

          log_fkr_yij <- pmin(log_fkr_yij, log(.Machine$double.xmax))
          log_fkr_yij <- pmax(log_fkr_yij, log(.Machine$double.xmin))
          fkr_yij <- exp(log_fkr_yij)


          # forwards backwards ( calcul de logProb(Yi)...)

          fb <- forwards_backwards(paramMixHMMR$pi_k[, k], paramMixHMMR$A_k[, , k], fkr_yij)

          gamma_ik <- fb$tau_tk
          xi_ik <- fb$xi_tk
          fwd_ik <- fb$alpha_tk
          backw_ik <- fb$beta_tk
          loglik_i <- fb$loglik

          #
          Li[i] <- loglik_i # loglik of the ith curve (logProb(Yi))

          #
          gamma_ikjr[(((i - 1) * modelMixHMMR$m + 1):(i * modelMixHMMR$m)), , k] <<- t(gamma_ik)#[n*m R K] : "segments" post prob for each cluster k
          #
          exp_num_trans_ck[, , i] <- apply(xi_ik, MARGIN = c(1, 2), sum) # [R R n]
          exp_num_trans_from_l_ck[, i] <- gamma_ik[, 1] # [R x n]
          #
        }

        exp_num_trans_from_l[, , k] <<- exp_num_trans_from_l_ck # [R n K]
        exp_num_trans[, , , k] <<- exp_num_trans_ck # [R R n K]

        # for computing the global loglik
        # w_k_fyi[, k] <- paramMixHMMR$w_k[k] * exp(Li)#[nx1]
        log_w_k_fyi[, k] <<- log(paramMixHMMR$w_k[k]) + Li
      }

      log_w_k_fyi <<- pmin(log_w_k_fyi, log(.Machine$double.xmax))
      log_w_k_fyi <<- pmax(log_w_k_fyi, log(.Machine$double.xmin))

      tau_ik <<- exp(log_w_k_fyi) / (apply(exp(log_w_k_fyi), 1, sum) %*% matrix(1, 1, modelMixHMMR$K)) # cluster post prob

      # log-likelihood for the n curves
      loglik <<- sum(log(apply(exp(log_w_k_fyi), 1, sum)))

    }
  )
)


StatMixHMMR <- function(modelMixHMMR) {
  tau_ik <- matrix(NA, modelMixHMMR$n, modelMixHMMR$K)
  gamma_ikjr <- array(NA, dim = c(modelMixHMMR$n * modelMixHMMR$m, modelMixHMMR$R, modelMixHMMR$K))
  log_w_k_fyi <- matrix(NA, modelMixHMMR$n, modelMixHMMR$K)
  exp_num_trans <- array(NA, dim = c(modelMixHMMR$R, modelMixHMMR$R, modelMixHMMR$n, modelMixHMMR$K))
  exp_num_trans_from_l <- array(NA, dim = c(modelMixHMMR$R, modelMixHMMR$n, modelMixHMMR$K))
  loglik <- -Inf
  stored_loglik <- list()
  cputime <- Inf
  klas <- matrix(NA, modelMixHMMR$n, 1) # klas: [nx1 double]
  z_ik <- matrix(NA, modelMixHMMR$n, modelMixHMMR$K) # z_ik: [nxK]
  smoothed <- matrix(NA, modelMixHMMR$m, modelMixHMMR$K)
  mean_curves <- array(NA, dim = c(modelMixHMMR$m, modelMixHMMR$R, modelMixHMMR$K))
  BIC <- -Inf
  AIC <- -Inf
  ICL1 <- -Inf
  # tau_tk <- matrix(NA, modelMixHMMR$m, modelMixHMMR$K) # tau_tk: smoothing probs: [nxK], tau_tk(t,k) = Pr(z_i=k | y1...yn)
  # alpha_tk <- matrix(NA, modelMixHMMR$m, ncol = modelMixHMMR$K) # alpha_tk: [nxK], forwards probs: Pr(y1...yt,zt=k)
  # beta_tk <- matrix(NA, modelMixHMMR$m, modelMixHMMR$K) # beta_tk: [nxK], backwards probs: Pr(yt+1...yn|zt=k)
  # xi_tkl <- array(NA, c(modelMixHMMR$m - 1, modelMixHMMR$K, modelMixHMMR$K)) # xi_tkl: [(n-1)xKxK], joint post probs : xi_tk\elll(t,k,\ell)  = Pr(z_t=k, z_{t-1}=\ell | Y) t =2,..,n
  # f_tk <- matrix(NA, modelMixHMMR$m, modelMixHMMR$K) # f_tk: [nxK] f(yt|zt=k)
  # log_f_tk <- matrix(NA, modelMixHMMR$m, modelMixHMMR$K) # log_f_tk: [nxK] log(f(yt|zt=k))
  # loglik <- -Inf # loglik: log-likelihood at convergence
  # stored_loglik <- list() # stored_loglik: stored log-likelihood values during EM
  # cputime <- Inf # cputime: for the best run
  # klas <- matrix(NA, modelMixHMMR$m, 1) # klas: [nx1 double]
  # z_ik <- matrix(NA, modelMixHMMR$m, modelMixHMMR$K) # z_ik: [nxK]
  # state_probs <- matrix(NA, modelMixHMMR$m, modelMixHMMR$K) # state_probs: [nxK]
  # BIC <- -Inf # BIC
  # AIC <- -Inf # AIC
  # regressors <- matrix(NA, modelMixHMMR$m, modelMixHMMR$K) # regressors: [nxK]
  # predict_prob <- matrix(NA, modelMixHMMR$m, modelMixHMMR$K) # predict_prob: [nxK]: Pr(zt=k|y1...y_{t-1})
  # predicted <- matrix(NA, modelMixHMMR$m, 1) # predicted: [nx1]
  # filter_prob <- matrix(NA, modelMixHMMR$m, modelMixHMMR$K) # filter_prob: [nxK]: Pr(zt=k|y1...y_t)
  # filtered <- matrix(NA, modelMixHMMR$m, 1) # filtered: [nx1]
  # smoothed_regressors <- matrix(NA, modelMixHMMR$m, modelMixHMMR$K) # smoothed_regressors: [nxK]
  # smoothed <- matrix(NA, modelMixHMMR$m, 1) # smoothed: [nx1]

  new(
    "StatMixHMMR",
    tau_ik = tau_ik,
    gamma_ikjr = gamma_ikjr,
    log_w_k_fyi = log_w_k_fyi,
    exp_num_trans = exp_num_trans,
    exp_num_trans_from_l = exp_num_trans_from_l,
    loglik = loglik,
    stored_loglik = stored_loglik,
    cputime = cputime,
    klas = klas,
    z_ik = z_ik,
    smoothed = smoothed,
    mean_curves = mean_curves,
    BIC = BIC,
    AIC = AIC,
    ICL1 = ICL1
  )
}
