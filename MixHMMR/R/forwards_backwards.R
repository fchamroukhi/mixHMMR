forwards_backwards <- function(prior, transmat, f_tk) {
  #[tau_tk, xi_ikl, alpha, beta, loglik] = forwards_backwards(prior, transmat, fik, filter_only)
  # forwards_backwards : calculates the E-step of the EM algorithm for an HMM
  # (Gaussian HMM)

  # Inputs :
  #
  #         prior(k) = Pr(z_1 = k)
  #         transmat(\ell,k) = Pr(z_t=k | z_{t-1} = \ell)
  #         f_tk(t,k) = Pr(y_t | z_y=k;\theta) %gaussian
  #
  # Outputs:
  #
  #        tau_tk(t,k) = Pr(z_t=k | X): post probs (smoothing probs)
  #        xi_tk\elll(t,k,\ell)  = Pr(z_t=k, z_{t-1}=\ell | Y) t =2,..,n
  #        with Y = (y_1,...,y_n);
  #        alpha_tk: [nxK], forwards probs: Pr(y1...yt,zt=k)
  #        beta_tk: [nxK], backwards probs: Pr(yt+1...yn|zt=k)
  #
  #
  #
  # Faicel Chamroukhi
  ##############################################################################
  N <- ncol(f_tk)
  K <- length(prior)
  if (nargs() < 6) {
    filter_only <- 0
  }

  if (length(prior) != 1) {
    prior <- t(prior)
  }

  scale <- matrix(1, nrow = 1, ncol = N) #pour que loglik = sum(log(scale)) part de zero

  prior <- prior
  tau_tk <- matrix(0, nrow = K, ncol = N)
  xi_tk <- array(0, dim = c(K, K, N - 1))
  xi_summed <- matrix(0, nrow = K, ncol = K)
  alpha_tk <- matrix(0, nrow = K, ncol = N)
  beta_tk <- matrix(0, nrow = K, ncol = N)

  ## forwards: calculation of alpha_tk
  t <- 1
  alpha_tk[, t] <- prior * f_tk[, t]
  alpha_tk[, t] <- normalize(alpha_tk[, t])[[1]]
  scale[t] <- normalize(alpha_tk[, t])[[2]]

  for (t in 2:N) {
    alpha_tk[, t] <- normalize((t(transmat) %*% alpha_tk[, t - 1]) * f_tk[, t])[[1]]
    scale[t] <- normalize((t(transmat) %*% alpha_tk[, t - 1]) * f_tk[, t])[[2]]
    #filtered_prob (t-1,:,:)= normalize((alpha(:,t-1) * fik(:,t)') .*transmat)
  }

  ##loglikehood (with the scaling technique) (see Rabiner's paper/book)
  loglik <- sum(log(scale))

  if (filter_only) {
    beta_tk <- c()
    xi_tk <- alpha_tk
  }
  ## backwards: calculation of beta_tk, tau_tk (and xi_tkl)
  #t<-T
  beta_tk[, N] <- matrix(1, 1, K)
  tau_tk[, N] <- normalize(alpha_tk[, N] * beta_tk[, N])[[1]]
  for (t in (N - 1):1) {
    beta_tk[, t] <- normalize(transmat %*% (beta_tk[, t + 1] * f_tk[, t + 1]))[[1]]
    # transmat * t(beta[t+1,] %*% fik[t+1,]) /scale[t]
    tau_tk[, t] <- normalize(alpha_tk[, t] * beta_tk[, t])[[1]]
    xi_tk[, , t] <- normalize(transmat * (alpha_tk[, t] %*% t(beta_tk[, t + 1] * f_tk[, t + 1])))[[1]]
    xi_summed <- xi_summed + xi_tk[, , t]
  }
  fb <-
    list(
      tau_tk = tau_tk,
      xi_tk = xi_tk,
      alpha_tk = alpha_tk,
      beta_tk = beta_tk,
      loglik = loglik
    )
  return(fb)
}
