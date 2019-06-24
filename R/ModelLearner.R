#' @export
emMixHMMR <- function(X, Y, K, R, p = 3, variance_type = c("heteroskedastic", "homoskedastic"), order_constraint = TRUE, n_tries = 1, max_iter = 1000, init_kmeans = TRUE, threshold = 1e-6, verbose = TRUE) {
    #   MixFHMMR =  seq_clust_MixFHMMR(data, K, R, p,fs, variance_type,...
    #               order_constraint, total_EM_tries, max_iter_EM, init_kmeans, threshold, verbose)
    # Learn a mixture of Hidden Markov Moedel Regression for curve clustering by EM
    #
    #
    # Inputs :
    #
    #          1. data :  n curves each curve is composed of m points : dim(Y)=[n m]
    #                * Each curve is observed during the interval [0,T]=[t_1,...,t_m]
    #                * t{j}-t_{j-1} = 1/fs (fs: sampling period)
    #          2. K: number of clusters
    #          3. R: Number of polynomial regression components (regimes)
    #          4. p: degree of the polynomials
    # Options:
    #          1. order_constraint: set to one if ordered segments (by default 0)
    #          2. variance_type of the poynomial models for each cluster (free or
    #          common, by defalut free)
    #          3. init_kmeans: initialize the curve partition by Kmeans
    #          4. total_EM_tries :  (the solution providing the highest log-lik is chosen
    #          5. max_iter_EM
    #          6. threshold: by defalut 1e-6
    #          7. verbose : set to 1 for printing the "complete-log-lik"  values during
    #          the EM iterations (by default verbose_EM = 0)
    #
    # Outputs :
    #
    #          MixFHMMR : structure containing the following fields:
    #
    #          1. param : a structure containing the model parameters
    #                       ({Wk},{alpha_k}, {beta_kr},{sigma_kr}) for k=1,...,K and k=1...R.
    #              1.1 Wk = (Wk1,...,w_kR-1) parameters of the logistic process:
    #                  matrix of dimension [(q+1)x(R-1)] with q the order of logistic regression.
    #              1.2 beta_k = (beta_k1,...,beta_kR) polynomial regression coefficient vectors: matrix of
    #                  dimension [(p+1)xR] p being the polynomial  degree.
    #              1.3 sigma_k = (sigma_k1,...,sigma_kR) : the variances for the R regmies. vector of dimension [Rx1]
    #              1.4 pi_jkr :logistic proportions for cluster g
    #
    #          2. paramter_vector: parameter vector of the model: Psi=({Wg},{alpha_k},{beta_kr},{sigma_kr})
    #                  column vector of dim [nu x 1] with nu = nbr of free parametres
    #          3. h_ik = prob(curve|cluster_k) : post prob (fuzzy segmentation matrix of dim [nxK])
    #          4. c_ik : Hard partition obtained by the AP rule :  c_{ik} = 1
    #                    if and only c_i = arg max_k h_ik (k=1,...,K)
    #          5. klas : column vector of cluster labels
    #          6. tau_ijkr prob(y_{ij}|kth_segment,cluster_k), fuzzy
    #          segmentation for the cluster g. matrix of dimension
    #          [nmxR] for each g  (g=1,...,K).
    #          7. Ex_k: curve expectation: sum of the polynomial components beta_kr ri weighted by
    #             the logitic probabilities pij_kr: Ex_k(j) = sum_{k=1}^R pi_jkr beta_kr rj, j=1,...,m. Ex_k
    #              is a column vector of dimension m for each g.
    #          8. loglik : at convergence of the EM algo
    #          9. stored_com-loglik : vector of stored valued of the
    #          comp-log-lik at each EM teration
    #
    #          10. BIC value = loglik - nu*log(nm)/2.
    #          11. ICL value = comp-loglik_star - nu*log(nm)/2.
    #          12. AIC value = loglik - nu.
    #          13. log_alphag_fg_xij
    #          14. polynomials
    #          15. weighted_polynomials
    #
    #
    #############################################################################################

    fData <- FData$new(X = X, Y = Y)

    try_EM <- 0
    best_loglik <- -Inf
    cputime_total <- c()

    while (try_EM < n_tries) {
      try_EM <- try_EM + 1
      if (n_tries > 1) {
        if (verbose) {
          cat(paste0("EM try number: ", try_EM + 1, "\n\n"))
        }
      }
      start_time <- Sys.time()

      # Initialization
      variance_type <- match.arg(variance_type)
      param <- ParamMixHMMR$new(fData = fData, K = K, R = R, p = p, variance_type = variance_type)
      param$initMixFHMMR(order_constraint, init_kmeans, try_EM)

      iter <- 0
      converged <- FALSE
      prev_loglik <- -Inf

      stat <- StatMixHMMR$new(paramMixHMMR = param)

      # EM
      while ((iter <= max_iter) & !converged) {

        # E-Step
        stat$EStep(param)

        # M-Step
        param$MStep(stat, order_constraint)

        iter <- iter + 1

        if (verbose) {
          cat(paste0("EM: Iteration : ", iter, " || log-likelihood : "  , stat$loglik, "\n"))
        }

        if (prev_loglik - stat$loglik > 1e-4) {
          warning(paste0("EM log-likelihood is decreasing from ", prev_loglik, "to ", stat$loglik, " !"))
        }

        converged <-
          (abs((stat$loglik - prev_loglik) / prev_loglik) < threshold)
        if (is.na(converged)) {
          converged <- FALSE
        } # Basically for the first iteration when prev_loglik is Inf

        prev_loglik <- stat$loglik
        stat$stored_loglik[iter] <- stat$loglik

      } # End of EM loop

      cputime_total <- cbind(cputime_total, Sys.time() - start_time)

      if (stat$loglik > best_loglik) {
        statSolution <- stat$copy()
        paramSolution <- param$copy()

        best_loglik <- stat$loglik
      }

      if (try_EM > 1) {
        if (verbose) {
          cat(paste0("Log-likelihood at convergence:", stat$loglik))
        }
      }

    }

    if (try_EM > 1) {
      if (verbose) {
        cat(paste0("Best value of the log-likelihood: ", statSolution$loglik, "\n"))
      }
    }

    # Finding the curve partition by using the MAP rule
    statSolution$MAP()

    # Finish computation of statSolution
    statSolution$computeStats(param, cputime_total)

    return(ModelMixHMMR(paramMixHMMR = paramSolution, statMixHMMR = statSolution))

  }
