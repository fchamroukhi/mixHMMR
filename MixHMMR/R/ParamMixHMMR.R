source("R/enums.R")
source("R/utils.R")
source("R/myKmeans.R")
source("R/mk_stochastic.R")

ParamMixHMMR <- setRefClass(
  "ParamMixHMMR",
  fields = list(
    w_k = "matrix", # prior = "matrix",
    pi_k = "matrix", # Initial distributions
    # trans_mat = "matrix",
    A_k = "array", # Transition matrices
    beta_kr = "array", # Polynomial regression coefficient vectors
    # beta = "matrix",
    sigma_kr = "matrix",
    # Standard deviations
    mask = "matrix"
  ),
  methods = list(
    init_MixFHMMR = function(modelMixHMMR, phi, order_constraint = TRUE, init_kmeans = TRUE, try_algo = 1) {
      # # 1. Initialization of cluster weights
      w_k <<- 1 / K * matrix(1, modelMixHMMR$K, 1)

      # Initialization of the model parameters for each cluster
      if (init_kmeans) {
        max_iter_kmeans <- 400
        n_tries_kmeans <- 20
        verbose_kmeans <- 0
        solution <- myKmeans(modelMixHMMR$Y, modelMixHMMR$K, n_tries_kmeans, max_iter_kmeans, verbose_kmeans)

        for (k in 1:modelMixHMMR$K) {
          Yk <- modelMixHMMR$Y[solution$klas == k , ] #if kmeans

          init_hmm_regression(Yk, k, modelMixHMMR$R, phi, modelMixHMMR$variance_type, order_constraint, try_algo)

        }

      } else{
        ind <- sample(1:modelMixHMMR$n, modelMixHMMR$n)
        for (k in 1:modelMixHMMR$K) {
          if (k < modelMixHMMR$K) {
            Yk <- modelMixHMMR$Y[ind[((k - 1) * round(modelMixHMMR$n / modelMixHMMR$K) + 1):(k * round(modelMixHMMR$n / modelMixHMMR$K))], ]
          }
          else{
            Yk <- modelMixHMMR$Y[ind[((k - 1) * round(modelMixHMMR$n / modelMixHMMR$K) + 1):modelMixHMMR$n], ]
          }

          init_hmm_regression(Yk, k, modelMixHMMR$R, phi, modelMixHMMR$variance_type, order_constraint, try_algo)

        }
      }
    },


    init_hmm_regression = function(Y, k, R, phi, variance_type, order_constraint = TRUE, try_algo) {
      # init_hmm_regression estime les paramètres initiaux d'un modèle de regression
      # à processus markovien cache où la loi conditionnelle des observations est une gaussienne
      #
      # Entrees :
      #
      #        data  = n sequences each sequence is of m points
      #        signaux les observations sont monodimentionnelles)
      #        K : nbre d'états (classes) cachés
      #        duree_signal :  duree du signal en secondes
      #        fs : fréquence d'échantiloonnage des signaux en Hz
      #        ordre_reg : ordre de regression polynomiale
      #
      # Sorties :
      #
      #         param : parametres initiaux du modele. structure
      #         contenant les champs: para: structrure with the fields:
      #         * le HMM initial
      #         1. initial_prob (k) = Pr(Z(1) = k) avec k=1,...,K. loi initiale de z.
      #         2. trans_mat(\ell,k) = Pr(z(i)=k | z(i-1)=\ell) : matrice des transitions
      #         *
      #         3.betak : le vecteur parametre de regression associe a la classe k.
      #         vecteur colonne de dim [(p+1)x1]
      #         4. sigmak(k) = variance de x(i) sachant z(i)=k; sigmak(j) =
      #         sigma^2_k.
      #
      # Faicel Chamroukhi, Novembre 2008
      ################################################################################

      # 1. Initialization of the HMM parameters
      if (order_constraint) {
        # # Initialisation en tenant compte de la contrainte:

        # Initialisation de la matrice des transitions
        maskM <- diag(R)#mask d'ordre 1
        for (r in 1:R - 1) {
          ind <- which(maskM[r, ] != 0)
          maskM[r, ind + 1] <- 1
        }

        # Initialisation de la loi initiale de la variable cachee
        pi_k[, k] <<- c(1, matrix(0, R - 1, 1))

        A_k[, , k] <<- normalize(maskM, 2)$M
        mask <<- maskM

      } else {

        # Initialisation de la loi initiale de la variable cachee
        pi_k[, k] <<- c(1, matrix(0, R - 1, 1)) # 1 / R * matrix(1, R, 1)
        A_k[, , k] <<- mk_stochastic(matrix(runif(R), R, R))
      }

      # 2.  Initialisation of regression coefficients and variances
      init_regression_param(Y, k, R, phi, variance_type, try_algo)
    },

    ###################################################################################
    init_regression_param = function(Y, k, R, phi, variance_type, try_algo) {

      n <- nrow(Y)
      m <- ncol(Y)

      if (variance_type == variance_types$homoskedastic) {
        s <- 0
      }

      if (try_algo == 1) {

        ##############################
        #decoupage de l'echantillon (signal) en K segments
        zi <- round(m / R) - 1
        for (r in 1:R) {
          i <- (r - 1) * zi + 1
          j <- r * zi

          Yij <- Y[, i:j]
          Yij <- matrix(t(Yij), ncol = 1, byrow = T)

          phi_ij = phi[i:j, ]
          Phi_ij = repmat(phi_ij, n, 1)

          bk = solve(t(Phi_ij) %*% Phi_ij) %*% t(Phi_ij) %*% Yij
          beta_kr[, r, k] <<- bk

          if (variance_type == variance_types$homoskedastic) {
            s <- s + sum((Yij - Phi_ij %*% bk) ^ 2)
            sigma_kr[, k] <<- s / (n * m)
          } else {
            mk <- j - i + 1
            z = Yij - Phi_ij %*% bk
            sk = t(z) %*% z / (n * mk)
            sigma_kr[r, k] <<- sk
          }
        }

      } else {# initialisation aléatoire
        Lmin <- round(m / (R + 1)) #nbr pts min dans un segments
        tr_init <- matrix(0, 1, R + 1)
        tr_init[1] <- 0
        R_1 <- R
        for (r in 2:R) {
          R_1 <- R_1 - 1
          temp <- seq(tr_init[r - 1] + Lmin, m - R_1 * Lmin)
          ind <- sample(1:length(temp), length(temp))
          tr_init[r] <- temp[ind[1]]
        }

        tr_init[R + 1] <- m
        for (r in 1:R) {
          i <- tr_init[r] + 1
          j <- tr_init[r + 1]
          Yij <- Y[, i:j]
          Yij <- matrix(t(Yij), ncol = 1, byrow = T)

          phi_ij = phi[i:j, ]
          Phi_ij = repmat(phi_ij, n, 1)

          bk = solve(t(Phi_ij) %*% Phi_ij) %*% t(Phi_ij) %*% Yij
          beta_kr[, r, k] <<- bk

          if (variance_type == variance_types$homoskedastic) {
            s <- s + sum((Yij - Phi_ij %*% bk) ^ 2)
            sigma_kr[, k] <<- s / (n * m)
          } else {
            mk <- j - i + 1
            z = Yij - Phi_ij %*% bk
            sk = t(z) %*% z / (n * mk)
            sigma_kr[r, k] <<- sk
          }
        }
      }
    },


    MStep = function(modelMixHMMR, statMixHMM, phi, order_constraint = TRUE) {

      # Maximization of Q1 w.r.t w_k

      w_k <<- matrix(apply(statMixHMM$tau_ik, 2, sum)) / modelMixHMMR$n

      exp_num_trans_k <- array(0, dim = c(modelMixHMMR$R, modelMixHMMR$R, modelMixHMMR$n))

      for (k in 1:modelMixHMMR$K) {

        if (modelMixHMMR$variance_type == variance_types$homoskedastic) {
          s <- 0
        }

        weights_cluster_k <- statMixHMM$tau_ik[, k]

        # Maximization of Q2 w.r.t \pi^g
        exp_num_trans_k_from_l <- (matrix(1, modelMixHMMR$R, 1) %*% t(weights_cluster_k)) * statMixHMM$exp_num_trans_from_l[, , k] # [R x n]

        pi_k[, k] <<- (1 / sum(statMixHMM$tau_ik[, k])) * apply(exp_num_trans_k_from_l, 1, sum) # sum over i

        # Maximization of Q3 w.r.t A^g

        for (r in 1:modelMixHMMR$R) {

          if (modelMixHMMR$n == 1) {
            exp_num_trans_k[r, , ] <- t(matrix(1, modelMixHMMR$R, 1) %*% weights_cluster_k) * drop(statMixHMM$exp_num_trans[r, , , k])
          } else{
            exp_num_trans_k[r, , ] <- (matrix(1, modelMixHMMR$R, 1) %*% t(weights_cluster_k)) * drop(statMixHMM$exp_num_trans[r, , , k])
          }
        }

        if (modelMixHMMR$n == 1) {
          temp <- exp_num_trans_k
        } else{
          temp <- apply(exp_num_trans_k, MARGIN = c(1, 2), sum) # sum over i
        }


        A_k[, , k] <<- mk_stochastic(temp)

        # if HMM with order constraints
        if (order_constraint) {
          A_k[, , k] <<- mk_stochastic(mask * A_k[, , k])
        }

        # Maximisation de Q4 par rapport aux betak et sigmak

        Ng <- apply(statMixHMM$tau_ik, 2, sum) # nbr of individuals within the cluster k ,k=1...K estimated at iteration q
        ng <- Ng # cardinal nbr of the cluster k
        # each sequence i (m observations) is first weighted by the cluster weights
        weights_cluster_k <- matrix(t(statMixHMM$tau_ik[, k]), nrow = modelMixHMMR$m, ncol = ncol(t(statMixHMM$tau_ik)), byrow = T)
        weights_cluster_k <- matrix(as.vector(weights_cluster_k), length(as.vector(weights_cluster_k)), 1)

        # secondly, the m observations of each sequence are weighted by the
        # weights of each segment k (post prob of the segments for each
        # cluster g)
        gamma_ijk <- statMixHMM$gamma_ikjr[, , k] # [n*m R]

        nm_kr <- apply(gamma_ijk, 2, sum) # cardinal nbr of the segments r,r=1,...,R within each cluster k, at iteration q

        for (r in 1:modelMixHMMR$R) {
          nmkr <- nm_kr[r] # cardinal nbr of segment k for the cluster k

          weights_seg_k <- matrix(gamma_ijk[, r])

          Xkr <- (sqrt(weights_cluster_k * weights_seg_k) %*% matrix(1, 1, modelMixHMMR$p + 1)) * repmat(phi, modelMixHMMR$n, 1) # [n*m x (p+1)]
          Ykr <- (sqrt(weights_cluster_k * weights_seg_k)) * modelMixHMMR$vecY # [n*m x 1]

          # Weighted least squares: maximization w.r.t beta_kr
          beta_kr[, r, k] <<- solve(t(Xkr) %*% Xkr) %*% t(Xkr) %*% Ykr # Maximisation par rapport aux betakr
          # W_kr = diag(weights_cluster_k.*weights_seg_k);
          # beta_kr(:,k) = inv(Phi'*W_kr*Phi)*Phi'*W_kr*X;

          # Maximization w.r.t sigmak :

          z <- sqrt(weights_cluster_k * weights_seg_k) * (modelMixHMMR$vecY - repmat(phi, modelMixHMMR$n, 1) %*% beta_kr[, r, k])

          if (modelMixHMMR$variance_type == variance_types$homoskedastic) {
            s <- s + (t(z) %*% z)
            ngm <- sum(sum(weights_cluster_k %*% matrix(1, 1, modelMixHMMR$R)) * gamma_ijk)

            # sigma_k <- s/ngm
            sigma_kr[k] <<- s / ngm
          }
          else{
            ngmk <- sum(weights_cluster_k * weights_seg_k)

            sigma_kr[r, k] <<- (t(z) %*% z) / ngmk
          }
        }
      }
    }
  )
)

ParamMixHMMR <- function(modelMixHMMR) {
  w_k <- matrix(NA, nrow = modelMixHMMR$K)
  # prior <- matrix(NA, ncol = modelMixHMMR$K - 1)
  pi_k <- matrix(NA, nrow = modelMixHMMR$R, ncol = modelMixHMMR$K)
  # trans_mat <- matrix(NA, modelMixHMMR$K, modelMixHMMR$K)
  A_k <- array(NA, dim = c(modelMixHMMR$R, modelMixHMMR$R, modelMixHMMR$K))
  # beta <- matrix(NA, modelMixHMMR$p + 1, modelMixHMMR$K)
  beta_kr <- array(NA, dim = c(modelMixHMMR$p + 1, modelMixHMMR$R, modelMixHMMR$K))
  if (modelMixHMMR$variance_type == variance_types$homoskedastic) {
    sigma_kr <- matrix(NA, ncol = modelMixHMMR$K)
  }
  else{
    sigma_kr <- matrix(NA, nrow = modelMixHMMR$R, ncol = modelMixHMMR$K)
  }
  mask <- matrix(NA, modelMixHMMR$R, modelMixHMMR$R)
  new(
    "ParamMixHMMR",
    w_k = w_k,
    pi_k = pi_k,
    A_k = A_k,
    beta_kr = beta_kr,
    sigma_kr = sigma_kr,
    mask = mask
  )
}
