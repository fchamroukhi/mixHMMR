#' @export
ParamMixHMMR <- setRefClass(
  "ParamMixHMMR",
  fields = list(
    fData = "FData",
    phi = "matrix",

    K = "numeric", # Number of clusters
    R = "numeric", # Number of regimes (HMM states)
    p = "numeric", # Dimension of beta (order of polynomial regression)
    variance_type = "character",
    nu = "numeric", # Degree of freedom

    alpha = "matrix", # Cluster weights
    prior = "matrix", # Initial distributions
    trans_mat = "array", # Transition matrices
    beta = "array", # Polynomial regression coefficient vectors
    sigma2 = "matrix", # Variances
    mask = "matrix"
  ),
  methods = list(
    initialize = function(fData = FData(numeric(1), matrix(1)), K = 2, R = 1, p = 3, variance_type = "heteroskedastic") {

      if (variance_type == "homoskedastic") {
        nu <<- (K - 1) + K * ((R - 1) + R * (R - 1) + R * (p + 1) + 1)
      }
      else{
        nu <<- (K - 1) + K * ((R - 1) + R * (R - 1) + R * (p + 1) + R)
      }

      fData <<- fData

      phi <<- designmatrix(x = fData$X, p = p)$XBeta

      K <<- K
      R <<- R
      p <<- p
      variance_type <<- variance_type

      alpha <<- matrix(NA, nrow = K)
      prior <<- matrix(NA, nrow = R, ncol = K)
      trans_mat <<- array(NA, dim = c(R, R, K))
      beta <<- array(NA, dim = c(p + 1, R, K))

      if (variance_type == "homoskedastic") {
        sigma2 <<- matrix(NA, ncol = K)
      } else {
        sigma2 <<- matrix(NA, nrow = R, ncol = K)
      }
      mask <<- matrix(NA, R, R)

    },

    initMixFHMMR = function(order_constraint = TRUE, init_kmeans = TRUE, try_algo = 1) {

      # 1. Initialization of cluster weights
      alpha <<- 1 / K * matrix(1, K, 1)

      # Initialization of the model parameters for each cluster
      if (init_kmeans) {
        max_iter_kmeans <- 400
        n_tries_kmeans <- 20
        verbose_kmeans <- 0
        solution <- kmeans(fData$Y, K, n_tries_kmeans, max_iter_kmeans, verbose_kmeans)

        for (k in 1:K) {
          Yk <- fData$Y[solution$klas == k ,] #if kmeans
          initHmmRegression(Yk, k, R, phi, variance_type, order_constraint, try_algo)
        }

      } else {
        ind <- sample(1:fData$n, fData$n)
        for (k in 1:K) {
          if (k < K) {
            Yk <- fData$Y[ind[((k - 1) * round(fData$n / K) + 1):(k * round(fData$n / K))],]
          } else {
            Yk <- fData$Y[ind[((k - 1) * round(fData$n / K) + 1):fData$n],]
          }

          initHmmRegression(Yk, k, R, phi, variance_type, order_constraint, try_algo)

        }
      }
    },

    initHmmRegression = function(Y, k, R, phi, variance_type, order_constraint = TRUE, try_algo) {
      # initHmmRegression estime les parametres initiaux d'un modele de regression
      # processus markovien cache ou la loi conditionnelle des observations est une gaussienne
      #
      # Entrees :
      #
      #        data  = n sequences each sequence is of m points
      #        signaux les observations sont monodimentionnelles)
      #        K : nbre d'etats (classes) caches
      #        duree_signal :  duree du signal en secondes
      #        fs : frequence d'echantiloonnage des signaux en Hz
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
      ################################################################################

      # 1. Initialization of the HMM parameters
      if (order_constraint) {

        # Initialization taking into account the constraint:

        # Initialization of the transition matrix
        maskM <- diag(R) # Mask of order 1
        for (r in 1:R - 1) {
          ind <- which(maskM[r,] != 0)
          maskM[r, ind + 1] <- 1
        }

        # Initialization of the initial distribution
        prior[, k] <<- c(1, matrix(0, R - 1, 1))

        trans_mat[, , k] <<- normalize(maskM, 2)$M
        mask <<- maskM

      } else {
        # Initialization of the initial distribution
        prior[, k] <<- c(1, matrix(0, R - 1, 1))
        trans_mat[, , k] <<- mkStochastic(matrix(runif(R), R, R))
      }

      # 2. Initialisation of regression coefficients and variances
      initRegressionParam(Y, k, R, phi, variance_type, try_algo)
    },

    initRegressionParam = function(Y, k, R, phi, variance_type, try_algo) {
      n <- nrow(Y)
      m <- ncol(Y)

      if (variance_type == "homoskedastic") {
        s <- 0
      }

      if (try_algo == 1) {
        zi <- round(m / R) - 1
        for (r in 1:R) {
          i <- (r - 1) * zi + 1
          j <- r * zi

          Yij <- Y[, i:j]
          Yij <- matrix(t(Yij), ncol = 1, byrow = T)

          phi_ij <- phi[i:j,]
          Phi_ij <- repmat(phi_ij, n, 1)

          bk <- solve(t(Phi_ij) %*% Phi_ij) %*% t(Phi_ij) %*% Yij
          beta[, r, k] <<- bk

          if (variance_type == "homoskedastic") {
            s <- s + sum((Yij - Phi_ij %*% bk) ^ 2)
            sigma2[, k] <<- s / (n * m)
          } else {
            mk <- j - i + 1
            z <- Yij - Phi_ij %*% bk
            sk <- t(z) %*% z / (n * mk)
            sigma2[r, k] <<- sk
          }
        }

      } else {
        Lmin <- round(m / (R + 1))
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

          phi_ij <- phi[i:j,]
          Phi_ij <- repmat(phi_ij, n, 1)

          bk <- solve(t(Phi_ij) %*% Phi_ij) %*% t(Phi_ij) %*% Yij
          beta[, r, k] <<- bk

          if (variance_type == "homoskedastic") {
            s <- s + sum((Yij - Phi_ij %*% bk) ^ 2)
            sigma2[, k] <<- s / (n * m)
          } else {
            mk <- j - i + 1
            z <- Yij - Phi_ij %*% bk
            sk <- t(z) %*% z / (n * mk)
            sigma2[r, k] <<- sk
          }
        }
      }
    },

    MStep = function(statMixHMM, order_constraint = TRUE) {

      # Maximization of Q1 w.r.t alpha
      alpha <<- matrix(apply(statMixHMM$tau_ik, 2, sum)) / fData$n

      exp_num_trans_k <- array(0, dim = c(R, R, fData$n))

      for (k in 1:K) {
        if (variance_type == "homoskedastic") {
          s <- 0
        }

        weights_cluster_k <- statMixHMM$tau_ik[, k]

        # Maximization of Q2 w.r.t \pi^g
        exp_num_trans_k_from_l <- (matrix(1, R, 1) %*% t(weights_cluster_k)) * statMixHMM$exp_num_trans_from_l[, , k] # [R x n]

        prior[, k] <<- (1 / sum(statMixHMM$tau_ik[, k])) * apply(exp_num_trans_k_from_l, 1, sum) # sum over i

        # Maximization of Q3 w.r.t A^g
        for (r in 1:R) {
          if (fData$n == 1) {
            exp_num_trans_k[r, ,] <- t(matrix(1, R, 1) %*% weights_cluster_k) * drop(statMixHMM$exp_num_trans[r, , , k])
          } else {
            exp_num_trans_k[r, ,] <- (matrix(1, R, 1) %*% t(weights_cluster_k)) * drop(statMixHMM$exp_num_trans[r, , , k])
          }
        }

        if (fData$n == 1) {
          temp <- exp_num_trans_k
        } else {
          temp <- apply(exp_num_trans_k, MARGIN = c(1, 2), sum) # sum over i
        }

        trans_mat[, , k] <<- mkStochastic(temp)

        # If HMM with order constraints
        if (order_constraint) {
          trans_mat[, , k] <<- mkStochastic(mask * trans_mat[, , k])
        }

        # Maximisation of Q4 w.r.t with betak et sigmak

        Ng <- apply(statMixHMM$tau_ik, 2, sum) # Nbr of individuals within the cluster k ,k=1...K estimated at iteration q
        ng <- Ng # Cardinal nbr of the cluster k
        # Each sequence i (m observations) is first weighted by the cluster weights
        weights_cluster_k <- matrix(t(statMixHMM$tau_ik[, k]), nrow = fData$m, ncol = ncol(t(statMixHMM$tau_ik)), byrow = T)
        weights_cluster_k <- matrix(as.vector(weights_cluster_k), length(as.vector(weights_cluster_k)), 1)

        # Secondly, the m observations of each sequence are weighted by the
        # weights of each segment k (post prob of the segments for each
        # cluster g)

        gamma_ijk <- statMixHMM$gamma_ikjr[, , k] # [n*m R]

        nm_kr <- apply(as.matrix(gamma_ijk), 2, sum) # Cardinal nbr of the segments r,r=1,...,R within each cluster k, at iteration q

        for (r in 1:R) {
          nmkr <- nm_kr[r] # Cardinal nbr of segment k for the cluster k

          weights_seg_k <- matrix(as.matrix(gamma_ijk)[, r])

          Xkr <- (sqrt(weights_cluster_k * weights_seg_k) %*% matrix(1, 1, p + 1)) * repmat(phi, fData$n, 1) # [n*m x (p+1)]
          Ykr <- (sqrt(weights_cluster_k * weights_seg_k)) * fData$vecY # [n*m x 1]

          # Weighted least squares: maximization w.r.t beta
          beta[, r, k] <<- solve(t(Xkr) %*% Xkr) %*% t(Xkr) %*% Ykr # Maximization w.r.t beta

          # Maximization w.r.t sigmak :
          z <- sqrt(weights_cluster_k * weights_seg_k) * (fData$vecY - repmat(phi, fData$n, 1) %*% beta[, r, k])

          if (variance_type == "homoskedastic") {
            s <- s + (t(z) %*% z)
            ngm <- sum(sum(weights_cluster_k %*% matrix(1, 1, R)) * gamma_ijk)

            sigma2[k] <<- s / ngm
          } else {
            ngmk <- sum(weights_cluster_k * weights_seg_k)

            sigma2[r, k] <<- (t(z) %*% z) / ngmk
          }
        }
      }
    }
  )
)
