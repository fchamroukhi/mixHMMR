myKmeans = function(X, K, nbr_runs, nbr_iter_max, verbose = FALSE) {
  #   function res = myKmeans(X, K, nbr_runs, nbr_iter_max, verbose)
  #
  #   Algorithme des K-means
  #
  #
  #
  #
  #
  #
  # distance euclidienne
  #
  ###########################################################################

  if (missing(nbr_iter_max)) {
    nbr_iter_max = 300
  }
  if (missing(nbr_runs)) {
    nbr_runs = 20
  }

  solution = list()
  solution$stored_err = numeric()

  n = nrow(X)
  p = ncol(X)
  # if one class
  global_mean = apply(X, 2, mean)

  if (K == 1) {
    dmin = rowSums((X - rep(1, n) %*% t(global_mean)) ^ 2)
    solution$muk = global_mean
    klas = rep(1, n)
    solution$klas = klas
    solution$err = sum(dmin)
    solution$Zik = rep(1, n)

    return(res)
  }

  nb_run = 0
  best_solution = list()
  best_solution$err = Inf

  while (nb_run < nbr_runs) {
    nb_run = nb_run + 1
    if (nbr_runs > 0 && verbose) {
      cat(sprintf("Kmeans : run n?? %1i", nb_run), "\n")
    }

    iter = 0
    converged = FALSE
    previous_err = -Inf
    Zik = matrix(0, nrow = n, ncol = K)

    ## 1. Initialization of the centers
    rnd_indx = sample(n)
    centres = X[rnd_indx[1:K], ]

    while ((iter < nbr_iter_max) && !(converged)) {
      iter = iter + 1
      old_centres = centres

      # The Euclidean distances
      dist_eucld = matrix(0, nrow = n, ncol = K)
      for (k in 1:K) {
        muk = centres[k,]
        dist_eucld[, k] = rowSums((X - rep(1, n) %*% t(muk)) ^ 2)
      }
      ## classification step

      dmin = apply(dist_eucld, 1, min)
      klas = apply(dist_eucld, 1, which.min)

      Zik = ((klas %*% t(rep(1, K))) == rep(1, n) %*% t(1:K)) * 1
      ## relocation step
      for (k in 1:K) {
        ind_ck = which(klas == k)
        # if empty classes
        if (length(ind_ck) == 0) {
          centres[k, ] = old_centres[k, ]
        } else  {
          # update the centres

          if (length(ind_ck) == 1) {

            centres[k,] = X[ind_ck, ]

          } else {

            centres[k,] = apply(X[ind_ck, ], 2, mean)
          }
        }
      }

      # test of convergence
      current_err = sum(rowSums(Zik * dist_eucld)) # the distorsion measure

      criteria = (abs(current_err - previous_err)) / previous_err < 1e-6

      if (iter == 1) {
        criteria = FALSE
      }

      if (criteria) {
        converged = TRUE
      }
      previous_err = current_err
      if (verbose) {
        cat(sprintf("Kmeans : Iteration  %1.1i, Objective %2.2f", iter, current_err), "\n")
      }
      solution$stored_err = cbind(solution$stored_err, current_err)
    } # one run

    solution$muk = centres
    solution$Zik = Zik
    solution$klas = klas
    solution$err = current_err

    if (current_err < best_solution$err) {
      best_solution = solution
    }
  } # end of the Kmeans runs

  return(best_solution)
}
