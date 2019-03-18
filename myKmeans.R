myKmeans = function(X, K, nbr_runs, nbr_iter_max, verbose)
{
  #   function res = myKmeans(X, K, nbr_runs, nbr_iter_max, verbose)
  #
  #   Algorithme des K-means
  #
  #
  #
  # Faicel CHAMROUKHI Septembre 2008 (mise a jour)
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
  
  n = nrow(X)
  p = ncol(X)
  
  # if one class
  global_mean = apply(X, 2, mean)

  solution = list()
  
  if (K==1)
  {
    dmin = rowSums((X - rep(1, n) %*% t(global_mean)) ^ 2)
    solution$muk = global_mean
    klas = rep(1, n)
    solution$klas = klas
    solution$err = sum(dmin)
    solution$Zik = rep(1, n)
    
    return(solution)
    
  } else {
    nbr_run = 0
    best_solution = list()
    best_solution$err = Inf
    
    while (nbr_run<nbr_runs) {
      nbr_run=nbr_run + 1
      
      if (nbr_runs>1 && verbose) {
        cat("Kmeans run n° : ",nbr_runs,"\n")  
      } 
      iter=0
      converged=0
      previous_err = -Inf
      Zik = matrix(0, nrow = n, ncol = K) #partition
      
      ## 1. Initialization of the centres
      rnd_indx = sample(n)
      centres = X[rnd_indx[1:K], ]
      
      while ((iter <= nbr_iter_max) && !(converged)) {
        iter = iter+1
        old_centres = centres
        
        # The Euclidean distances
        eucld_dist = matrix(0, nrow = n, ncol = K)
          
        for (k in 1:K) {
          muk = centres[k,]
          eucld_dist[, k] = rowSums((X - rep(1, n) %*% t(muk)) ^ 2)
        }
        ## classification step
        
        dmin = apply(eucld_dist, 1, min)
        klas = apply(eucld_dist, 1, which.min)
        Zik = ((klas %*% t(rep(1, K))) == rep(1, n) %*% t(1:K)) * 1
        
        for (k in 1:K) {
          ind_ck = which(klas == k)
          #if empty classes
          if (length(ind_ck) == 0) {
            centres[k,]= old_centres[k,]
          } else if (length(ind_ck) == 1) {
            centres[k,]= X[ind_ck, ]
          } else{
            # update the centres
            centres[k,]= apply(X[ind_ck, ], 2, mean)
          }
        }
        
        # test of convergence
        current_err = sum(rowSums(Zik * eucld_dist)) # the distorsion measure
        
        test1= (abs(current_err-previous_err))/previous_err <1e-6
        if (is.na(test1)) {
          test1 = FALSE
        }
        
        if (test1) {
          converged=TRUE
        }
        previous_err = current_err
        if (verbose) {
          cat("Kmeans : Iteration ",iter ," Objective: ", current_err,"\n")  
        }
        solution$stored_err[iter] = current_err
      } #one run
      
      solution$muk = centres
      solution$Zik = Zik
      solution$klas = klas
      solution$err = current_err
      
      if (current_err < best_solution$err) {
        best_solution = solution
      }
    }
    solution = best_solution
    return(solution)
  }
}