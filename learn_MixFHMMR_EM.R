learn_MixFHMMR_EM = function(data, K, R, p, variance_type, ordered_states, total_EM_tries, max_iter_EM, init_kmeans, threshold, verbose)
{
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
  ################################## Faicel Chamroukhi (septembre 2009) ###########################################################

  n=nrow(Y)
  m=ncol(Y) #n  curves of m observations
  
  
  exp_num_trans_k = array(data = 0, dim = c(R, R, n))    # intialization
  stored_loglik = vector()
  
  # regression matrix
  t=seq(0,m-1)
  
  phi = designmatrix(t,p) # pour 1 courbe
  X= apply(phi,2,rep,n) # pour les n courbes
  
  Y=matrix(t(data),ncol=1)
  
  ## main algorithm
  try_EM = 0; 
  best_loglik = -Inf
  cputime_total = vector()
  
  while (try_EM < total_EM_tries) {
    try_EM = try_EM +1
    cat("EM_MixFHMMR try n° ",try_EM,"\n")
    time = Sys.time()
    
    ###################
    #  Initialization #
    ###################
    
    param = init_MixFHMMR(data, K, R, X, variance_type, ordered_states, init_kmeans, try_EM)
    
    browser()
    
    iter = 0
    converge = 0
    loglik = 0
    prev_loglik=-Inf
  
  # # EM ####  
  while ((!converge) && iter <= max_iter_EM) {
    #
    exp_num_trans_ck = array(data = 0, dim = c(R, R, n))
    exp_num_trans_from_l_cg = matrix(0,R,n)
    #
    exp_num_trans = array(data = 0, dim = c(R, R, n, K))
    exp_num_trans_from_l = array(data = 0, dim = c(R, n, K))
    #
    #w_k_fyi = zeros(n,K)
    log_w_k_fyi = matrix(0,n,K) 
    
    ##########
    # E-Step #
    ##########
    gamma_ikjr = array(data = 0, dim = c(n*m, R, K))

    
    for (k in 1:K) {
      # run a hmm for each sequence
      #     
      
      Li = matrix(0,n,1)
      for (i in 1:n) {
        log_fkr_yij = matrix(0,R,m)
        Y_i = data[i,] # ith curve
      
        for (r in 1:R) {
          beta_kr = param$beta_kr[,r,k]
          if (variance_type == "common") {
            sigma_kr = param$sigma_k[k]
            sk = sigma_kr
          } else {
            sigma_kr = param$sigma_kr[,k]
            sk = sigma_kr[r]
          }
          z=((Y_i - t(phi %*% beta_kr)) ^ 2) /sk
          log_fkr_yij[r,] = -0.5 * rep(1,3) * (log(2 * pi) + log(sk)) - 0.5 * z     
        }
        
        log_fkr_yij  = pmin(log_fkr_yij, log(.Machine$double.xmax))
        log_fkr_yij = pmax(log_fkr_yij ,log(.Machine$double.xmin))
        fkr_yij =  exp(log_fkr_yij)
        
        # forwards backwards ( calcul de logProb(Yi)...)
        param_fb = forwards_backwards(param$pi_k[,k], param$A_k[,,k], fkr_yij)
        
        Li[i] = param_fb$loglik   # loglik for the ith curve  ( logProb(Yi))
        gamma_ikjr[((i-1)*m+1) : (i*m),,k] = t(param_fb$tau_tk)   # [n*m R K] : "segments" post prob for each cluster k

        exp_num_trans_ck[,,i] = rowSums(param_fb$xi_tkl, dims=2)      # [R R n]
        exp_num_trans_from_l_cg[,i] = param_fb$tau_tk[,1]  # [R x n]
      }
      
      exp_num_trans_from_l[,,k] = exp_num_trans_from_l_cg   #[R n K]
      exp_num_trans[,,,k] = exp_num_trans_ck             #[R R n K]
      
      # for computing the global loglik
      log_w_k_fyi[,k] = log(param$w_k[k]) + Li   #[nx1]
    }
    
    log_w_k_fyi = pmin(log_w_k_fyi,log(.Machine$double.xmax))
    log_w_k_fyi = pmax(log_w_k_fyi,log(.Machine$double.xmin)) 
          
    tau_ik = exp(log_w_k_fyi) / (rowSums(exp(log_w_k_fyi)) %*% t(rep(1,K)))    #cluster post prob
    
    ## log-likelihood for the n curves
    loglik = sum(log(rowSums(exp(log_w_k_fyi))))
    
    ###########
    # M-Step  #
    ###########
    
    # Maximization of Q1 w.r.t w_k 
    param$w_k = t(colSums(tau_ik)) / n
    
    for (k in 1:K){
      if (variance_type == "common"){
        s=0
      }
      
      weights_cluster_k = tau_ik[,k]
      # Maximization of Q2 w.r.t \pi^g
      exp_num_trans_k_from_l =(rep(1,R) %*% t(weights_cluster_k)) * exp_num_trans_from_l[,,k]   #[R x n]
      param$pi_k[,k] = (1/sum(tau_ik[,k])) %*% rowSums(exp_num_trans_k_from_l)    # sum over i
      
      # Maximization of Q3 w.r.t A^g (the trans mat)
      
      for (r in 1:R){
        if (n == 1){
          exp_num_trans_k[r,,] = (rep(1,R) %*% t(matrix(weights_cluster_k))) * (exp_num_trans[r,,,k])
        } else {
          #exp_num_trans_k(k,:,:,g)
          exp_num_trans_k[r,,] = (rep(1,R) %*% t(matrix(weights_cluster_k))) * (exp_num_trans[r,,,k])
        }
      }
      if (n == 1){
        temp = exp_num_trans_k
      }else{
        temp = rowSums(exp_num_trans_k, dims=2)   # sum over i
      }

      
      param$A_k[,,k] = mk_stochastic(temp)
      # if HMM with order constraints
      if (ordered_states == 1){
        param$A_k[,,k] = mk_stochastic(param$mask * param$A_k[,,k])
      }
      
      # Maximisation de Q4 par rapport aux betak et sigmak
      Ng = colSums(tau_ik)  # nbr of individuals within the cluster k ,k=1...K estimated at iteration q
      ng = Ng[k] #cardinal nbr of the cluster k
      #each sequence i (m observations) is first weighted by the cluster weights

      weights_cluster_k =  apply(t(tau_ik[,k]),2,rep,m)
      weights_cluster_k = as.vector(weights_cluster_k)
      # secondly, the m observations of each sequance are weighted by the
      # wights of each segment k (post prob of the segments for each cluster g)
      gamma_ijk = gamma_ikjr[,,k]   # [n*m R]  
      nm_kr=colSums(gamma_ijk)   # cardinal nbr of the segments r,r=1,...,R within each cluster k, at iteration q              
      sigma_kr = matrix(0,R,1)
      
      beta_kr = cbind(array(data=0,dim =c(nrow(param$beta_kr),R)))
      
      for (r in 1:R){
        nmkr = nm_kr[r]    #cardinal nbr of segment r for the cluster k
        weights_seg_k = gamma_ijk[,r]
        Xkr = (sqrt(weights_cluster_k * weights_seg_k) %*% matrix(1,1,p+1)) * X   #[n*m x (p+1)]
        Ykr = (sqrt(weights_cluster_k * weights_seg_k)) * Y     #[n*m x 1] 
        #  Weighted least squares: maximization w.r.t beta_kr
        beta_kr[,r] = solve(t(Xkr) %*% Xkr , tol = 1e-30) %*% t(Xkr) %*% Ykr    # Maximisation par rapport aux betakr

         # # Maximization w.r.t sigmak :
        z = sqrt(weights_cluster_k * weights_seg_k) * (Y-X %*% beta_kr[,r])
        if (variance_type == "common"){
          s = s + t(z) %*% z
          ngm = sum(sum((weights_cluster_k %*%  matrix(1,1,R)) * gamma_ijk))
          sigma_k = s/ngm
        } else {
          ngmk = sum(weights_cluster_k * weights_seg_k)
          sigma_kr[r]=  t(z) %*%z /(ngmk)
        }
      }

      param$beta_kr[,,k] = beta_kr
      if (variance_type == "common"){
        param$sigma_k[k] = sigma_k
      } else {
        param$sigma_kr[,k] = sigma_kr
      }
    }

    iter=iter+1
    
    if (prev_loglik - loglik > threshold) {
      cat(sprintf("EM loglik is decreasing from %6.4f to %6.4f!\n",15.2,45.2),"\n")
    }
    if (verbose) {
      cat(sprintf("EM_MixFHMMR : Iteration : %d   log-likelihood : %f \n",  iter,loglik), "\n")
    }
    converge =  abs((loglik-prev_loglik)/prev_loglik) <= threshold
    
    if (is.na(converge)) {
      converge = FALSE
    }
    
    prev_loglik = loglik         
    stored_loglik[iter] = loglik
  } # end of EM  loop 
  
  cputime_total = cbind(cputime_total, Sys.time() - time)
  
  MixFHMMR =list()
  MixFHMMR$model = param
  if (variance_type == "common"){
    MixFHMMR$stats$paramter_vector = c(as.vector(param$w_k), as.vector(param$A_k), as.vector(param$pi_k),as.vector(param$beta_kr), as.vector(param$sigma_k))
  } else {
    MixFHMMR$stats$paramter_vector = c(as.vector(param$w_k), as.vector(param$A_k), as.vector(param$pi_k),as.vector(param$beta_kr), as.vector(param$sigma_kr))
  }
  MixFHMMR$stats$tau_ik = tau_ik
  MixFHMMR$stats$gamma_ikjr = gamma_ikjr
  MixFHMMR$stats$loglik = loglik
  MixFHMMR$stats$stored_loglik = stored_loglik 
  MixFHMMR$stats$log_w_k_fyi = log_w_k_fyi
  
  if (MixFHMMR$stats$loglik > best_loglik){
    best_loglik = MixFHMMR$stats$loglik
    best_MixFHMMR = MixFHMMR
  }
  
  if (try_EM >=1){
    cat(sprintf("log-lik at convergence: %2.2f", MixFHMMR$stats$loglik), "\n")
  }
  

  }   # Fin de la boucle sur les essais EM
  
  MixFHMMR$stats$loglik =   best_loglik
  if (try_EM>1){
    cat(sprintf("log-lik max: %2.2f", MixFHMMR$stats$loglik), "\n")
  }
  
  MixFHMMR = best_MixFHMMR
  # Finding the curve partition by using the MAP rule
  param_MAP =list()
  
  param_MAP = MAP(MixFHMMR$stats$tau_ik)    # MAP partition of the n sequences
  klas = param_MAP$klas
  Cig = param_MAP$Z
  MixFHMMR$stats$klas = klas
  
  # cas où on prend la moyenne des gamma_ijkr
  smoothed =  matrix(0,m,K)
  for (k in 1:K){
    betakr = MixFHMMR$model$beta_kr[,,k]
    weighted_segments = rowSums(MixFHMMR$stats$gamma_ikjr[,,k] * (X %*% betakr))
    
    dim(weighted_segments) = c(m,n)
    weighted_clusters = (matrix(1,m,1) %*% t(MixFHMMR$stats$tau_ik[,k])) * weighted_segments
    smoothed[,k] = (1 / sum(MixFHMMR$stats$tau_ik[,k])) %*% rowSums(weighted_clusters)    
  }
  MixFHMMR$stats$smoothed = smoothed 
  MixFHMMR$stats$cputime = mean(cputime_total)
  
  nu = length(MixFHMMR$stats$paramter_vector)
  
  # BIC AIC et ICL*
    MixFHMMR$stats$BIC = MixFHMMR$stats$loglik - (nu * log(n)/2)
    MixFHMMR$stats$AIC = MixFHMMR$stats$loglik - nu
  # ICL*             
  # Compute the comp-log-lik 
  cig_log_w_k_fyi = (Cig) * (MixFHMMR$stats$log_w_k_fyi)
  comp_loglik = sum(rowSums(cig_log_w_k_fyi))
  MixFHMMR$stats$ICL1 = comp_loglik - nu * log(n) / 2    
    
  return(MixFHMMR)
}
  
  
