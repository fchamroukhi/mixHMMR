forwards_backwards = function(prior, transmat, f_tk)
{
  ####################################################################################################
    # function [tau_tk, xi_tlk, alpha_tk, beta_tk, loglik, xi_summed] = forwards_backwards(prior, transmat, f_tk)
    # forwards_backwards : calculates the E-step of the EM algorithm for an HMM
    # (Gaussian HMM)
    
    # Inputs :
      #
    #         prior(k) = Pr(z_1 = k)
    #         transmat(\ell,k) = Pr(z_t=k | z_{t-1} = \ell)
    #         f_tk(t,k) = Pr(y_t | z_y=k;\theta) #gaussian
    #
    # Outputs:
      #
    #        tau_tk(t,k) = Pr(z_t=k | X): post probs (smoothing probs)
    #        xi_tk\elll(t,k,\ell)  = Pr(z_t=k, z_{t-1}=\ell | Y) t =2,..,n
    #        with Y = (y_1,...,y_n);
    #        alpha_tk(k,t): [Kxn], forwards probs: Pr(y1...yt, zt=k)
    #        beta_tk(k,t): [Kxn], backwards probs: Pr(yt+1...yn|zt=k)
    #
    #
    #
    # Faicel Chamroukhi
    #######################################################################################################  
  
  filter_only = 0
  T = dim(f_tk)[2]
  K = length(prior) 
  
  if (dim(t(t(prior)))[2] != 1) {
    prior = t(prior)
  } 
  scale1 = rep(1,T) # pour que loglik = sum(log(scale)) part de zero

  tau_tk = matrix(0,K,T)
  xi_tkl = array(data = 0, dim = c(K,K,T-1)) 
  xi_summed = matrix(0,K,K)
  alpha_tk = matrix(0,K,T)
  beta_tk = matrix(0,K,T)  
    
  # forwards: calculate the alpha_tk's
  t=1
  alpha_tk[,1] = prior * f_tk[,t]
  param_norm = normalize(alpha_tk[,t])
  alpha_tk[,t] = param_norm$M
  scale1[t]= param_norm$z
  
  for (t in 2:T){
    param_norm = normalize((t(transmat) %*% alpha_tk[,t-1]) * f_tk[,t])
    alpha_tk[,t] = param_norm$M
    scale1[t]= param_norm$z
  }
  loglik = sum(log(scale1))
  
  beta_tk[,T] = rep(1,K)
  param_norm = normalize(alpha_tk[,T] * beta_tk[,T])
  tau_tk[,T]=param_norm$M
  
  for (t in seq(T-1,1)){
    param_norm1 = normalize(transmat %*% (beta_tk[,t+1] * f_tk[,t+1]))
    beta_tk[,t] = param_norm1$M
    param_norm2 = normalize(alpha_tk[,t] * beta_tk[,t])
    tau_tk[,t] = param_norm2$M
    param_norm3 = normalize((transmat * (alpha_tk[,t] %*% t(beta_tk[,t+1] * f_tk[,t+1]))))
    xi_tkl[,,t] =  param_norm3$M
    xi_summed = xi_summed + sum(xi_tkl[,,t],3)
  }
  
  to_return = list()
  to_return$tau_tk = tau_tk
  to_return$xi_tkl = xi_tkl
  to_return$alpha_tk = alpha_tk
  to_return$beta_tk = beta_tk
  to_return$loglik = loglik
  to_return$xi_summed = xi_summed
  
  return(to_return)
}