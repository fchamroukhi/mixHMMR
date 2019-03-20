show_MixRHLP_results = function(data, mixFHMMR) {
  

  n = nrow(data)
  m = ncol(data)
  t = 0:(m - 1)
  
  G = length(mixFHMMR$model$w_k)
  
  couleur = rainbow(G)
  matplot(t, t(data), type = "l", lty = "dotted", xlab = "Time", ylab = "y")
  
  for (k in 1:K) {
    cluster_k = data[mixFHMMR$stats$klas==k,]
    matplot(t, t(cluster_k), col = couleur[k], type = "l", lty = "dotted", xlab = "Time", ylab = "y")
    lines(t, mixFHMMR$stats$smoothed[,k], lty = "solid", lwd = 2, col = "black")
    title(main = sprintf("Cluster %1.1i", k))
  }
  
}