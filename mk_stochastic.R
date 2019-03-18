mk_stochastic = function(T)
  {
  # MK_STOCHASTIC Ensure the argument is a stochastic matrix, i.e., the sum over the last dimension is 1.
  # T = mk_stochastic(T)
  #
  # If T is a vector, it will sum to 1.
  # If T is a matrix, each row will sum to 1.
  # If T is a 3D array, then sum_k T(i,j,k) = 1 for all i,j.
  
  # Set zeros to 1 before dividing
  # This is valid since S(j) = 0 iff T(i,j) = 0 for all j
  
  if ((length(dim(T)) == 2) & ( dim(T)[1] == 1 | dim(T)[2] == 1)) {    #isvector
    T = normalize(T)
  } else if (length(dim(T)) == 2){
    S = colSums(T)
    S = S +(S==0)
    norm1 = matrix(rep(S, dim(T)[2]), nrow = length(S))
    T = T / norm1
  }   
  else {    # multi-dimensional array
    ns = dim(T)
    T = matrix(T, nrow = prod(ns(1 : (length(ns)-1))), ncol = ns(length(ns)))
    S = rowSums(T)
    S = S +(S==0)
    norm1 = matrix(rep(S, ns(length(ns))), nrow = nrow(S))
    T = T / norm1
    dim(T) = c(ns)
  }
  return(T)
}