designmatrix = function(x,p) 
{

  # function X = designmatrix(x,p)
  # constructs the design matrix of a polynomial regression of degree p
  #
  #
  # Faicel Chamroukhi
  ###########################################################################  

  
  if (!is.vector(x)) {
    x = as.vector(t(x))
  }
n = length(x)
X= matrix(0,nrow=n,ncol=p+1)

for (i in seq(1, p+1)) 
{
  print(i)
  X[,i] = x^(i-1)
}

return(X)
}
           
  