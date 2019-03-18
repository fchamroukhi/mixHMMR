normalize =function(A, dimension)
{
  # NORMALISE Make the entries of a (multidimensional) array sum to 1
  # [M, c] = normalise(A)
  # c is the normalizing constant
  #
  # [M, c] = normalise(A, dimension)
  # If dimension is specified, we normalise the specified dimension only,
  # otherwise we normalise the whole array.  
  
  if (missing(dimension)) {
    z=sum(A)
    # Set any zeros to one before dividing
    # This is valid, since c=0 => all i. A(i)=0 => the answer should be 0/1=0
    s = z + (z==0)
    M = A / s
    
  }else if (dimension==1) {
    z = sum(A)
    s = z + (z==0)
    M = A / apply(s, 2, rep, dim(A)[1])
  }
  else {
    z=rowSums(A)
    s = z + (z==0)

    L=dim(A)[dimension]
    d=length(dim(A))
    v=rep(1,d)
    v[dimension]=L
    c=t(apply(t(s),2,rep,v[2]))
    M=A/c
  }
  
  to_return = list()
  #to_return$c = c
  to_return$M = M
  to_return$z = z
  return(to_return)
}