genSizeSurv = function(n,p0,alphaXsi) {
  # dum = 1-alphaXsi*(1-p0)
  R = runif(n,min=p0,max=1)    # conditional on non-extinction, so n >=1
  out = ceiling(log((1-R)/(1-p0))/log(alphaXsi))    # marginal distribution (mu=0)
  return(sum(out))
}
