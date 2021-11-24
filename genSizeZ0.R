genSizeZ0 = function(n,p0,alphaXsi) {
  if(n < 1) {out = 0} else {
    R = runif(n)
    out = ceiling(log(R/(1-p0))/log(alphaXsi))    # sample sizes cond. on no new mutations
    out[out<0] = 0    # if negative, a clone has gone to extinction
  }
  return(sum(out))
}
