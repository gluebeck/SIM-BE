TC1 = function(n,p,q,mu) {
  muR = mu/runif(n)
  ftimes = log((muR-q)/(muR-p))/(p-q)
  ftimes[is.na(ftimes)] = Inf
  # return(list(ftime = min(ftimes,na.rm=T),survivor=!is.na(ftimes)))
  return(list(ftimes = ftimes,survivor=is.finite(ftimes)))
}
