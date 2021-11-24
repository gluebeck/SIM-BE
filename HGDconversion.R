HGDconversion = function(n,p,q,rho) {
  rhoR = rho/runif(n,min=0,max=rho/q)
  ftimes = log((rhoR-q)/(rhoR-p))/(p-q)
  ftimes[is.na(ftimes)] = Inf
  # return(list(ftime = min(ftimes,na.rm=T),survivor=!is.na(ftimes)))
  return(ftimes)
}
