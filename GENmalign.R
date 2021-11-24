GENmalign = function(x=Malign,s=1,p=p0,aX=alphaXsi) {
  if(is.null(x) | s > 1) {
    warning('set is empty or time increment > 1')
    return()
  }
  s1  =   as.numeric(x[3])+s     # stop time (start + time increment)
  n0  =   as.numeric(x[4])       # no. of premalignant cells in clone
  # m   =   as.numeric(x[5])     # malignant tumor size (no. of cells)
  # s0  =   as.numeric(x[6])       # time of persistent malignant cell birth
  
  # n = no. of progenitor cells that give rise to clones with FPTs
  if(n0 < theta) {
    n1 = genSizeSurv(n=n0,p0=p,alphaXsi=aX)
    if(s < 1) {
      n1 = n0+round(s*(n1-n0))   # approximate simulated size for s < 1
    }
  } 
  else {                         # for n0 > theta, proceed deterministically
    n1=round(n0*exp(g.I*s))
  }
  
  # update DF at time s1 (includes deterministic malignant growth)
  x[3] = s1
  x[4] = n1
  x[5] = round(max(1,exp(g.M*(s1-as.numeric(x[6])))))-1 # deterministic growth of malign clone for speed 
  
  return(x)
}
