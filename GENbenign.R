GENbenign = function(x=Benign,s=1,p=p0,aX=alphaXsi) {
  if(is.null(x) | s > 1) {
    warning('set is empty or time increment > 1')
    return()
  }
  s1  =   as.numeric(x[3])+s     # stop time (start + time increment)
  n0  =   as.numeric(x[4])       # no. of premalignant cell in clone
  # n = no. of progenitor cells that DONOT give rise to clones with FPTs
  
  n1 = genSizeZ0(n0,p,aX)
  if(s < 1) {
    n1 = n0+round(s*(n1-n0))
  }
  
  x[3] = s1
  x[4] = n1
  return(x)
}
