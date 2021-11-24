byear = 1950
sex   = "M"

#### parameters for regularized APC model for EAC
nu =  0.27405E-04
y0 =  0.20327E+04
b1 =  0.37547E-03
b2 =  0.70670E-02
w1 =  0.22038E-02
w2 =  0.15532E-01
scl = 32

g0 =  0.34586E+00
g1 = -0.41138E-01 
g2 = -0.69307E-02  
y1 =  0.19417E+04

#### cell divisions
a0 = 10
alpha.M = 50

#### initiation
Xorig = 1.e6
mu0 = mu1 = 0.31889E-02

#### cell division rates (cohort adjusted)
dby = byear-y1
gfac = exp(g1*dby*(1+g2*dby))
alpha.I = alpha = a0*gfac # cell division rate adjusted for birth cohort

#### malignant transformation
mu2 = mu = 0.32405E-06*gfac
tlag = 5

