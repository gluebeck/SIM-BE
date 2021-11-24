byear = 1950
sex   = "M"

#### parameters for regularized APC model for EAC
nu =  0.78316E-04
y0 =  0.20087E+04
b1 =  0.92691E-03
b2 =  0.11170E-01
w1 =  0.15052E-02
w2 = -0.11711E-01
scl = 0.32849E+02

g0 =  0.16249E+00
g1 = -0.49090E-02
g2 =  0.12270E+00
y1 =  0.19417E+04

#### cell divisions
a0 = 10
alpha.M = 50

#### initiation
Xorig = 1.e6
mu0 = mu1 = 0.10703E-01

#### cell division rates (cohort adjusted)
dby = byear-y1
gfac = exp(g1*dby*(1+g2*dby))
alpha.I = alpha = a0*gfac # cell division rate adjusted for birth cohort

#### malignant transformation
mu2 = mu = 0.29868E-06*gfac
tlag = 5

