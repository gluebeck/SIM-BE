#### stochastic-to-determisistic threhold
theta = 1000
Threshold.InSitu = 1000000
Threshold.HGD = 50000

#### HGD conversion process with per cell rate (per year) rho.HGD applied to 1 year increments
rho.HGD = 1 - (0.5)^(1/Threshold.HGD)  # 50% conversion prob when Y-clone is Threshold.HGD cells

#### Cancer observation process with per cell rate rho
rho.InSitu = 1 - (0.5)^(1/Threshold.InSitu)   # for 50% prob when malignant clone is 1.e9 cells

#### Cancer observation process with per cell rate rho
rho.EAC = 1 - (0.5)^(1/1.e9)   # for 50% prob when malignant clone is 1.e9 cells

#### promotion of intermediate cells
g.I = g0*gfac
gamma.I = 1-(g.I+mu)/alpha.I  # ratio  death/birth (d/b)

#### malignant growth
# g.M = log(1.e9)/5  # 10^9 cells on average in clinical tumor - 5 year lag time
g.M = -log(alpha.M*rho.EAC)/tlag
for (i in 1:10) g.M = -log(alpha.M*rho.EAC/g.M/g.M)/tlag
g.M = 1.*g.M  #empirical adjustment for agreement with clinical cancers at FPT+tlag
gamma.M = 1-g.M/alpha.M  #assuming g.M = 1


#### scale up /down X to for Dysplasia prevalence on biopsy screening
#### invariant slope: X*mu0*mu1 
Xfac = 1
X = Xfac * Xorig
mu0 = mu1 = sqrt(1/Xfac)*mu0 # again, assuming mu0=mu1

#### identifiable parameter combinations
rootarg = g.I*g.I+4*alpha.I*mu
p = 0.5*(-g.I-sqrt(rootarg))
q = 0.5*(-g.I+sqrt(rootarg))

#### combinations for stochastic conversion from LGD to HGD
rootarg = g.I*g.I+4*alpha.I*rho.HGD
p.HGD = 0.5*(-g.I-sqrt(rootarg))
q.HGD = 0.5*(-g.I+sqrt(rootarg))
