#### generates BE onsets for N individuals
# BEonset.max = defined upfront for conditioning (death age or screening age)
GERD = TRUE

yrsmax = byear+BEonset.max
years = byear:yrsmax # period
age = years-byear #    s=age

#### GERD prevalence
if(sex=="F") {
  a1 = 0.18012E-03
  a2 = 0.71326E-02
  a3 = 0.25808E+02
  rrgerd = 5.0E0} else {
    a1 = 0.61422E-03
    a2 = 0.70447E-02
    a3 = 0.26002E+02
    rrgerd = 5.0E0
  }

gerdpr = numeric(length(age))
#### GERD prevalence
for (i in 1:length(age)) {
  gerdpr[i] =   1 - exp( -a1 * min(a3, age[i]) - a2 * max(0, age[i] - a3))
}
gpop = (1 - gerdpr + gerdpr * rrgerd)

#### generate N BE onsets #################################
#### birth cohort effect on nu
dby = byear-y0
byear.effect = b1*dby*dby*(1+b2*dby); byear.effect[dby>0]=0

x = years-y0; x[x>0]=0
periodfac = exp(log(age/scl)*w1*x*x*(1+w2*x)+byear.effect) 
#periodfac = exp(w1*x*x*(1+w2*x)+byear.effect) 

#### GERD factor
if(GERD==TRUE) {gerdfac = rrgerd*gerdpr/gpop} else if(GERD==FALSE) {gerdfac = (1-gerdpr)/gpop}
#### uncomment the next line for entire population BE prob
gerdfac = rrgerd/(1-gerdpr + rrgerd*gerdpr)

BErate = nu*gerdfac*periodfac   #conditional rate of BEonset <= BEonset.max 
BEprob = 1-exp(-cumsum(BErate))

BEonsets = approx(x=BEprob,y=age,xout=sort(runif(N,0,max(BEprob)), decreasing = TRUE),
                  yleft=0,yright=BEonset.max)

# uncomment next line if testing BE onset==0
# BEonsets$y = rep(0,N)

plot(age,BEprob,type='l',cex.lab=1.3,lwd=3,lty=1,ylab='BE prevalence',cex.lab=1.3,ykim=c(0,.08))
# lines(age,BEprob,lwd=1,lty=1)
# legend('topleft',c('BE prevalence (general pop)','BE prevalnce (GERD)',
#                   'BE prevalence (non-GERD)','GERD prevalence/25'),bty='n',
#                  lwd=rep(2,4),lty=c(1,1,3,1),col=c(1,1,1,2),y.intersp = 1.2)