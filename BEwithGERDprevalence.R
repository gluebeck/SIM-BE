parm = as.data.frame(read.table(file='../bhat.inp',header=F))
parm = as.data.frame(read.table(file='../bhat_conv3s_EAC2_AM.out',header=F))
parm = as.data.frame(read.table(file='../bhat_conv3s_EAC2_sfunc1_AM.out',header=F))
dimnames(parm)[[2]] = c('name','value','status','dum','lower','upper')
dum = dimnames(parm)[[1]] = as.character(parm$name)
parm = parm[,-c(1,3:6)]
names(parm) = dum

nu = parm['nu']
b1 = parm['b1']
b2 = parm['b2']
# b3 = parm['b3']
g0 = parm['g0']
g1 = parm['g1']
g2 = parm['g2']
# g3 = parm['g3']
yh = y0 = parm['yh']
c0 = parm['c0']

scl = parm['scl']
w1 = parm['w1']
w2 = parm['w2']
g3=b3=0

BEonset.max = 90  # used for conditioning on T_BE < t_screening
GERD = FALSE
byear = 1950

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

#### birth cohort effect on nu
# byear = 1900:2000

dby = byear-y0
byear.effect = b1*dby*dby + b2*dby^3 + b3*dby^4 ; byear.effect[dby>0]=0
periodfac = exp(byear.effect) 

dum = years-y0; dum[dum>0]=0
sfunc = log(age/scl)
periodfac = periodfac*exp(sfunc*w1*dum*dum*(1+w2*dum))

plot(age,nu*periodfac,type='l',lwd=2)
lines(age,nu*periodfac,lty=5,lwd=2,col=4)

# g
dum = byear - c0
gfac = g0*exp(g1*dum+g2*dum*dum+g3*dum*dum*dum)
plot(byear,gfac,type='l',lwd=2,ylim=c(0,1.2))
lines(byear,gfac,lwd=2,lty=1,col=2)

# x = years-y0
# #periodfac = exp(log(age/scl)*w1*x*x*(1+w2*x)+byear.effect) 
# #periodfac = exp(w1*x*x*(1+w2*x)+byear.effect) 
# periodfac[x>0]=exp(byear.effect)

#### generate N BE onsets

#### GERD factor
# if(GERD==TRUE) {gerdfac = rrgerd*gerdpr/gpop} else if(GERD==FALSE) {gerdfac = (1-gerdpr)/gpop}
#### uncomment the next line for entire population BE prob
gerdfac0 = 1/(1-gerdpr + rrgerd*gerdpr)
gerdfac = 1

BErate0 = nu*gerdfac0*periodfac   #conditional rate of BEonset <= BEonset.max 
BErate = nu*periodfac   #conditional rate of BEonset <= BEonset.max 
BEprob0 = 1-exp(-cumsum(BErate0))
BEprob = 1-exp(-cumsum(BErate))

BEonsets = approx(x=BEprob,y=age,xout=sort(runif(N,0,max(BEprob)), decreasing = TRUE),
                  yleft=0,yright=BEonset.max)


#### compute BE prevalence among patients with GERD
len = length(age)
fgerd = c(0,gerdpr[2:len]-gerdpr[1:(len-1)])

BEprobGERD = numeric(len-1)
for(s in age[-1]) {
  dum = BErate0[1:s]
  x =0 
  for(w in 1:s) {
    dum[w:s]=BErate0[w:s]*rrgerd
    x=x+fgerd[w]*(sum(dum[1:(w-1)])+sum(dum[w:s]))
  }
  BEprobGERD[s]=x/gerdpr[s]; BEprobGERD[1]=0
}

plot(age,rrgerd*BErate0,type='l',cex.lab=1.3,lwd=3,lty=1,ylab='BE rate',col=2)
lines(age,BErate,col=4,lty=1,lwd=2)
lines(age,BErate0,col=1,lty=1,lwd=2)
# lines(1:90,BEprobGERD,lwd=1,col=4)
legend('topleft',c('projected background','population','with GERD'),col=c(1,4,2),lwd=c(2,2,2),bty='n',y.intersp = 1.2)

plot(age,BEprob0,type='l',cex.lab=1.3,lwd=3,lty=1,ylab='BE prevalence',ylim=c(0,.08))
lines(age,BEprob,col=4,lty=1,lwd=2)
lines(1:90,BEprobGERD,lwd=1,col=4)
legend('topleft',c('projected background','population','with GERD'),col=c(1,4,2),lwd=c(2,2,2),bty='n',y.intersp = 1.2)

plot(1:40,fgerd[1:40]/gerdpr[40],type='l',cex.lab=1.3,xlab='age',ylab='GERD density function',lwd=2,lty=1,xlim=c(0,70))
lines(1:50,fgerd[1:50]/gerdpr[50],lwd=2,lty=2)
lines(1:60,fgerd[1:60]/gerdpr[60],lwd=2,lty=3)
lines(1:70,fgerd[1:70]/gerdpr[70],lwd=2,lty=4)

