#### 3-stage hazard function

age = seq(tlag, MaxAge, .1)
haz_3s = mu0*X*(1 - ((q-p)/(q*exp(-p*(age-tlag)) -p*exp(-q*(age-tlag))))^(mu1/alpha.I))
ptu_3s = 1-exp(-.1*cumsum(haz_3s))
plot(age,ptu_3s,type='l',lwd=2)
plot(age,haz_3s,type='l',lwd=2)
# lines(age,ptu_3s,lty=3,lwd=2)

dum = as.data.frame(Status)

f10 = sum(dum$EACpreScreen[dum$age==10]==T)/N
f20 = sum(dum$EACpreScreen[dum$age==20]==T)/N
f30 = sum(dum$EACpreScreen[dum$age==30]==T)/N
f25 = sum(dum$EACpreScreen[dum$age==25]==T)/N
f40 = sum(dum$EACpreScreen[dum$age==40]==T)/N
points(c(10,20,25,30,40),c(f10,f20,f25,f30,f40))
points(c(10,20,30),c(f10,f20,f30))

##### evaluate findings LGD, HGD etc
dum = as.data.frame(Status)

cond = dum$age==firstScreen
sum(dum$LGD[cond]==T)/N
sum(dum$HGD[cond]==T)/N
sum(dum$IMC[cond]==T)/N
sum(dum$InSitu[cond]==T)/N
sum(dum$EAC[cond]==T)/N
sum(dum$EACpreScreen[cond]==T)/N

cond = dum$age==65
N60 = sum(cond)
sum(dum$LGD[cond]==T)/N60
sum(dum$HGD[cond]==T)/N60
sum(dum$IMC[cond]==T)/N60
sum(dum$InSitu[cond]==T)/N60
sum(dum$EAC[cond]==T)/N60
sum(dum$EACpreScreen[cond]==T)/N60

cond = dum$age==70
N70 = sum(cond)
sum(dum$LGD[cond]==T)/N70
sum(dum$HGD[cond]==T)/N70
sum(dum$IMC[cond]==T)/N70
sum(dum$InSitu[cond]==T)/N70
sum(dum$EAC[cond]==T)/N70
sum(dum$EACpreScreen[cond]==T)/N70

cond = dum$age==80
N80 = sum(cond)
sum(dum$LGD[cond]==T)/N80
sum(dum$HGD[cond]==T)/N80
sum(dum$IMC[cond]==T)/N80
sum(dum$InSitu[cond]==T)/N80
sum(dum$EAC[cond]==T)/N80
sum(dum$EACpreScreen[cond]==T)/N80

#### Ymalign and Ybenign trajectories
dumY = as.matrix(Ymalign)
dumS = as.matrix(Smalign)
plot(dumS[1,],dumY[1,],type='s',lwd=1,col='pink',xlab ='BE age',ylab='size',log='y',cex.lab=1.3,
     xlim=c(0,42),ylim=c(1,100000))
for (i in 2:nrow(dumS)) lines(dumS[i,],dumY[i,],col='pink',lwd=1,type='s')

ids = sample(1:nrow(dumS),20)
for (i in 1:20) lines(dumS[ids[i],],dumY[ids[i],],col=2,lwd=2,type='s')

dumY = as.matrix(Ybenign)
dumS = as.matrix(Sbenign)
dumY[dumY==0] = .5
ids = sample(1:nrow(dumS),20)
for (i in 1:20) lines(dumS[ids[i],],dumY[ids[i],],col=1,lwd=2,type='s')



x = dum$BEonset
y = dum$ageClinicalEAC

plot(x,y,pch=19,cex.lab=1.3,xlim=c(20,MaxAge),xlab='age of BE onset',ylab='age EAC is detected')
abline(0,1)

points(x,y,pch=19,col=2)

####  BE at age 60 stats
dum = as.data.frame(Status)
dum1 = dum[which(dum$age==60),]
sum(dum1$ageClinicalEAC < 60)/N

n60 =  sum(dum1$ageClinicalEAC<60 | dum1$LGD==T | dum1$HGD == T | dum1$IMC == T | dum1$InSitu == T)
n70 =  sum(dum$age==70 & dum1$LGD==F & dum1$HGD == F & dum1$IMC == F & dum1$InSitu == F)
n60 = sum(BEex60)
(sum(dum1$ageClinicalEAC < 70) - n60)/(N-n60)/10  #in 10 years
(sum(dum1$ageClinicalEAC < 80) - n60)/(N-n60)/20  #in 20 years



n70 = sum(dum1$ageClinicalEAC < 70)
(sum(dum1$ageClinicalEAC < 80) - n70)/(N-n70)



dum1 = dum[which(dum$age==70),]
sum(dum1$ageClinicalEAC<70)

dum1 = dum[which(dum$age==80),]
sum(dum1$ageClinicalEAC<80)

sum(dum$LGD == T)
sum(dum$HGD == T)
sum(dum$IMC == T)
sum(dum$InSitu == T)
sum(dum$EAC == T)

sum(dum$ageClinicalEAC<60)
sum(dum$ageClinicalEAC<70)
sum(dum$ageClinicalEAC<80)

##########################################################################################
#### plot Y and Z sizes against Y start times
ids = order(Benign$start)
x0 = Benign$start[ids]; x0 = x0+runif(length(x0)) 
Y0 = Benign[ids,"Y"]

ids = order(Malign$start)
x = Malign$start[ids]; x = x+runif(length(x)) 
Y=Malign[ids,"Y"]
Z=Malign[ids,"Z"]

plot(x,Y,pch=19,xlim=c(0,45),cex=.8,cex.lab=1.3,xlab='origin of clone (years)',col=2,log='y')
points(x0,Y0,cex=.8)
legend('bottomleft',c('non-progressive (mean: 464)','progressive (mean: 1088)'),pch=c(1,19),col=c(1,2),bty='n',y.intersp = 1.2)

ids = which(Z>0)
points(x[ids],Y[ids],cex=2,col=4)
points(x[ids],Y[ids],cex=2.2,col=4)

##########################################################################################
#### size distribution of Y and Z at firstScreen
plot(density(Y.initiated.malignant,bw = 10),main='',lwd=2,cex.lab=1.3,xlab='clone size',log='x',col=2)
lines(density(Y.initiated.benign,bw=10),col=1,lwd=2)

##########################################################################################
#### prevalences over time for BE at age 40 - no interventions

LGD = HGD = IMC = InSitu = EAC = numeric(7)
tmp = as.data.frame(Status)
N.AR =2000
for (i in 1:7) {
  dat=tmp[tmp$age==BEscreenings[i],]
  LGD[i] = sum(unlist(dat$LGD))/N.AR
  HGD[i] = sum(unlist(dat$HGD))/N.AR
  IMC[i] = sum(unlist(dat$IMC))/N.AR
  InSitu[i] = sum(unlist(dat$InSitu))/N.AR
  EAC[i] = sum(unlist(dat$EACpreScreen))/N.AR
}

plot(smooth.spline(BEscreenings,LGD),lwd=2,lty=3,type='l',cex.lab=1.3,xlab='age',ylab='prevalence',
     xlim=c(40,72),ylim=c(0,.6))
lines(smooth.spline(BEscreenings,HGD),lwd=2,lty=2,)  
lines(smooth.spline(BEscreenings,IMC+InSitu),lwd=2,lty=1,)  
lines(smooth.spline(BEscreenings,EAC),lwd=3,lty=1,col=2)  

legend('topleft',c('EAC','IMC+InSitu','HGD',"LGD"),lwd=c(3,2,2,2),lty=c(1,1,2,3),col=c(2,1,1,1),bty='n',y.intersp = 1.2)

##########################################################################################
#### prevalences over time for GERD patients

len = length(BEscreenings)
EACprop = EACinc = numeric(len)
for (i in 1:len) {
  dat=Status[Status$age==BEscreenings[i],]
  
  EACprop[i]=sum(dat$EACpreScreen)/N.AR
  if(i > 1) {
    EACinc[i] = 1.e5*(EACprop[i] -EACprop[i-1])/5 } else {
      EACinc[1] = 1.e5*EACprop[1]/2/BEscreenings[1]
    }
}

plot(smooth.spline(BEscreenings,popfracBE*EACprop,df=5),lwd=2,lty=1,type='l',cex.lab=1.3,xlab='age',ylab='EAC prevalence',
     xlim=c(40,82))

plot(smooth.spline(BEscreenings,popfracBE*EACinc,df=5),lwd=2,lty=1,type='l',cex.lab=1.3,xlab='age',ylab='EAC prevalence',
     xlim=c(40,82))


done:  1000 
done:  2000 
system time info 156.771 15.967 172.817 0 0 

MaxAge:  30   BEScreen:  30 
LGD:  0.0785 
HGD:  0.0135 
IMC:  0.012 
InSitu:  0.067 
EAC:  0.2495 
preScrEAC:  0.252 
postScrEAC:  0.748 