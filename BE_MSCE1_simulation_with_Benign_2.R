#### number of BE individuals to simulate
# library(dqrng)
dqRNGkind("Xoroshiro128+")
dqset.seed(43)

N=1
MaxAge = 40
Tdeath = 40
BEonset.max = min(MaxAge,Tdeath) # should be Tdeath unless conditioning on BE presence at a given screening age

# screening times (array with integer years)
BEscreenings = 1:40
firstScreen = min(BEscreenings)
lastScreen = max(BEscreenings)

AtRisk = rep(N,length(BEscreenings))
Status =  data.frame(ID=NULL,BEonset=NULL,age=NULL,ageClinicalEAC=NULL,LGD=NULL,
      HGD=NULL,IMC=NULL,InSitu=NULL,EAC=NULL,EACpreScreen=NULL,EACpostSCreen=NULL)

source('BE_MSCE1_modelParameters_A.R')
source('BE_MSCE1_Parameters.R')
source('BE_MSCE1_BEonset.R')
source('BE_MSCE1_shared.R')
source('BE_MSCE1_initiation_with_Benign.R')
source('BE_MSCE1_status.R')
source('BE_MSCE1_detectProb.R')
source('TC1.R')
source('HGDconversion.R')

info = system.time(
  for (n in 1:N){
    
    BEonset = BEonsets$y[n]   #needs to be integer for synchronization
    BEage = MaxAge-BEonset
    # BEage = lastScreen-BEonset

    #############################################################################
    # define data frame that includes all intended screens for this individual 
    ns = length(BEscreenings)
    datf = data.frame(ID=rep(n,ns),BEonset=rep(BEonset,ns),age=BEscreenings,ageClinicalEAC=rep(Inf,ns),
          LGD=rep(FALSE,ns),HGD=rep(FALSE,ns),IMC=rep(FALSE,ns),InSitu=rep(FALSE,ns),EAC=rep(FALSE,ns),
          EACpreScreen = rep(FALSE,ns),EACpostScreen=rep(FALSE,ns))

    
    #############################################################################
    # TwoMutation Initiation from BEonset to BEonset.max
    firstScreen = min(BEscreenings[BEscreenings > BEonset])
    # only update Status if BE is found at a screening
    if(firstScreen < Inf) {

      tmax = ceiling(lastScreen - firstScreen)
      Init = TwoMutationInitiation(BEage,BEonset,firstScreen)

      s.initiated.benign = Init$BenignS
      N.initiated.benign = Init$BenignN
      Y.initiated.benign = Init$BenignY
      
      s.initiated.malignant = Init$MalignS
      N.initiated.malignant = Init$MalignN
      Y.initiated.malignant = Init$MalignY
      Z.malignant = rep(0,N.initiated.malignant)
      s.malignant = Init$FPTs

      time.to.firstScreen = firstScreen-(BEonset+s.malignant) # for initiated clones that turn malignant
      cond = time.to.firstScreen > 0
      Z.malignant[cond] = exp(g.M*time.to.firstScreen[cond])  # give malignancies deterministic sizes / can be modified 
    
      #############################################################################
      # dataframe for Y clones that DO NOT give rise to malignant clones (extinction) - setup
      Benign = data.frame(BEage=rep(tmax,N.initiated.benign),
                          start=s.initiated.benign,s=s.initiated.benign,
                          Y=Y.initiated.benign,Z=rep(0,N.initiated.benign),
                          fpT=rep(Inf,N.initiated.benign))
      #############################################################################
      # dataframe for Y clones that give rise to a malignant clone (non-extinction) - setup
      Malign = data.frame(BEage=rep(tmax,N.initiated.malignant),
                        start=s.initiated.malignant,s=s.initiated.malignant,
                        Y=Y.initiated.malignant,Z=Z.malignant,fpT=s.malignant,
                        LGD  = rep(TRUE,N.initiated.malignant),
                        HGD  = rep(FALSE,N.initiated.malignant),
                        IMC  = rep(FALSE,N.initiated.malignant),
                        EAC  = rep(FALSE,N.initiated.malignant))

      #############################################################################
      # starting at firstScreen, determine true (underlying) state of each FTP and non-FTP clone 
      cond2 = dqrunif(length(Malign$Z)) < 1 - (1-rho.EAC)^Malign$Z
      cond3 = dqrunif(length(Malign$Z)) < 1 - (1-rho.InSitu)^Malign$Z
      cond4 = Init$MalignU < firstScreen         # if True LGD status has changed to HGD before firstScreen
      Malign$LGD[cond4] = FALSE
      Malign$HGD[cond4] = TRUE
      Malign$EAC[cond2] = TRUE
      Malign$IMC[Malign$Z > 0 & !cond3] = TRUE  # any cancer less than InSitu transition. Perhaps redefine!!!!
    
      # cond1 = dqrunif(length(Malign$Y)) < 1 - (1-rho.HGD)^Malign$Y
      # Malign$LGD[cond1] = FALSE
      # Malign$HGD[cond1] = TRUE
      #############################################################################
    
      # InitialStatus = BEstatus(BEonset=BEonset,runT=0,Benign=Benign,Malign=Malign)
      cond.stop = FALSE
      IndicatorFirstScreen = TRUE
      runT = 0
    
      # loop over all initiated clones and evolve them incrementally, step size s
      # start time: first screen
      # can be stopped when clinical EAC occurs
    
      while(runT <= tmax) { #} & cond.stop == FALSE) {

        age   = firstScreen  + runT

        # record BE status when screening - and include other actions (treatment ...)
        if(age%in%BEscreenings==TRUE) {
        
          #### update true (underlying states) for subsequent screens
          if(IndicatorFirstScreen == FALSE) {
            sinceBEonset = age-BEonset
            cond2 = Malign$s <= sinceBEonset  # present or coming up
            cond3 = Init$MalignU < age        # conversion occured before age
            cond4 = dqrunif(length(Malign$Z[cond2])) < 1 - (1-rho.EAC)^Malign$Z[cond2]
          
            # cond1 = Benign$s <= sinceBEonset  # present or coming up
            # Benign$LGD[cond1][cond]  = FALSE  #### approximate all Benign lesions as LGD and unobservable
            # Benign$HGD[cond1][cond]  = TRUE 

            Malign$LGD[cond3]         = FALSE
            Malign$HGD[cond3]         = TRUE
            Malign$EAC[cond2][cond4]  = TRUE

            if (n == 1) {
              Ybenign = cbind(Ybenign,Benign$Y)
              Sbenign = cbind(Sbenign,Benign$s)   
              Ymalign = cbind(Ymalign,Malign$Y)  
              Smalign = cbind(Smalign,Malign$s)
            }
        } else {
            if (n == 1) {
              Ybenign = Benign$Y   # Y.initiated.benign
              Sbenign = Benign$s   # s.initiated.benign
              Ymalign = Malign$Y   # Y.initiated.malignant
              Smalign = Malign$s   # s.initiated.malignant
            }
          }
          
          #### update Status (kept in datf)
          tmp = BEstatusByBiopsy(ID=n,age=age,BEonset=BEonset,Malign=Malign)
          datf[datf$age ==age,] = unname(as.data.frame(tmp))

          # optional condition to remove LGD, HGD, IMC/InSitu and EACs patients from next screen
          cond.stop = tmp$LGD == TRUE | tmp$HGD == TRUE | tmp$IMC == TRUE | tmp$EAC == TRUE
          # Malign and Benign only need true status update after first Screen
          IndicatorFirstScreen = FALSE
        }
      
        # time increment s is defined in shared.R
        # update running time from BE onset to MaxAge
        sinceBEonset = runT + firstScreen-BEonset
        cond1 = Benign$s <= sinceBEonset  # present or coming up
        cond2 = Malign$s <= sinceBEonset  # present or coming up

        ################################################################################              
        # go through Benign/Malign data frames using apply function on all clones of age < runT
        Benign[cond1,] = as.data.frame(t(apply(Benign[cond1,],1,GENbenign,s=s,p=p0,aX=alphaXsi)))
        # Benign = Benign[Benign$Y > 0,] # only keep benign clones with Y>0 in the dataframe
        Malign[cond2,] = as.data.frame(t(apply(Malign[cond2,],1,GENmalign,s=s,p=p0,aX=alphaXsi)))
        ################################################################################              
        
        # advance by time s
        runT  = runT+s
      }
    }
    # concatenate all individual status reports including BE negative screens 
    Status = rbind(Status,datf)
    if(n%%1000==0) cat('done: ',n,'\n')
  }
)
cat('system time info',info,'\n\n')

################################################################################
# post simulation processing and stats
################################################################################


# fraction of BEonsets not occuring in lifetime (Tdeath)
popfracBE = BEprob[length(BEprob)]


# tmp=tmp[tmp$age==BEscreenings[length(BEscreenings)],]
tmp=Status[Status$age==max(BEscreenings),]; N.AR = nrow(tmp)
cat('MaxAge: ',MaxAge,'  BEScreen: ',unique(tmp$age),'\n')
cat('LGD: ',sum(tmp$LGD)/N.AR,'\n')
cat('HGD: ',sum(tmp$HGD)/N.AR,'\n')
cat('IMC: ',sum(tmp$IMC)/N.AR,'\n')
cat('InSitu: ',sum(tmp$InSitu)/N.AR,'\n')
cat('EAC: ',sum(tmp$EAC)/N.AR,'\n')
cat('preScrEAC: ',sum(tmp$EACpreScreen)/N.AR,'\n')
cat('postScrEAC: ',sum(tmp$EACpostScreen)/N.AR,'\n')

