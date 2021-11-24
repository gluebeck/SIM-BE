BEstatusByBiopsy = function(ID=0,age=age,BEonset=BEonset,Malign=Malign) {

  # (LGD,HGD,IMC,InSitu) evaluator via quadrant biopsies

  Probs.LGD = detectProb(Malign$Y[Malign$LGD==T],scaleDetection)
  Probs.HGD = detectProb(Malign$Y[Malign$HGD==T],scaleDetection)
  Probs.IMC = detectProb(Malign$Z,scaleDetection)
  Probs.InSitu = as.numeric(Malign$Z > Threshold.InSitu)  # deterministic 0,1
  Probs.EAC =  1-(1-rho.EAC)^Malign$Z
  
  BEstate.BE = BEonset <= age
  BEstate.lgd = any(rbinom(length(Probs.LGD),1,Probs.LGD) == 1) # stochastic 0,1 
  BEstate.hgd = any(rbinom(length(Probs.HGD),1,Probs.HGD) == 1) # stochastic 0,1
  # if(BEstate.hgd==TRUE) BEstate.lgd=FALSE # use highest dysplasia grade for diagnsostic status of patient at this exam
  
  if(BEstate.lgd == TRUE & HGD == FALSE)  LGD=TRUE  # permanent until HGD is found
  if(BEstate.hgd == TRUE) {LGD=FALSE; HGD=TRUE} # this is permanent

  BEstate.EAC = any(rbinom(length(Probs.EAC),1,Probs.EAC) == 1)
  BEstate.IMC = BEstate.InSitu = FALSE

  # BEstate.EAC = any(runif(length(Malign$Z)) < Probs.EAC)
  dum = rbinom(length(Probs.IMC),1,Probs.IMC) # if == 1 then biopsy shows clone that can be considered at least IMC or more advanced
  if(sum(dum)>0) {
    BEstate.IMC =    all(Malign$Z[dum == 1] < Threshold.InSitu)  # all malignant biopsies are sub-InSitu threshold
    BEstate.InSitu = any(Malign$Z[dum == 1] > Threshold.InSitu) & !BEstate.EAC  # but not yet clinical cancer
  }
  # by fpT and tlag
  TclinicalCancer = BEonset + min(Malign$fpT)+tlag
  BEstate.EACpreScreen =   any(TclinicalCancer <= age)
  BEstate.EACpostScreen =  any(TclinicalCancer >  age & TclinicalCancer < Tdeath)
  
  # NDBE
  BEstate.NDBE = !LGD & !HGD & !BEstate.IMC & !BEstate.InSitu & !BEstate.EACpreScreen
  
  return(list(ID=ID,
              age=age,
              ageClinicalEAC=BEonset+min(Malign$fpT)+tlag,
              BEonset = BEonset,
              BE = BEstate.BE,
              NDBE=BEstate.NDBE,
              LGD=LGD,lgd=BEstate.lgd, #lower case = current obs
              HGD=HGD,hgd=BEstate.hgd, #lower case = current obs
              IMC=BEstate.IMC,
              InSitu=BEstate.InSitu,
              EAC=BEstate.EAC,
              EACpreScreen =BEstate.EACpreScreen,
              EACpostScreen=BEstate.EACpostScreen))
}
