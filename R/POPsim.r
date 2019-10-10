POPsim <- function(survey1,survey2,survey3,
                   paa.catch.female,paa.catch.male,
                   n.trips.paa.catch,
                   paa.survey1.female,paa.survey1.male,
                   n.trips.paa.survey1,
                   catch,paa.mature,weight.female,weight.male,
                   misc.fixed.param=NULL,
                   theta.ini=NULL,
                   lkhd.paa="normal",
                   var.paa.add=TRUE,
                   enable.priors=TRUE){
  #/////////////////////////////////////////////////////////////////////////////
  #### Documentation ####
  #/////////////////////////////////////////////////////////////////////////////
  #' @title Simulate data according to POP model.
  #' @description
  #' Uses fixed covariates and parameters to generate random effects and
  #' response variables. Outputs a list of class "POPobj" that can be directly
  #' supplied to \code{\link{POPfit}}.
  #' @param Ct Vector of length TC; total catch, fixed covariate.
  #' @details To be used with the POP cpp template. Compatible with v0.4.
  #' @return A list of class POPobj properly formatted to be fed to POPfit,
  #' with the following objects:
  #' \itemize{
  #'   \item datalist: a list containing all DATA inputs for the TMB POP
  #'   template;
  #'   \item parlist: a list containing all PARAMETER inputs for the TMB POP
  #'   template;
  #' }
  # TODO: finish roxygen2 doc
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Functions ####
  #/////////////////////////////////////////////////////////////////////////////
  
  # bound11 <- function(x){(1-exp(-x))/1+exp(-x)}
  invlogitSelectivity <- function(a,mu,ups){
    tmp <- exp((a-(mu-ups/2))/ups*10)
    return(tmp/(1+tmp))
  } # steepness=10 mimics well selectivitiy curve in (F.7) and (F.8) with ups<5
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Setup ####
  #/////////////////////////////////////////////////////////////////////////////
  
  if (is.null(theta.ini)){ # then default starting values
    # fixed param, not in user-supplied data
    R0 <- 5000 # 4000 in Andy's priors
    M1 <- 0.07
    M2 <- 0.07
    muC <- 10.5
    deltaC <- 0
    upsilonC <- 5 # different meaning than in Andy's code, smaller than his
    muS1 <- 13.3
    deltaS1 <- 0.22
    upsilonS1 <- 10  # different meaning than in Andy's code, smaller than his
    h <- 0.674
    qS1 <- 1
    qS2 <- 1
    qS3 <- 1
  }
  
  yearsvec <- catch[,1]
  length.theta <- 13 # better way?
  
  if (is.null(misc.fixed.param)){ # then default misc param values
    muS2 <- 13.3
    deltaS2 <- 0.22
    upsilonS2 <- 10
    muS3 <- 13.3
    deltaS3 <- 0.22
    upsilonS3 <- 10
    sigmaR <- 0.05 # 0.9 in Andy's code
  } else {
    muS2 <- as.numeric(misc.fixed.param['muS2'])
    deltaS2 <- as.numeric(misc.fixed.param['deltaS2'])
    upsilonS2 <- as.numeric(misc.fixed.param['upsilonS2'])
    muS3 <- as.numeric(misc.fixed.param['muS3'])
    deltaS3 <- as.numeric(misc.fixed.param['deltaS3'])
    upsilonS3 <- as.numeric(misc.fixed.param['upsilonS3'])
    sigmaR <- as.numeric(misc.fixed.param['sigmaR'])
  }
  
  A <- length(weight.female)
  TC <- dim(catch)[[1]]
  TS1 <- dim(survey1)[[1]]
  TS2 <- dim(survey2)[[1]]
  TS3 <- dim(survey3)[[1]]
  UC <- dim(paa.catch.female)[[1]]
  US1 <- dim(paa.survey1.female)[[1]]
  
  # R indices, converted for C++ in Outputs
  tS1 <- which(catch[,1]%in%survey1[,1])
  tS2 <- which(catch[,1]%in%survey2[,1])
  tS3 <- which(catch[,1]%in%survey3[,1])
  tUC <- which(catch[,1]%in%paa.catch.female[,1])
  tUS1 <- which(catch[,1]%in%paa.survey1.female[,1])
  
  Ct <- catch[,2]
  ntC <- n.trips.paa.catch[,2]
  ntS1 <- n.trips.paa.survey1[,2]
  wa1 <- weight.female
  wa2 <- weight.male
  ma <- paa.mature
  
  kappaS1 <- survey1[,3]
  kappaS2 <- survey2[,3]
  kappaS3 <- survey3[,3]
  
  if (lkhd.paa=="normal"){
    lkhdpropatage <- 1L
  } else if (lkhd.paa=="binomial"){
    lkhdpropatage <- 2L
  } else {stop('lkhd.paa can only be "normal" or "binomial".')}
  
  varweight <- as.integer(var.paa.add)
  enablepriors <- as.integer(enable.priors)
  
  G <- 4 # number of fleets # TODO: config, allow for more than 4?
  indC <- 1 # index for catch data
  indS1 <- 2 # g index for S1
  indS2 <- 3 # g index for S2
  indS3 <- 4 # g index for S3
  
  Nat1 <- matrix(NA_real_,A,TC) # abundance females
  Nat2 <- Nat1 # abundance males
  sag1 <- matrix(NA_real_,A,G) # selectivity females
  sag2 <- sag1 # selectivity males
  Bt <- rep(NA_real_,TC) # 1<=t<=T, B_0 is separate
  Rt <- Bt # number of recruits, only strict randeff
  Vt <- Bt # vulnerable biomass
  ut <- Bt # alternative to fishing mortality
  uat1 <- Nat1
  uat2 <- Nat1
  
  patC1 <- matrix(NA_real_,A,UC) # prop-at-age catch females
  patC2 <- patC1 # prop-at-age catch males
  patS11 <- matrix(NA_real_,A,US1) # prop-at-age S1 females
  patS12 <- patS11 # prop-at-age S1 males
  
  expnegM1 <- exp(-M1)
  expnegM2 <- exp(-M2)
  sqrtexpnegM1 <- sqrt(expnegM1) # exp(-M1*0.5)
  sqrtexpnegM2 <- sqrt(expnegM2) # exp(-M2*0.5)
  
  
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Generate iid process errors ####
  #/////////////////////////////////////////////////////////////////////////////
  
  errRt <- rlnorm(n=TC,meanlog=-sigmaR^2/2,sdlog=sigmaR) # proc error Rt
  errS1t <- rlnorm(n=TS1,meanlog=-kappaS1^2/2,sdlog=kappaS1) # obs error S1t
  errS2t <- rlnorm(n=TS2,meanlog=-kappaS2^2/2,sdlog=kappaS2) # obs error S2t
  errS3t <- rlnorm(n=TS3,meanlog=-kappaS3^2/2,sdlog=kappaS3) # obs error S3t
  
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### States dynamics ####
  #/////////////////////////////////////////////////////////////////////////////
  
  ### Selectivities, all fleets
  
  for (a in 1:A){
    sag1[a,indC] <- invlogitSelectivity(a,muC,upsilonC) # catch females (F.7)
    sag2[a,indC] <- invlogitSelectivity(a,muC+deltaC,upsilonC) # catch males (F.8)
    sag1[a,indS1] <- invlogitSelectivity(a,muS1,upsilonS1) # S1 females (F.7)
    sag2[a,indS1] <- invlogitSelectivity(a,muS1+deltaS1,upsilonS1) # S1 males (F.8)
    sag1[a,indS2] <- invlogitSelectivity(a,muS2,upsilonS2) # S2 females (F.7)
    sag2[a,indS2] <- invlogitSelectivity(a,muS2+deltaS2,upsilonS2) # S2 males (F.8)
    sag1[a,indS3] <- invlogitSelectivity(a,muS3,upsilonS3) # S3 females (F.7)
    sag2[a,indS3] <- invlogitSelectivity(a,muS3+deltaS3,upsilonS3) # S3 males (F.8)
  }
  
  ### init Nat, Bt, Rt, Vt, ut and uat
  
  # Nat, t=1, s=1,2
  for (a in 1:(A-1)){
    Nat1[a,1] <- 0.5*R0*exp(-M1*a) # (F.4)
    Nat2[a,1] <- 0.5*R0*exp(-M2*a) # (F.4)
  }
  Nat1[A,1] <- 0.5*R0*exp(-M1*(A-1))/(1-expnegM1) # (F.5)
  Nat2[A,1] <- 0.5*R0*exp(-M2*(A-1))/(1-expnegM2) # (F.5)
  
  # Bt, t=0,1
  Bt[1] <- as.numeric((wa1*ma)%*%Nat1[,1]) # t=1 (F.6)
  B0 <- Bt[1] # t=0 (F.6)
  
  # Rt, t=1
  meanR0 <- R0 # (F.10)
  Rt[1] <- meanR0*errRt[1] # (F.17)
  
  # Vt and ut, requires Ct[1], t=1
  sumwsN1 <- (wa1*sag1[,indC])%*%Nat1[,1] # (F.11)
  sumwsN2 <- (wa2*sag2[,indC])%*%Nat2[,1] # (F.11)
  Vt[1] <- sqrtexpnegM1*sumwsN1+sqrtexpnegM2*sumwsN2 # (F.11)
  ut[1] <- Ct[1]/Vt[1] # Ct must be a fixed covariate (F.12)
  
  # uat, t=1
  uat1[,1] <- sag1[,indC]*ut[1] # (F.13)
  uat2[,1] <- sag2[,indC]*ut[1] # (F.13)
  
  ### dynamics for Rt, Nat, Bt, Vt, ut and uat, 2<=t<=T
  
  for (t in 2:TC){
    # Rt, based on Bt[t-1]
    meanRt <- 4*h*R0*Bt[t-1]/((1-h)*B0+(5*h-1)*Bt[t-1]) # (F.10)
    Rt[t] <- meanRt*errRt[t] # (F.17)
    
    # Nat, based on Rt[t], uat[,t-1] and Nat[,t-1]
    Nat1[1,t] <- 0.5*Rt[t] # (F.1)
    Nat2[1,t] <- 0.5*Rt[t] # (F.1)
    for (a in 2:(A-1)){
      Nat1[a,t] <- expnegM1*(1-uat1[a-1,t-1])*Nat1[a-1,t-1] # (F.2)
      Nat2[a,t] <- expnegM2*(1-uat2[a-1,t-1])*Nat2[a-1,t-1] # (F.2)
    }
    Nat1[A,t] <- expnegM1*(1-uat1[A-1,t-1])*Nat1[A-1,t-1]+
      expnegM1*(1-uat1[A,t-1])*Nat1[A,t-1] # (F.3)
    Nat2[A,t] <- expnegM2*(1-uat2[A-1,t-1])*Nat2[A-1,t-1]+
      expnegM2*(1-uat2[A,t-1])*Nat2[A,t-1] # (F.3)
    
    # Bt, based on Nat[,t]
    Bt[t] <- as.numeric((wa1*ma)%*%Nat1[,t]) # (F.9)
    
    # Vt and ut, based on Nat[,t] an Ct[t]
    sumwsN1 <- (wa1*sag1[,indC])%*%Nat1[,t] # (F.11)
    sumwsN2 <- (wa2*sag2[,indC])%*%Nat2[,t] # (F.11)
    Vt[t] <- sqrtexpnegM1*sumwsN1+sqrtexpnegM2*sumwsN2 # (F.11)
    ut[t] <- Ct[t]/Vt[t] # Ct must be a fixed covariate (F.12)
    
    # uat, based on ut
    uat1[,t] <- sag1[,indC]*ut[t] # (F.13)
    uat2[,t] <- sag2[,indC]*ut[t] # (F.13)
  }
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Observation equations ####
  #/////////////////////////////////////////////////////////////////////////////
  
  # survey index S1, t in tS1
  sumuwsN1 <- as.numeric(
    sag1[,indS1]%*%((1-0.5*uat1[,tS1])*Nat1[,tS1]*matrix(rep(wa1,TS1),A,TS1))) # (F.14)
  sumuwsN2 <- as.numeric(
    sag2[,indS1]%*%((1-0.5*uat2[,tS1])*Nat2[,tS1]*matrix(rep(wa2,TS1),A,TS1))) # (F.14)
  meanS1t <- qS1*(sqrtexpnegM1*sumuwsN1+sqrtexpnegM2*sumuwsN2) # (F.14)
  S1t <- meanS1t*errS1t # (F.20)
  
  # survey index S2, t in tS2
  sumuwsN1 <- as.numeric(
    sag1[,indS2]%*%((1-0.5*uat1[,tS2])*Nat1[,tS2]*matrix(rep(wa1,TS2),A,TS2))) # (F.14)
  sumuwsN2 <- as.numeric(
    sag2[,indS2]%*%((1-0.5*uat2[,tS2])*Nat2[,tS2]*matrix(rep(wa2,TS2),A,TS2))) # (F.14)
  meanS2t <- qS2*(sqrtexpnegM1*sumuwsN1+sqrtexpnegM2*sumuwsN2) # (F.14)
  S2t <- meanS2t*errS2t # (F.20)
  
  # survey index S3, t in tS3
  sumuwsN1 <- as.numeric(
    sag1[,indS3]%*%((1-0.5*uat1[,tS3])*Nat1[,tS3]*matrix(rep(wa1,TS3),A,TS3))) # (F.14)
  sumuwsN2 <- as.numeric(
    sag2[,indS3]%*%((1-0.5*uat2[,tS3])*Nat2[,tS3]*matrix(rep(wa2,TS3),A,TS3))) # (F.14)
  meanS3t <- qS3*(sqrtexpnegM1*sumuwsN1+sqrtexpnegM2*sumuwsN2) # (F.14)
  S3t <- meanS3t*errS3t # (F.20)
  
  # prop at age C, t in tUC
  if (lkhdpropatage==1){ # Gaussian lkhd
    if (varweight==0){ # no variance inflation
      varcst <- 0 # no effect
    } else if (varweight==1){ # inflate prop-at-age variance by adding varcst
      varcst <- 1/(10*A) # F.5.3, Stanley et al. (2009) # (F.19)
    } else {
      stop('varweight can only be "0" (disabled) or "1" (variance inflation).');
    }
    for (t in 1:UC){
      ind <- tUC[t] # scalar ind
      usN1 <- (1-0.5*uat1[,ind])*sag1[,indC]*Nat1[,ind] # (F.15)
      usN2 <- (1-0.5*uat2[,ind])*sag2[,indC]*Nat2[,ind] # (F.15)
      sumusN1 <- sum(usN1) # (F.15)
      sumusN2 <- sum(usN2) # (F.15)
      sumusN <- (sqrtexpnegM1*sumusN1+sqrtexpnegM2*sumusN2) # (F.15)
      meanpat1 <- sqrtexpnegM1*usN1/sumusN # (F.15)
      meanpat2 <- sqrtexpnegM2*usN2/sumusN # (F.15)
      # meanpat instead of observed patC[,t], cannot follow exact Coleraine
      sdpat1 <- sqrt((meanpat1*(1-meanpat1)+varcst)/ntC[t]) # (F.19)
      sdpat2 <- sqrt((meanpat2*(1-meanpat2)+varcst)/ntC[t]) # (F.19)
      patC1[,t] <- pmax(rnorm(n=A,mean=meanpat1,sd=sdpat1),0) # (F.19)
      patC2[,t] <- pmax(rnorm(n=A,mean=meanpat2,sd=sdpat2),0) # (F.19)
    }
  } else if (lkhdpropatage==2){ # binomial lkhd
    for (t in 1:UC){
      ind <- tUC[t] # scalar ind
      usN1 <- (1-0.5*uat1[,ind])*sag1[,indC]*Nat1[,ind] # (F.15)
      usN2 <- (1-0.5*uat2[,ind])*sag2[,indC]*Nat2[,ind] # (F.15)
      sumusN1 <- sum(usN1) # (F.15)
      sumusN2 <- sum(usN2) # (F.15)
      sumusN <- (sqrtexpnegM1*sumusN1+sqrtexpnegM2*sumusN2) # (F.15)
      meanpat1 <- sqrtexpnegM1*usN1/sumusN # (F.15)
      meanpat2 <- sqrtexpnegM2*usN2/sumusN # (F.15)
      patC1[,t] <- rbinom(n=A,prob=meanpat1,size=ntC[t])/ntC[t]
      patC2[,t] <- rbinom(n=A,prob=meanpat2,size=ntC[t])/ntC[t]
    }
  } else {stop('lkhdpropatage can only be "1" (Gaussian) or "2" (binomial).')}
  
  
  # prop at age S1, t in tUS1
  if (lkhdpropatage==1){ # Gaussian lkhd
    if (varweight==0){ # no variance inflation
      varcst <- 0 # no effect
    } else if (varweight==1){ # inflate prop-at-age variance by adding varcst
      varcst <- 1/(10*A) # F.5.3, Stanley et al. (2009) # (F.19)
    } else {
      stop('varweight can only be "0" (disabled) or "1" (variance inflation).');
    }
    for (t in 1:US1){
      ind <- tUS1[t] # scalar ind
      usN1 <- (1-0.5*uat1[,ind])*sag1[,indS1]*Nat1[,ind] # (F.15)
      usN2 <- (1-0.5*uat2[,ind])*sag2[,indS1]*Nat2[,ind] # (F.15)
      sumusN1 <- sum(usN1) # (F.15)
      sumusN2 <- sum(usN2) # (F.15)
      sumusN <- (sqrtexpnegM1*sumusN1+sqrtexpnegM2*sumusN2) # (F.15)
      meanpat1 <- sqrtexpnegM1*usN1/sumusN # (F.15)
      meanpat2 <- sqrtexpnegM2*usN2/sumusN # (F.15)
      # meanpat instead of observed patS1[,t], cannot follow exact Coleraine
      sdpat1 <- sqrt((meanpat1*(1-meanpat1)+varcst)/ntS1[t]) # (F.19)
      sdpat2 <- sqrt((meanpat2*(1-meanpat2)+varcst)/ntS1[t]) # (F.19)
      patS11[,t] <- pmax(rnorm(n=A,mean=meanpat1,sd=sdpat1),0) # (F.19)
      patS12[,t] <- pmax(rnorm(n=A,mean=meanpat2,sd=sdpat2),0) # (F.19)
    }
  } else if (lkhdpropatage==2){ # binomial lkhd
    for (t in 1:US1){
      ind <- tUS1[t] # scalar ind
      usN1 <- (1-0.5*uat1[,ind])*sag1[,indS1]*Nat1[,ind] # (F.15)
      usN2 <- (1-0.5*uat2[,ind])*sag2[,indS1]*Nat2[,ind] # (F.15)
      sumusN1 <- sum(usN1) # (F.15)
      sumusN2 <- sum(usN2) # (F.15)
      sumusN <- (sqrtexpnegM1*sumusN1+sqrtexpnegM2*sumusN2) # (F.15)
      meanpat1 <- sqrtexpnegM1*usN1/sumusN # (F.15)
      meanpat2 <- sqrtexpnegM2*usN2/sumusN # (F.15)
      patS11[,t] <- rbinom(n=A,prob=meanpat1,size=ntS1[t])/ntS1[t]
      patS12[,t] <- rbinom(n=A,prob=meanpat2,size=ntS1[t])/ntS1[t]
    }
  } else {stop('lkhdpropatage can only be "1" (Gaussian) or "2" (binomial).')}
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Outputs ####
  #/////////////////////////////////////////////////////////////////////////////
  
  # convert R indices into C++ indices
  tS1 <- as.integer(tS1-1)
  tS2 <- as.integer(tS2-1)
  tS3 <- as.integer(tS3-1)
  tUC <- as.integer(tUC-1)
  tUS1 <- as.integer(tUS1-1)
  
  true.theta <- c(R0,M1,M2,muC,deltaC,upsilonC,
                  muS1,deltaS1,upsilonS1,h,qS1,qS2,qS3)
  
  parlist <- list('logR0'=log(R0),
                  'logM1'=log(M1),'logM2'=log(M2),
                  'logmuC'=log(muC),'deltaC'=deltaC,'logupsilonC'=log(upsilonC),
                  'logmuS1'=log(muS1),'deltaS1'=deltaS1,
                  'logupsilonS1'=log(upsilonS1),
                  'logh'=log(h),'logqS1'=log(qS1),'logqS2'=log(qS2),
                  'logqS3'=log(qS3),
                  'logRt'=rep(log(R0),TC))
  
  datalist <- list('S1t'=S1t,'S2t'=S2t,'S3t'=S3t, # response
                   'patC1'=patC1,'patC2'=patC2, # response
                   'patS11'=patS11,'patS12'=patS12, # response
                   'Ct'=Ct,'ntC'=ntC,'ntS1'=ntS1, # fixed covariate
                   'wa1'=wa1,'wa2'=wa2,'ma'=ma, # fixed covariate
                   'tS1'=tS1,'tS2'=tS2,'tS3'=tS3,'tUC'=tUC,'tUS1'=tUS1, # indices
                   'kappaS1'=kappaS1,'kappaS2'=kappaS2,'kappaS3'=kappaS3, # fixed
                   'muS2'=muS2,'deltaS2'=deltaS2,'upsilonS2'=upsilonS2, # fixed
                   'muS3'=muS3 ,'deltaS3'=deltaS3,'upsilonS3'=upsilonS3, # fixed
                   'sigmaR'=sigmaR, # fixed
                   'lkhdpropatage'=lkhdpropatage,'varweight'=varweight, # options
                   'enablepriors'=as.integer(enable.priors) # options
  )
  
  res <- list('parlist'=parlist,'datalist'=datalist,
              'true.theta'=true.theta,
              'true.R'=Rt,'true.B'=Bt,'true.V'=Vt,'true.u'=ut,
              'true.uaa.female'=uat1,'true.uaa.male'=uat2,
              'true.select.female'=sag1,'true.select.male'=sag2,
              'true.abund.female'=Nat1,'true.abund.male'=Nat2,
              'length.theta'=length.theta,'years'=yearsvec,
              'A'=A,'TC'=TC,'TS1'=TS1,'TS2'=TS2,'TS3'=TS3,'UC'=UC,'US1'=US1)
  class(res) <- 'POPobj' # useful?
  return(res)
  
}
# END POPsim
