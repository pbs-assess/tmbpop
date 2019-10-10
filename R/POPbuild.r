POPbuild <- function(survey1,survey2,survey3,
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
  #' @title Formats data files for POPfit.
  #' @description
  #' to do later
  #' @param catch Data frame (or matrix) of commercial trawl data, dimensions
  #' are TC rows by 2 columns. The first column is the year, while the second
  #' column is the observed catch.
  #' @details To be used with the POP cpp template.
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
  } else {
    R0 <- as.numeric(theta.ini['R0'])
    M1 <- as.numeric(theta.ini['M1'])
    M2 <- as.numeric(theta.ini['M2'])
    muC <- as.numeric(theta.ini['muC'])
    deltaC <- as.numeric(theta.ini['deltaC'])
    upsilonC <- as.numeric(theta.ini['upsilonC'])
    muS1 <- as.numeric(theta.ini['muS1'])
    deltaS1 <- as.numeric(theta.ini['deltaS1'])
    upsilonS1 <- as.numeric(theta.ini['upsilonS1'])
    h <- as.numeric(theta.ini['h'])
    qS1 <- as.numeric(theta.ini['qS1'])
    qS2 <- as.numeric(theta.ini['qS2'])
    qS3 <- as.numeric(theta.ini['qS3'])
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
  
  tS1 <- which(catch[,1]%in%survey1[,1])-1L # C++ indices
  tS2 <- which(catch[,1]%in%survey2[,1])-1L # C++ indices
  tS3 <- which(catch[,1]%in%survey3[,1])-1L # C++ indices
  tUC <- which(catch[,1]%in%paa.catch.female[,1])-1L # C++ indices
  tUS1 <- which(catch[,1]%in%paa.survey1.female[,1])-1L # C++ indices
  
  kappaS1 <- survey1[,3]
  kappaS2 <- survey2[,3]
  kappaS3 <- survey3[,3]
  
  muS2 <- as.numeric(MiscFixedParam['muS2'])
  deltaS2 <- as.numeric(MiscFixedParam['deltaS2'])
  upsilonS2 <- as.numeric(MiscFixedParam['upsilonS2'])
  muS3 <- as.numeric(MiscFixedParam['muS3'])
  deltaS3 <- as.numeric(MiscFixedParam['deltaS3'])
  upsilonS3 <- as.numeric(MiscFixedParam['upsilonS3'])
  sigmaR <- as.numeric(MiscFixedParam['sigmaR'])
  
  if (lkhd.paa=="normal"){
    lkhdpropatage <- 1L
  } else if (lkhd.paa=="binomial"){
    lkhdpropatage <- 2L
  } else {stop('lkhd.paa can only be "normal" or "binomial".')}
  
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Outputs ####
  #/////////////////////////////////////////////////////////////////////////////
  
  parlist <- list('logR0'=log(R0),
                  'logM1'=log(M1),'logM2'=log(M2),
                  'logmuC'=log(muC),'deltaC'=deltaC,'logupsilonC'=upsilonC,
                  'logmuS1'=log(muS1),'deltaS1'=deltaS1,
                  'logupsilonS1'=log(upsilonS1),
                  'logh'=log(h),'logqS1'=log(qS1),'logqS2'=log(qS2),
                  'logqS3'=log(qS3),
                  'logRt'=rep(log(R0),TC))
  
  datalist <- list('Ct'=catch[,2],
                   'S1t'=survey1[,2],'S2t'=survey2[,2],'S3t'=survey3[,2],
                   'patC1'=t(as.matrix(unname(paa.catch.female[,-1]))), # AxUC
                   'patC2'=t(as.matrix(unname(paa.catch.male[,-1]))), # AxUC
                   'ntC'=n.trips.paa.catch[,-1],
                   'patS11'=t(as.matrix(unname(paa.survey1.female[,-1]))), # AxUS1
                   'patS12'=t(as.matrix(unname(paa.survey1.female[,-1]))), # AxUS1
                   'ntS1'=n.trips.paa.survey1[,-1],
                   'tS1'=tS1,'tS2'=tS2,'tS3'=tS3,'tUC'=tUC,'tUS1'=tUS1,
                   'wa1'=weight.female,'wa2'=weight.male,'ma'=paa.mature,
                   'kappaS1'=kappaS1,'kappaS2'=kappaS2,'kappaS3'=kappaS3,
                   'muS2'=muS2,'deltaS2'=deltaS2,'upsilonS2'=upsilonS2,
                   'muS3'=muS3 ,'deltaS3'=deltaS3,'upsilonS3'=upsilonS3,
                   'sigmaR'=sigmaR,
                   'lkhdpropatage'=lkhdpropatage,
                   'varweight'=as.integer(var.paa.add),
                   'enablepriors'=as.integer(enable.priors))
  
  res <- list('parlist'=parlist,'datalist'=datalist,
              'length.theta'=length.theta,'years'=yearsvec,
              'A'=A,'TC'=TC,'TS1'=TS1,'TS2'=TS2,'TS3'=TS3,'UC'=UC,'US1'=US1)
  class(res) <- 'POPobj' # useful?
  return(res)
}
# END POPbuild
