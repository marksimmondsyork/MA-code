# simple MA simulators for use in testing MA methods

#################################################
# generic data

ma_simul_gen <- function(N.studies,true.eff,hetgty=0,Nlow=20,Nhigh=100,SD=1,method="fixed"){

  N <- floor(runif(N.studies, Nlow, Nhigh))
  SE <- ifelse(method=="fixed", SD/sqrt(N), rnorm(N.studies, SD/sqrt(N), 0.25*min(SD/sqrt(N))))
  true.y <- rnorm(N.studies, true.eff, sqrt(hetgty))
  Eff.est <- rnorm(N.studies, true.y, SE)
  out <- data.frame(N, Eff.est, SE, true.y)
  return(out)
}

#####################################################
# binary data

ma_simul_bin <- function(N.studies,true.LOR,hetgty=0,Nlow=20,Nhigh=100,PC=0.2){

  N <- floor(runif(N.studies, Nlow, Nhigh))
  N.arm <- floor(N/2)

  t.LOR <- rnorm(N.studies, true.LOR, sqrt(hetgty))
  K <- exp(t.LOR)*PC/(1-PC)
  PE <- K/(1+K)

  ev.E <- rbinom(N.studies, N.arm, PE)
  #ev.E <- ifelse(ev.E==0, 1, ev.E)
  ev.C <- rbinom(N.studies, N.arm, PC)
  #ev.C <- ifelse(ev.C==0, 1, ev.C)
  LOR <- log((ev.E/(N.arm-ev.E)) / (ev.C/(N.arm-ev.C)))
  SE <- sqrt(1/ev.E + 1/ev.C + 1/(N.arm-ev.E) + 1/(N.arm-ev.C))

  out <- data.frame(N.E=N.arm, ev.E, N.C=N.arm, ev.C, LOR, SE, t.LOR)
  return(out)
}

##########################################################
# continuous data

ma_simul_cont <- function(N.studies,true.MD,hetgty=0,Nlow=20,Nhigh=100,SD=1){

  N <- floor(runif(N.studies, Nlow, Nhigh))
  N.arm <- floor(N/2)

  t.MD <- rnorm(N.studies, true.MD, sqrt(hetgty))

  M.E <- rnorm(N.studies, t.MD, SD/N.arm)
  M.C <- rnorm(N.studies, 0, SD/N.arm)

  MD <- M.E - M.C
  SE <- 2 * SD / N.arm

  out <- data.frame(N.E=N.arm, M.E, SD.E=SD, N.C=N.arm, M.C, SD.C=SD, MD, SE, t.MD)
  return(out)
}


#######################################################
# cont. data with specified total sample size

ma_simul_cont_2 <- function(N,N.studies,true.MD,M.C.low=0,M.C.high=0,hetgty=0,SD=1,sSD=0){

  N.trial <- stick_tsize(N,N.studies)
  N.arm <- floor(N.trial/2)

  t.MD <- rnorm(N.studies, true.MD, hetgty)

  SD <- rnorm(N.studies,SD,sSD*SD)

  M.C <- runif(N.studies,M.C.low,M.C.high)
  M.E <- rnorm(N.studies, M.C + t.MD, SD/N.arm)

  MD <- M.E - M.C
  SE <- sqrt(2 * SD^2 / N.arm)

  out <- data.frame(N.E=N.arm, M.E, SD.E=SD, N.C=N.arm, M.C, SD.C=SD, MD, SE, t.MD)
  return(out)
}


#####################################################
# binary data with specified sample size

ma_simul_bin_2 <- function(N,N.studies,true.LOR,hetgty=0,PC=0.2,sPC=0){

  N.trial <- stick_tsize(N,N.studies)
  N.arm <- floor(N.trial/2)

  t.LOR <- rnorm(N.studies, true.LOR, hetgty)
  PC.t = runif(N.studies,PC-PC*sPC,PC+PC*sPC)
  K <- exp(t.LOR)*PC.t/(1-PC.t)
  PE.t <- K/(1+K)

  ev.E <- rbinom(N.studies, N.arm, PE.t)
  #ev.E <- ifelse(ev.E==0, 1, ev.E)
  ev.C <- rbinom(N.studies, N.arm, PC.t)
  #ev.C <- ifelse(ev.C==0, 1, ev.C)
  LOR <- log((ev.E/(N.arm-ev.E)) / (ev.C/(N.arm-ev.C)))
  SE <- sqrt(1/ev.E + 1/ev.C + 1/(N.arm-ev.E) + 1/(N.arm-ev.C))

  out <- data.frame(N.E=N.arm, ev.E, N.C=N.arm, ev.C, LOR, SE, t.LOR)
  return(out)
}


###############################################
# Simulate Normal-Normal data (optional small sample effect)

ma_normnorm <- function(Nt,Nmin,Nmax,theta,tau2,SD=1,trend=0){
  
  trialsize <- floor(runif(Nt, Nmin, Nmax))
  
  theta.true <- rnorm(Nt, theta, sqrt(tau2))
  trialsize.c <- trialsize - mean(trialsize)
  theta.i.bias <- theta.true - trend * trialsize.c
  
  se.theta <- SD / sqrt(trialsize)
  theta.est <- rnorm(Nt, theta.i.bias, se.theta)
  out <- data.frame(Trial=1:Nt, N=trialsize, theta=theta.est, se=se.theta)
  return(out)
}


#################################################################
# Simulate Binomial-Normal data (optional small sample effect)

ma_binnorm <- function(Nt,Nmin,Nmax,theta,tau2,p.c=0.5,trend=0){
  
  trialsize <- floor(runif(Nt, Nmin, Nmax))
  armsize  <- floor(trialsize/2)
  
  theta.true <- rnorm(Nt, theta, sqrt(tau2))
  trialsize.c <- trialsize - mean(trialsize)
  theta.i.bias <- theta.true - trend * (0.01*trialsize.c)
  
  ev.c <- rbinom(Nt,armsize,p.c)
  odds.c <- p.c /(1-p.c)
  odds.e <- exp(theta.i.bias) * odds.c
  p.e <- odds.e / (1+odds.e)
  ev.e <- rbinom(Nt,armsize,p.e)
  res <- metabin(ev.e,armsize,ev.c,armsize)
  out <- data.frame(Trial=1:Nt, N=2*armsize, ev.e, N.e=armsize, ev.c, 
                    N.c=armsize,theta=res$TE,se=res$seTE)
  return(out)
}


######################################
# stick-breaking trial size calculator

stick_tsize <- function(N,N.studies,alpha=NULL,beta=1){


  alpha <- ifelse(is.null(alpha),N.studies,alpha)
  betas <- rbeta(N.studies-1, beta, alpha)
  stick.to.right <- c(1, cumprod(1 - betas))[1:(N.studies-1)]
  weights <- stick.to.right * betas
  N.trial  <- floor(N*weights)
  N.trial <- c(N.trial,N-sum(N.trial))
  N.trial <- ifelse(N.trial<10,10,N.trial)
  return(N.trial)
}




