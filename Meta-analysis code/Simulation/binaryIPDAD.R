# simulate binary meta-analysis data: IPD and AD
# allow for ecological bias and "publication" bias

binarydata <- 

function(N=10000,Nt=10,theta=0.1,mu=0,gamma=0,delta=0,het=0,phi=0,m=0,w=10,bias=0){

# generate IPD

  # trial sizes
  trialsize <- rep(0,Nt)
  Nrem <- N 
  for (i in 1:(Nt-1)){
    if (Nrem<50)
      trialsize[i] <- floor(Nrem/2)
    else
      trialsize[i] <- floor(runif(1,25,(0.5*Nrem)))
    Nrem <- Nrem-trialsize[i]
  }
  trialsize[Nt] <- Nrem
  N <- sum(trialsize)

  # small study effect (pub. bias)
  theta.trial <- rep(theta,Nt)
  theta.trial <- theta - (bias/1000)*(trialsize-max(trialsize))

  # hetgty
  ran.eff <- rnorm(Nt,0,het)
  theta.trial <- theta.trial + ran.eff
  theta.pers <- rep(theta.trial,trialsize)
  
  phi.i <- rnorm(Nt,phi,0.1)
  phi.pers <- rep(phi.i,trialsize)

  # treat and trial
  treat <- rep(c(0,1),length.out=N)
  treathalf <- rep(c(-0.5,0.5),length.out=N)
  trial <- rep(1:Nt,trialsize)

  # covariate
  m1 <- m
  if(m==0)
    m1 <- rep(0,Nt)
  covar.cent <- rep(m1,trialsize)
  covar <- runif(N,covar.cent-w,covar.cent+w)
  covar <- covar - mean(covar)
  meancovar <- rep(0,Nt)
  for (i in 1:Nt)
    meancovar[i] <- mean(covar[trial==i])
  meancovar <- rep(meancovar,trialsize)
  covar.cent <- covar - meancovar

  #LOR and events
  p.e <- expit( phi.pers + theta.pers*treat + mu*covar + gamma*covar*treat + delta*treat*meancovar )
  

  event <- rbinom(N,1,p.e)

  IPD <- data.frame(trial=letters[trial],trialn=trial,treat=ifelse(treat==0,"Control","Experimental"),treatn=treat,treathalf=treathalf,event=ifelse(event==0,"Event","No event"),eventn=event,covar=round(covar,2),meancovar=round(meancovar,2),covar.cent=round(covar.cent,2))

# generate AD from IPD

  n.t <- n.c <- e.t <- e.c <- rep(0,Nt)
  LOR <- se.LOR <- meancovar.ad <- rep(0,Nt)

  for (i in 1:Nt){

    n.t[i] <- length(treat[treat==1 & trial==i])
    n.c[i] <-  trialsize[i] - n.t[i]
    
    e.t[i] <- length(treat[event==1 & treat==1 & trial==i])
    e.c[i] <- length(treat[event==1 & treat==0 & trial==i])

    meancovar.ad[i] <- mean(subset(covar,trial==i))
  }

  LOR <- log( ( e.t * (n.c-e.c) ) / ( e.c * (n.t-e.t )) )
  se.LOR <- 1/e.t + 1/(n.t-e.t) + 1/e.c + 1/(n.c-e.c)

  AD <- data.frame(trial=1:Nt,trialsize,n.t,n.c,e.t,e.c,LOR,se.LOR,true.LOR=theta.trial,meancovar=meancovar.ad)

  out <- list(IPD=IPD,AD=AD)
  return(out)

}


logit <-
  function(x){
    return( log(x / (1-x)) )
  }

expit <-
  function(x){
    return( exp(x)/(1+exp(x)) )
  }

