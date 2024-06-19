# simulate continuously distributed IPD

cont.ipd = function(N=10000,Nt=10,MD=0.1,SD=1,het=0,phi=0,tsize="random"){

  # trial sizes
  if(tsize=="random"){
    trialsize = rep(0,Nt)
    Nrem = N
    for (i in 1:(Nt-1)){
      if (Nrem<50)
        trialsize[i] = floor(Nrem/2)
      else
        trialsize[i] = floor(runif(1,25,(0.5*Nrem)))
      Nrem = Nrem-trialsize[i]
    }
    trialsize[Nt] = Nrem
    N = sum(trialsize)
  } else{
    trialsize = rep(floor(N/Nt),Nt)
    N = sum(trialsize)
  }

  # hetgty
  ran.eff = rnorm(Nt,0,sqrt(het))
  MD.trial = MD + ran.eff
  MD.pers = rep(MD.trial,trialsize)

  # baseline
  phi.i = rnorm(Nt,phi,0.1)
  phi.pers = rep(phi.i,trialsize)

  # treat and trial
  treat = rep(c(0,1),length.out=N)
  trial = rep(1:Nt,trialsize)

  # IPD effects
  response.mean = phi.pers + MD.pers * treat
  SD.trial = rnorm(Nt,SD,0.05)
  SD.pers = rep(SD.trial,trialsize)
  response = rnorm(N,response.mean,SD.pers)

  IPD = data.frame(trial=letters[trial],trialn=trial,treat=ifelse(treat==0,"Control","Experimental"),treatn=treat,response=response)

  # convert IPD to AD
  n.gp = t(matrix(tapply(response,factor(trial):factor(treat),length),2,Nt))
  mean.gp = t(matrix(tapply(response,factor(trial):factor(treat),mean),2,Nt))
  sd.gp = t(matrix(tapply(response,factor(trial):factor(treat),sd),2,Nt))

  AD = data.frame(trial=letters[1:Nt],trialn=1:Nt,N.exp=n.gp[,2],Mean.exp=mean.gp[,2],SD.exp=sd.gp[,2],N.cont=n.gp[,1],Mean.cont=mean.gp[,1],SD.cont=sd.gp[,1])

  out = list(IPD=IPD,AD=AD)
  return(out)

}
