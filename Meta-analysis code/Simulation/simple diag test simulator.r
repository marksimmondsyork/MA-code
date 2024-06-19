source("C:/mark/R/gamma dist.r")
source("C:/mark/r/logit.r")

data.simulate =

  function(n.studies=20,N=1000,DR50=0.5,sd.ratio=1,FPR=0.05,hetgty=0,p=0.1,distn="normal"){
    
    # study size
    #N.s = floor(runif(n.studies,N/2,2*N))
    N.s = rep(N,n.studies)
    N.tot = sum(N.s)
    study = rep(1:n.studies,N.s)

    # No. aff and unaff
    affected = rbinom(N.tot,1,p)
    unaffected = 1 - affected
    n.aff = as.numeric( tapply(affected, factor(study), sum) )
    n.unaff = N.s - n.aff

    # parameter from DR50
    if(distn=="logis")      
      param = logit(DR50)-logit(0.5)
    if(distn=="exp")
      param =  log(DR50) / log(0.5)
    if(distn=="gamma")
      param = gamma.dr(DR50,1,1)
    if(distn=="normal")
      param = qnorm(DR50)-qnorm(0.5)
      
    # hetgty in parameters
    param.1 = rnorm(n.studies,param,hetgty)
    par.pers = rep(param.1,N.s)

    # test score and test positives
    if(distn=="logis")      
      score = ifelse(affected==0, rlogis(N.tot,0,1), rlogis(N.tot,par.pers/sd.ratio,sd.ratio))
    if(distn=="exp")
      score = ifelse(affected==0, rexp(N.tot,1), rexp(N.tot,par.pers))
    if(distn=="gamma")
      score = ifelse(affected==0, rgamma.m(N.tot,1,1), rgamma.m(N.tot,par.pers,sd.ratio))
    if(distn=="normal")
      score = ifelse(affected==0, rnorm(N.tot,0,1), rnorm(N.tot,par.pers/sd.ratio,sd.ratio))

    # generate cutoffs
    if (length(FPR)==1)
      FPR = rep(FPR,n.studies)
    if(distn=="logis")
      cutoff = rep(qlogis(1-FPR,0,1),N.s)
    if(distn=="exp")
      cutoff = rep(qexp(1-FPR,1),N.s)
    if(distn=="gamma")
      cutoff = rep(qgamma.m(1-FPR,1,1),N.s)
    if(distn=="normal")
      cutoff = rep(qnorm(1-FPR,0,1),N.s)

    # generate data
    positive = ifelse(score > cutoff, 1, 0)
    true.pos = as.numeric( tapply(positive*affected, factor(study), sum) )
    false.pos = as.numeric( tapply(positive*unaffected, factor(study), sum) )

    # FPR, DR and DOR
    FPR.d = false.pos / n.unaff
    se.FPR = sqrt(FPR.d*(1-FPR.d) / n.unaff)
    DR = true.pos / n.aff
    se.DR = sqrt(DR*(1-DR) / n.unaff)
    DOR = (true.pos*(n.unaff-false.pos)) / (false.pos*(n.aff-true.pos))
    se.logDOR = sqrt(1/true.pos + 1/(n.aff-true.pos) + 1/false.pos + 1/(n.unaff-false.pos))

    # true values
    fpr.set = seq(0.01,0.99,0.01)
    if(distn=="logis"){
      DR.t = 1 - plogis(qlogis(1-FPR,0,1),sd.ratio*param,sd.ratio)
      DR.ta = 1 - plogis(qlogis(1-fpr.set,0,1),sd.ratio*param,sd.ratio)}
    if(distn=="exp"){
      DR.t = 1 - pexp(qexp(1-FPR,1),param)
      DR.ta = 1 - pexp(qexp(1-fpr.set,1),param)}
    if(distn=="gamma"){
      DR.t = 1 - pgamma.m(qgamma.m(1-FPR,1,1),param,sd.ratio)
      DR.ta = 1 - pgamma.m(qgamma.m(1-fpr.set,1,1),param,sd.ratio)}
    if(distn=="normal"){
      DR.t = 1 - pnorm(qnorm(1-FPR,0,1),sd.ratio*param,sd.ratio)
      DR.ta = 1 - pnorm(qnorm(1-fpr.set,0,1),sd.ratio*param,sd.ratio)}

    # output
    true = data.frame(FPR,DR=DR.t)
    true.all = data.frame(FPR=fpr.set,DR=DR.ta)
    ipd = data.frame(study,affected,unaffected,positive,cutoff,score)
    summdata = data.frame(N.s,n.aff,n.unaff,true.pos,false.pos,FPR=FPR.d,se.FPR,DR,se.DR,DOR,se.logDOR)
    out = list(true=true,true.all=true.all,summdata=summdata,ipd=ipd)
    return(out)
      
  }



# find gamma mean from DR50

gamma.dr =

  function(DR50,mean.unaff=1,sd.ratio=1){

    FPR50.cutoff = qgamma.m(0.5,mean.unaff,1)
    mean.aff.set = seq(mean.unaff,(mean.unaff+5),0.1)
    DR = 1-pgamma.m(FPR50.cutoff,mean.aff.set,sd.ratio)
    mean.aff = mean.aff.set[DR>=DR50]
    mean.aff = mean.aff[1]
    return(mean.aff)
  }
    


# very simple simulator, constant DOR

vsimp.diag = function(N=1000,DR50=0.5,FPR=0.05,Ns=20,p=0.1){
  
  log.DOR = logit(DR50)-logit(0.5)
  if (length(FPR)==1)
    FPR = runif(Ns,FPR/2,2*FPR)
  DR = expit(logit(FPR) + log.DOR)
  
  tsize = rep(N,Ns)
  n.aff = floor(p*tsize)
  n.unaff = tsize - n.aff
  true.pos = floor(DR*n.aff)
  false.pos = floor(FPR*n.unaff)
  
  out = data.frame(n.aff,n.unaff,true.pos,false.pos,FPR=false.pos/n.unaff,DR=true.pos/n.aff)
  return(out)
  
}
