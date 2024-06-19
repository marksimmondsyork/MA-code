survivalweibull <-

  function (trialsize=rep(200,10),p=2,theta=0,phi=0,mu=0,gamma=0,w=10,het=0.00001) {

# generates a set of trials (from trialsize) with weibull survival times with treatment effect theta and a covariate age
# p is shape parameter, theta treatment effect, phi trial effect, mu age effect, gamma age-treatment interaction
# w width of each age distrbution, het heterogeneity of treatment effect

    n <- sum(trialsize)
    numtrials <- length(trialsize)
    trial <- rep(1:numtrials,trialsize)
    meanage <- rep(0,numtrials)
    id <- 1:n

#generate treatment and age variables

    treat <- rep(c(0,1),length.out=n)
    lower <- meanage-w
    upper <- meanage+w
    meanage1 <- rep(meanage,trialsize)
    age <- runif(n,rep(lower,trialsize),rep(upper,trialsize))
    agec <- age-mean(age)
    meanagec <- meanage1-mean(age)

    thetar <- rep(rnorm(numtrials,theta,sqrt(het)),trialsize)
    phir <- rep(rnorm(numtrials,phi,1),trialsize)

# generate failure times

    scale<-exp((1/p)*(phir+thetar*treat+mu*agec+gamma*agec*treat))
    s <- rweibull(n,p,1/scale)

#cmax is maximum time (end of trial), vis is visibility (1 observed, 0 censored), censtime is observed time

    cmax <- rep(2*mean(s),n)
    event <- eventtime <- rep(0,n)

#generate observed times and random censoring

    rand <- runif(n,0,1)
    for (i in 1:n){
      if (s[i]<cmax[i]){
        event[i] <- 1
        eventtime[i] <- s[i]
      }
      else{
        event[i] <- 0
        eventtime[i] <- cmax[i]
      }
      if (rand[i]<0.05){
        eventtime[i] <- runif(1,0,eventtime[i])
        event[i]<-0
      }
    }
    
    out <- data.frame(id,eventtime,event,treat,trial,age=agec,meanage=meanagec)
    return(out)

  }
