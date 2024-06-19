############################################

# Code to apply all avaialble sequential methods to a single data set

##########################################

# Inputs:

# Est, SE, N : effect estimate, std error and trials size from each trial
# alpha, beta: Type I and II error
# delta: desired true effect (from data if NULL)

seq_ma_analysis <-

  function(Est, SE, alpha=0.05, beta=0.1, delta=NULL){

  ma.res <- metagen(Est,SE)
  N.studies <- length(Est)

  # if no desired true effect given use result from MA
  if (is.null(delta)==T)
    delta <- ma.res$TE.random
  # adjust if near zero
  delta <- ifelse(abs(delta)>0.01, delta, 0.1)

  # Apply sequential methods
  # Naive cumul MA
  naive.MA <- naive.cumul(Est,SE)

  # SMA
  ab.prior.50 <- I2.prior(SE,0.5,1)
  ab.prior.90 <- I2.prior(SE,0.9,1)
  SMA.res <- SMA(Est,SE,alpha,beta,delta,NULL,NULL)
  SMA.AB50.res <- SMA(Est,SE,alpha,beta,delta,ab.prior.50$lambda,ab.prior.50$eta)
  SMA.AB90.res <- SMA(Est,SE,alpha,beta,delta,ab.prior.90$lambda,ab.prior.90$eta)

  # TSA
  TSA.res <- TSA(Est,SE,alpha,beta,delta,NULL)

  # Combine all results
  method.names <- c("Naive","SMA","SMA I250","SMA I290","TSA")
  all.res <- bind_rows(naive.MA,SMA.res,SMA.AB50.res,SMA.AB90.res,TSA.res) %>%
    mutate(method=rep(method.names,each=N.studies))

  # Extract results at stopping or last trial
  stop.res <- all.res %>%
    filter(!status  %in% c("Inconclusive","Too few studies")) %>%
    group_by(method) %>%
    slice_head(n=1)

  out <- list(all.res=all.res, stop.res=stop.res)
  return(out)

  }

################################################################
# analyse one simulated meta-analysis
# using code above

# data: full simulated data set
# case: specific simulated meta-analysis to apply code above to
# alpha, beta, delta: desired Type I and II errors and treatment effect
# out: "stop" for results at time of stopping, "full" for complete results


one.case.analysis <-function(data, case, alpha=0.05, beta=0.1, delta=0.1, out="stop"){

  data.k <- subset(data,macase==case)

  res.k <- seq.ma.analysis(data.k$Eff.est,data.k$SE,alpha=alpha,beta=beta,delta=delta)

  full.k <- data.frame(macase=rep(case,length(res.k$all.res$method)),res.k$all.res)
  stop.k <- data.frame(macase=rep(case,length(res.k$stop.res$method)),res.k$stop.res)

  if (out=="full")
    return(full.k)
  if (out=="stop")
    return(stop.k)

  # realistic updating, max 4 updates
  if (out=="realistic"){
    ns <- max(res.k$all.res$study)
    if (ns==5)
      real.k <- subset(full.k,study %in% c(3,4,5))
    if (ns==10)
      real.k <- subset(full.k,study %in% c(5,7,9,10))
    if (ns==20)
      real.k <- subset(full.k,study %in% c(10,15,18,20))
    if (ns==50)
      real.k <- subset(full.k,study %in% c(25,40,45,50))
    if (!ns %in% c(5,10,20,50)){
      s.set <- c(floor(0.5*ns),floor(0.75*ns),floor(0.9*ns),ns)
      real.k <- subset(full.k,study %in% s.set)
    }
    real.k <- subset(real.k,!status  %in% c("Inconclusive","Too few studies"))
    real.k <- data.frame(do.call("rbind", by(real.k,real.k$method,function(x) head(x,1))))
    row.names(real.k) <- NULL

    return(real.k)
  }

}



