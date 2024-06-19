######################################################
#
# Collection of all finalised code for
# updating / sequential meta-analyses
#
#####################################################

# NOTE: The meta package (for meta-analysis) and ldbounds package (for alpha-spending boundaries) are required


###############################################
# Naive cumulative meta-analysis
###############################################

# Inputs:
# Est , SE: effect estimate, std error and sample size for each study

naive.cumul <-

  function(Est,SE){

    Nt <- length(Est)

    # meta-analysis at one update time
    one.naive <- function(i,Est,SE){
      MAi <- metagen(Est[1:i],SE[1:i],control=list(stepadj=0.5,maxiter=1000))
      out.one <- data.frame(theta=MAi$TE.random,se.theta=MAi$seTE.random,ci.low=MAi$lower.random,ci.high=MAi$upper.random,tau2=MAi$tau^2)
      return(out.one)
    }

    # all MAs
    MA.res <- map_df(1:Nt, ~one.naive(.,Est,SE))

    # conclusion of MA
    study <- 1:Nt
    status <- case_when(
      study<3 ~ "Too few studies",
      MA.res$ci.low>0 ~ "Favourable",
      MA.res$ci.high<0 ~ "Unfavourable",
      study==Nt ~ "Does not stop",
      TRUE ~ "Inconclusive"
    )

    out <- data.frame(study,MA.res,status)
    return(out)


}


###############################################
# Trial Sequential Analysis
###############################################

# Version based on statistical information only

# Inputs:
# Est, SE: effect estimate and std error from each trial
# alpha, beta: Type I and II error
# delta: desired effect size (NULL = generated from data)
# D2: estimate of D2 used to calculate sample size (NULL = generated from data)

TSA <-

  function(Est,SE,alpha=0.05,beta=0.1,delta=NULL,D2=NULL){

    Nt <- length(Est)

    # summary result
    MAall <- metagen(Est,SE,control=list(stepadj=0.5,maxiter=1000))
    delta <- ifelse(is.null(delta), MAall$TE.random, delta)

    # Meta-analysis adding one study at a time
    one.TSA <- function(i,Est,SE){

      MAi <- metagen(Est[1:i],SE[1:i],control=list(stepadj=0.5,maxiter=1000))
      cumul.V <- sum(1/(SE[1:i]^2 + MAi$tau^2))

      out.one <- data.frame(theta=MAi$TE.random,se.theta=MAi$seTE.random,tau2=MAi$tau^2,cumul.V=cumul.V)
      return(out.one)
    }

    # all MAs
    MA.res <- map_df(1:Nt, ~one.TSA(.,Est,SE))

    Zscore <- MA.res$theta / MA.res$se.theta

    # information size / fraction
    IS.stat <- inf.size(delta=delta,alpha=alpha,beta=beta,sigma=1)$IS.stat
    IF <- MA.res$cumul.V / IS.stat

    # O'Brien Fleming bounds (using lbounds package)
    obf.find <- function(i,IF,alpha){
      obf.bound <- ifelse(IF[i]<1,
                          ldBounds(IF[i],iuse=c(1,1),alpha=rep(alpha/2,2))$upper.bounds,
                          ldBounds(1,iuse=c(1,1),alpha=rep(alpha/2,2))$upper.bounds)
      return(obf.bound)
    }

    obf.bound <- map_dbl(1:Nt, ~obf.find(.,IF,alpha))

    # Determime if a stopping boundary is crossed
    study <- 1:Nt
    status <- case_when(
      study<3 ~ "Too few studies",
      Zscore > obf.bound ~ "Favourable",
      Zscore < (-obf.bound) ~ "Unfavourable",
      IF>1 ~ "No effect",
      study==Nt ~ "Does not stop",
      TRUE ~ "Inconclusive"
    )

    # Adjusted CIs
    ci.low <- MA.res$theta - obf.bound * MA.res$se.theta
    ci.high <- MA.res$theta + obf.bound * MA.res$se.theta

    out <- data.frame(study,theta=MA.res$theta,se.theta=MA.res$se.theta,ci.low,ci.high,tau2=MA.res$tau2,status)

    return(out)

}


########################################################
# Approximate Bayes Sequential Meta-Analysis
########################################################

# Inputs:
# Est, SE: effect estimate and std error from each trial
# alpha, beta: Type I and II error
# delta : desired effect size (NULL = generated from data)
# eta, lambda: parameters of IG(eta,lambda) heterogeneity prior

SMA <-

  function(Est,SE,alpha=0.05,beta=0.1,delta=NULL,lambda=NULL,eta=NULL){

    Nt <- length(Est)

    # summary result
    MAall <- metagen(Est,SE,control=list(stepadj=0.5,maxiter=1000))
    delta <- ifelse(is.null(delta)==T, MAall$TE.random, delta)

    # SMA boundaries (using support function below)
    bounds <- calc.obf.bounds(delta=delta,alpha=alpha,beta=beta)

    # for each trial in the sequential meta-analysis
    one.SMA <- function(i,Est,SE,lambda,eta){

      MAi <- metagen(Est[1:i],SE[1:i],control=list(stepadj=0.5,maxiter=1000))

      # Find approx. Bayes estimate of heterogeneity
      tau2 <- MAi$tau^2
      tau2.AB <-  ifelse(is.null(lambda)==F, max((2 * lambda + i * tau2)/(2 * eta + i - 2), 0), tau2)

      # Find cumulative Z and V
      cumul.V <- sum(1/(SE[1:i]^2 + tau2.AB))
      cumul.Z <- sum(Est[1:i]/(SE[1:i]^2 + tau2.AB))

      out.one <- data.frame(cumul.Z,cumul.V,tau2,tau2.AB)
      return(out.one)
    }

    # all MAs
    MA.res <- map_df(1:Nt, ~one.SMA(.,Est,SE,lambda,eta))

    # MA results by step
    theta <- MA.res$cumul.Z / MA.res$cumul.V
    se.theta <- sqrt(1 / MA.res$cumul.V)

    # Find repeated confidence intervals
    ci.low <- (MA.res$cumul.Z - bounds$H)/MA.res$cumul.V
    ci.high <- (MA.res$cumul.Z + bounds$H)/MA.res$cumul.V

    # Determime if a stopping boundary is crossed
    study <- 1:Nt
    status <- case_when(
      study<3 ~ "Too few studies",
      ci.low>0 ~ "Favourable",
      ci.high<0 ~ "Unfavourable",
      MA.res$cumul.V>bounds$Vmax ~ "No effect",
      study==Nt ~ "Does not stop",
      TRUE ~ "Inconclusive"
    )

    out <- data.frame(study,theta,se.theta,ci.low,ci.high,tau2=MA.res$tau2,status)
    return(out)
  }


######################################################
#Support functions
######################################################

#  Calcualte Information size for TSA ("generic" data only)

# Inputs:
# delta: desired effect size
# alpha, beta: Type I and II error
# sigma: std error (not required, default=1)

inf.size <- function(delta,alpha=0.05,beta=0.1,sigma=1,H=0){

  qfunc <- qnorm(1-alpha/2) + qnorm(1-beta)

  IS <- 4 * qfunc^2 *(sigma^2 / delta^2)
  IS.stat <- qfunc^2 / delta^2

  # adjust for heterogeneity
  AF <- 1 / (1-H)
  IS.random <- AF * IS

  out <- data.frame(IS,IS.stat,IS.random)
  return(out)
}


#####################################################
# find diversity, D2

# Inputs as above

D2 <- function(Est,SE){

  require(meta)
  Wf <- 1/SE^2
  tau2 <- metagen(Est,SE)$tau^2
  Wr <- 1 / (SE^2 + tau2)
  D2 <- 1 - (sum(Wr) / sum(Wf))
  return(D2)
}



###################################################
# calcualte O'Brien-Fleming rectangular bounds
# data taken from Higgins et al paper

# Inputs:
# delta: desired true treatment effect to detect (generic data)
# alpha, beta: Type I and II error rates

calc.obf.bounds <- function(delta,alpha=0.05,beta=0.1){

  # core data (from Higgins paper)
  H.matrix <- matrix(c(14.576,16.12,17.394,9.779,11.029,12.061,6.457,7.461,8.228),3,3)
  Vmax.matrix <- matrix(c(17.535,21.447,24.972,12.138,15.438,18.461,8.299,11.079,13.673),3,3)
  alpha.set <- c(0.001,0.01,0.05)
  beta.set <- c(0.2,0.1,0.05)

  # select correct boundaries and adjust for effect size
  H <- H.matrix[match(beta,beta.set),match(alpha,alpha.set)]
  Vmax <- Vmax.matrix[match(beta,beta.set),match(alpha,alpha.set)]

  H <- abs(H / delta)
  Vmax <- Vmax / delta^2

  out <- data.frame(H,Vmax)
  return(out)

}



######################################################
# Approx Bayes prior based on reasonable I2 from current data

# Inputs:
# SE: standard erros of effect in trials
# I2.prior: prior I2 required
# ntrials: numbers of "pseudotrials"

I2.prior <- function(SE,I2.p=0.9,ntrials=1){

  Nt <- length(SE)
  if (Nt>1){
    wt <- 1 / (SE^2)
    typvar <- (sum(wt)*(Nt-1)) / (sum(wt)^2 - sum(wt^2))
  } else{
    typvar <- SE^2
  }

  tau2.prior <- (I2.p * typvar)  / (1 - I2.p)
  eta <- ntrials/2 +1
  lambda <- ntrials*tau2.prior / 2
  out <- data.frame(eta,lambda,I2.p,typvar,tau2.prior)
  return(out)
}

