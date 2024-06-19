library(tidyverse)
library(magrittr)
library(meta)

##############################################
# MA of single estimate (cont.)

simple_MA_data <- function(Nt=10,Nmin=100,Nmax=1000,mu=0,tau2=1,epsilon=1,thetain=NULL){
  
  # parameters
  if (is.null(thetain))
    theta.i <- rnorm(Nt, mu, sqrt(tau2))
  else
    theta.i <- theta.in

  # IPD
  trialsize <- floor(runif(Nt, Nmin, Nmax))
  epsilon.i <- rnorm(Nt, epsilon, 0.1*epsilon)
  Y <- pmap(list(trialsize,theta.i,epsilon.i), rnorm)
  
  # summary data
  Ybar <- map_dbl(Y,mean)
  Ysd <- map_dbl(Y,sd)
  Ysem <- Ysd / sqrt(trialsize)

  out <- data.frame(N=trialsize,Ybar,Ysd,Ysem,mu.true=mu,tau2.true=tau2,epsilon)
  return(out)
}



#################################################
# MA of proportions (binary)

simple_binary_MA_data <- function(Nt=10,Nmin=100,Nmax=1000,mu=0.1,tau2=0.05,LOin=NULL){

  if (is.null(LOin))
    theta.i <- rnorm(Nt, mu, sqrt(tau2))
  else
    theta.i <- LOin

  p.i <- exp(theta.i) / (1+exp(theta.i))
  trialsize <- floor(runif(Nt, Nmin, Nmax))
  Y <- map2_dbl(trialsize, p.i, ~rbinom(1,.x,.y))
  
  # summaries (Log odds)
  p.est <- Y / trialsize
  Ybar <- log(p.est / (1-p.est))
  Ysem <- sqrt( 1 / (trialsize * p.est * (1-p.est)) )

  out <- data.frame(N=trialsize,Ybar,Ysem,mu.true=mu,tau2.true=tau2)
  return(out)
}

#########################################################
# version with "small sample" effects

bias_MA_data <- function(Nt=10,Nmin=100,Nmax=1000,mu=0,tau2=1,epsilon=1,trend=0.001){
  
  require(purrr)

  theta.i <- rnorm(Nt, mu, sqrt(tau2))

  # IPD
  trialsize <- floor(runif(Nt, Nmin, Nmax))
  epsilon.i <- rnorm(Nt, epsilon, 0.1*epsilon)
  trialsize.c <- trialsize - mean(trialsize)
  theta.i.bias <- theta.i - trend * trialsize.c
 
  Y <- pmap(list(trialsize,theta.i.bias,epsilon.i), rnorm)
  
  # summary data
  Ybar <- map_dbl(Y,mean)
  Ysd <- map_dbl(Y,sd)
  Ysem <- Ysd / sqrt(trialsize)
  
  out <- data.frame(N=trialsize,Ybar,Ysd,Ysem,mu.true=mu,tau2.true=tau2,epsilon)
  return(out)
}
