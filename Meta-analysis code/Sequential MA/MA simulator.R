#####################################################
# simulating meta-analyses for large-scale study
#####################################################

# generic one MA simulator with specified total sample size

ma_simul_gen <- function(N,N.trials,true.eff,I2=0,SD=1,tsize="different"){

  Trial.size <- rep(floor(N/N.trials),N.trials)
  if(tsize=="different")
    Trial.size <- stick_tsize(N,N.trials)

  SE <- 2 * SD / sqrt(Trial.size)
  tau2 <- I2tohet(I2,SE)
  true.y <- rnorm(N.trials, true.eff, sqrt(tau2))
  Eff.est <- rnorm(N.trials, true.y, SE)
  out <- data.frame(N,N.trials,true.eff,I2,SD,tsize,Trial.size, Eff.est, SE,tau2=rep(tau2,N.trials))
  return(out)
}


# stick-breaking trial-size calculator

stick_tsize <- function(N,N.trials,alpha=NULL,beta=NULL){

  alpha <- ifelse(is.null(alpha)==T, N.trials, alpha)
  beta <- ifelse(is.null(beta)==T, 1, beta)
  prop.trial <- rbeta(N.trials, beta, alpha)
  prop.trial <- prop.trial / sum(prop.trial)
  Trial.size <- floor(N * prop.trial)
  leftover <- N-sum(Trial.size)
  Trial.size <- Trial.size + ceiling(leftover / N.trials)
  Trial.size <- ifelse(Trial.size<10,10,Trial.size)
  return(Trial.size)
}


# convert I2 to heterogeneity

I2tohet <- function(I2,SE){
  Nt <- length(SE)
  w <- 1/SE^2
  typvar <- ((Nt-1) * sum(w)) / ((sum(w))^2 - sum(w^2))
  tau2 <- (typvar * I2) / (1 - I2)
  return(tau2)
}


