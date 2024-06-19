# "one stage" diagnostic MA models

# Inputs:
# TP, FP, TN, FN: true positive value etc. one row per study

diagnostic_MA <-

  function(TP,FP,TN,FN){

    require(lme4)
    require(msm)

    # convert data
    data2 <- data.frame(n=as.numeric(t(cbind(TP+FN,FP+TN))), 
                        true=as.numeric(t(cbind(TP,TN))), 
                        study=rep(1:length(TP),each=2), 
                        sens=rep(c(1,0),length(TP)),
                        spec=rep(c(0,1),length(TP)))

    # run logistic regression model (this does the MA)
    res <- glmer(cbind(true,n-true) ~ 0 + sens + spec + (0+sens + spec|study),
                 data=data2,family=binomial(link="logit"))

    # extract diagnostic accuracy parameters
    # logit-transformed sens/spec and their SEs
    logit.sens <- summary(res)$coefficients[1,1]
    logit.spec <- summary(res)$coefficients[2,1]
    logit.sens.se <- summary(res)$coefficients[1,2]
    logit.spec.se <- summary(res)$coefficients[2,2]
    
    # covariance and heterogeneity parameters
    cov.mat <- summary(res)$vcov
    cov.sens.spec <- cov.mat[1,2]
    tau.sens <- as.numeric(attr(summary(res)$varcor$study,"stddev")[1])
    tau.spec <-  as.numeric(attr(summary(res)$varcor$study,"stddev")[2])
    rho <-  as.numeric(attr(summary(res)$varcor$study,"correlation")[2,1])
    
    # output data frame
    params <- data.frame(logit.sens,logit.sens.se,
                         logit.spec,logit.spec.se,
                         cov.sens.spec,tau.sens,tau.spec,rho)
    
    # Summary sensitivity and specificity (what you usually need)
    cithresh <- qnorm(0.975)
    sens.spec <- data.frame(Sens=plogis(logit.sens),
                          Sens.low=plogis(logit.sens-cithresh*logit.sens.se),
                          Sens.high=plogis(logit.sens+cithresh*logit.sens.se),
                          Spec=plogis(logit.spec),
                          Spec.low=plogis(logit.spec-cithresh*logit.spec.se),
                          Spec.high=plogis(logit.spec+cithresh*logit.spec.se))
    
    # DOR and likelihood ratios
    logDOR <- logit.sens+logit.spec
    logLRpos <- plogis(logit.sens) / (1-plogis(logit.spec))
    logLRneg <- ((1-plogis(logit.sens)) / plogis(logit.spec))
    
    se.logDOR <- deltamethod (~ (x1+x2), c(logit.sens,logit.spec),cov.mat)
    se.logLRpos = deltamethod (~ log( (exp(x1)/(1+exp(x1)))/(1-(exp(x2)/ (1+exp(x2))))), c(logit.sens,logit.spec),cov.mat)
    se.logLRneg = deltamethod (~ log( (1-(exp(x1)/(1+exp(x1))))/ (exp(x2)/(1+exp(x2)))), c(logit.sens,logit.spec),cov.mat)
    
    # output data frame
    other.results <- data.frame(DOR=exp(logDOR),
                                DOR.low=exp(logDOR-cithresh*se.logDOR),
                                DOR.high=exp(logDOR+cithresh*se.logDOR),
                                LRpos=exp(logLRpos),
                                LRpos.low=exp(logLRpos-cithresh*se.logLRpos),
                                LRpos.lhigh=exp(logLRpos+cithresh*se.logLRpos),
                                LRneg=exp(logLRneg),
                                LRneg.low=exp(logLRneg-cithresh*se.logLRneg),
                                LRneg.lhigh=exp(logLRneg+cithresh*se.logLRneg))
    
    # Rutter-Gatsonis model parameters (using Harbord 2007)
    beta <- log(tau.spec/tau.sens)
    alpha <- (tau.spec/tau.sens)^0.5 * logit.sens + (tau.sens/tau.spec)^0.5 * logit.spec
    
    # HSROC curve data from R+G model
    fpr.set <- c(seq(0.001,0.049,0.001),seq(0.05,0.99,0.01))
    sens.HSROC <- plogis( qlogis(fpr.set)*exp(-beta) + alpha*exp(-beta/2) )
    hsroc <- data.frame(Spec=(1-fpr.set),Sens=sens.HSROC)
    
    out <- list(res=summary(res),
                params=params,
                sens.spec=sens.spec,
                other.results=other.results,
                hsroc=hsroc)
    
    return(out)
  }



