#############################################
# General code for running two-stage
# IPD meta-analyses
#############################################


# Inputs:
# outcome, must be coded 0/1 for events
# trial, any coding
# arm, treatment arm: 0=control, 1=experimental
# measure, one of RR, OR, MD, SMD
# method, method.tau: MA options, see metabin code in meta library
# filename, filename for forest plot (NULL if no forest plot required)


#####################################################
# Two-stage using data tables

two_stage_MA <- function(outcome,trial,arm,measure="RR",method="Inverse",method.tau="DL",filename=NULL){

  require(meta)
  trialnames <- as.character(unique(trial))
  trialnames <- trialnames[order(trialnames)]
  ntrials <- length(trialnames)

  num.patients <- tapply(outcome, factor(arm):factor(trial), function(x) length(x[!is.na(x)]))
  num.patients <- ifelse(is.na(num.patients),0,num.patients)
  num.patients <- matrix(num.patients,ntrials,2)

  # binary data
  if (measure %in% c("RR","OR")){

    num.events <- tapply(outcome, factor(arm):factor(trial), function(x) sum(x[!is.na(x)]))
    num.events <- ifelse(is.na(num.events),0,num.events)
    num.events <- matrix(num.events,ntrials,2)

    ma <- metabin(num.events[,2],num.patients[,2],num.events[,1],num.patients[,1],sm=measure,method=method,method.tau=method.tau,studlab=trialnames)
  }

  # continuous data
  if (measure %in% c("MD","SMD")){

    mean.outcome <- tapply(outcome, factor(arm):factor(trial), function(x) mean(x[!is.na(x)]))
    mean.outcome <- matrix(mean.outcome,ntrials,2)
    sd.outcome <- tapply(outcome, factor(arm):factor(trial), function(x) sd(x[!is.na(x)]))
    sd.outcome <- matrix(sd.outcome,ntrials,2)

    ma <- metacont(num.patients[,2],mean.outcome[,2],sd.outcome[,2],num.patients[,1],mean.outcome[,1],sd.outcome[,1],sm=measure,method.tau=method.tau,studlab=trialnames)
  }

  # forest plot
  if(!is.null(filename)){
    filename <- paste(filename,"forest plot.png",sep=" ")
    png(filename,width=800,height=500)
    forest(ma,studlab=trialnames,xlab=measure,rightcols=c("effect", "ci"))
    dev.off()
  }

  return(ma)

}

########################################
# two-stage subgroup analysis
# patient or trial-level analyses possible

# inputs as above except:
# covar: subgroup identifier, can be any text/number

two_stage_subgroups <- function(outcome,trial,arm,covar,measure="RR"){

  data1 <- data.frame(outcome,trial,arm,covar)
  data1s <- split(data1,factor(data1$covar))
  ma.sbg <- lapply(data1s,function(x) two_stage_MA(x$outcome,x$trial,x$arm,measure))

  # extract simplified results
  one_analy <- function(mai){
    data.frame(NTrials=mai$k,Npats=sum(mai$n.e+mai$n.c),Nevents=sum(mai$event.e + mai$event.c),
               Est=exp(mai$TE.random),CIlow=exp(mai$lower.random),CIhigh=exp(mai$upper.random),I2=mai$I2)
  }

  results <- lapply(ma.sbg, one_analy)
  results <- do.call(rbind,results)

  out <- list(MA=ma.sbg,results=results)
  return(out)
}


################################################################
# two-stage meta-analysis of interactions

# Inputs as above

two_stage_MAoI <- function(outcome,trial,arm,covar,measure="RR"){

  require(meta)

  # regression analysis in one trial
  one_trial_int_model <- function(trial.num,data,measure){
    datai <- subset(data,trial==trial.num)

    if (measure=="OR")
      res <- glm(outcome ~ factor(arm)*covar,family=binomial,data=datai)
    if (measure=="RR")
      res <- glm(outcome ~ factor(arm)*covar,family=binomial(link="log"),data=datai)
    if (measure=="MD")
      res <- lm(outcome ~ factor(arm)*covar,data=datai)

    res.coeffs <- data.frame(Param=c("Treatment","Interaction"),
                             summary(res)$coefficients[c(2,4),])
    return(res.coeffs)
  }

  # run model in all trials
  all.data <- data.frame(outcome,trial,arm,covar)
  trial.names <- unique(all.data$trial)
  all.res <- lapply(trial.names, function(x) one_trial_int_model(x,all.data,measure))
  all.res <- do.call(rbind,all.res)
  all.res <- data.frame(Trial=rep(trial.names,each=2),all.res)

  # MA of treatment effect
  ma.data <- subset(all.res,Param=="Treatment")
  ma.res <- metagen(Estimate,Std..Error,data=ma.data,studlab=Trial)

  # MA of interactions
  ma.int.data <- subset(all.res,Param=="Interaction")
  ma.int.res <- metagen(Estimate,Std..Error,data=ma.int.data,studlab=Trial)

  out <- list(MA.res=ma.res, MAoI.res=ma.int.res)
  return(out)

}


##################################################################
# convert log odds/risks into standard form with CI

log_convert <- function(result){

  coeffs <- result$coefficients
  out <- data.frame(Estimate=exp(coeffs[,1]),CIlow=exp(coeffs[,1]-1.96*coeffs[,2]),CIhigh=exp(coeffs[,1]+1.96*coeffs[,2]))
  out <- round(out,3)

  return(out)
}
