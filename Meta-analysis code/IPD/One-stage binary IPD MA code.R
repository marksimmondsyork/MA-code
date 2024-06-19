#############################################
# General code for running one-stage
# IPD meta-analyses
#############################################

# Inputs:
# outcome, must be coded 0/1 for events
# trial, any coding
# arm, treatment arm: 0=control, 1=experimental
# measure: effect estimate type, one of "RR", "OR", "MD"
# ma.type:
# "FE" all fixed effect
# "RE_treat" for random treatment effects only
# "RE_uncorr" for uncorrelated trial and treatment REs
# "RE_corr" for correlated trial/treatment REs

#############################################
# one-stage MA with no interactions

one_stage_MA <- function(outcome,trial,arm,measure="RR",ma.type="RE_treat"){

  require(lme4)

  if (measure=="OR"){
    if (ma.type=="FE")
      ma <- glm(outcome ~ factor(trial) + factor(arm),family=binomial)
    if (ma.type=="RE_treat")
      ma <- glmer(outcome ~ factor(trial) + factor(arm) + (arm-1|trial),family=binomial)
    if (ma.type=="RE_uncorr")
      ma <- glmer(outcome ~ factor(arm) + (1|trial) + (arm-1|trial),family=binomial)
    if (ma.type=="RE_corr")
      ma <- glmer(outcome ~ factor(arm) + (1+arm|trial),family=binomial)

  }

  if (measure=="RR"){
    if (ma.type=="FE")
      ma <- glm(outcome ~ factor(trial) + factor(arm),family=binomial(link="log"))
    if (ma.type=="RE_treat")
      ma <- glmer(outcome ~ factor(trial) + factor(arm) + (arm-1|trial),family=binomial(link="log"))
    if (ma.type=="RE_uncorr")
      ma <- glmer(outcome ~ factor(trial) + factor(arm) + (arm-1|trial),family=binomial(link="log"))
    if (ma.type=="RE_both")
      ma <- glmer(outcome ~ factor(arm) + (1+arm|trial),family=binomial(link="log"))
  }

  if (measure=="MD"){
    if (ma.type=="FE")
      ma <- lm(outcome ~ factor(trial) + factor(arm))
    if (ma.type=="RE")
      ma <- lmer(outcome ~ factor(trial) + factor(arm) + (arm-1|trial))
    if (ma.type=="RE_both")
      ma <- lmer(outcome ~ factor(arm) + (1+arm|trial))
  }

  return(summary(ma))

}



###############################################################
# one-stage meta-analysis with treatment-covariate interaction

# Inputs as above except:
# covar: covariate value, NOTE binary covariates must be coded 0/1
# centre: should continuous covariates be centred around the average?
#         generally set as TRUE for continuous and FALSE for binary covariates
# ma.type has new options:
#         "RE_inter", for adding random interaction effects
#         "EcoBias", for dividing within/between trials data for covariate

one_stage_MA_inter <- function(outcome,trial,arm,covar,measure="RR",ma.type="RE_treat",centre=TRUE){

  require(lme4)

  if(centre==T)
    covar <- covar - mean(covar,na.rm=T)
  
  if (ma.type=="EcoBias"){
    trial.size <- tapply(trial,trial,length)
    covar.between <- tapply(covar,trial,mean,na.rm=TRUE)
    covar.between <- rep(covar.between,trial.size)
    covar.within <- covar - covar.between
  }

  if (measure=="OR"){
    if (ma.type=="FE")
      ma <- glm(outcome ~ factor(trial) + factor(arm)*covar,family=binomial)
    if (ma.type=="RE_treat")
      ma <- glmer(outcome ~ factor(trial) + factor(arm)*covar + (arm-1|trial),family=binomial)
    if (ma.type=="RE_uncorr")
      ma <- glmer(outcome ~ factor(arm)*covar + (1|trial) + (arm-1|trial),family=binomial)
    if (ma.type=="RE_corr")
      ma <- glmer(outcome ~ factor(arm)*covar + (1+arm|trial),family=binomial)
    if (ma.type=="RE_inter")
      ma <- glmer(outcome ~ factor(arm)*covar + (1|trial) + (arm-1|trial) + (arm:covar-arm-1|trial),family=binomial)
    if (ma.type=="EcoBias")
      ma <- glmer(outcome ~ factor(arm)*(covar.within+covar.between) + (1+arm|trial),family=binomial)
  }

  if (measure=="RR"){
    if (ma.type=="FE")
      ma <- glm(outcome ~ factor(trial) + factor(arm)*covar,family=binomial(link="log"))
    if (ma.type=="RE_treat")
      ma <- glmer(outcome ~ factor(trial) + factor(arm)*covar + (arm-1|trial),family=binomial(link="log"))
    if (ma.type=="RE_uncorr")
      ma <- glmer(outcome ~ factor(arm)*covar + (1|trial) + (arm-1|trial),family=binomial(link="log"))
    if (ma.type=="RE_corr")
      ma <- glmer(outcome ~ factor(arm)*covar + (1+arm|trial),family=binomial(link="log"))
    if (ma.type=="RE_inter")
      ma <- glmer(outcome ~ factor(arm)*covar + (1|trial) + (arm-1|trial) + (arm:covar-arm-1|trial),family=binomial(link="log"))
    if (ma.type=="EcoBias")
      ma <- glmer(outcome ~ factor(arm)*(covar.within+covar.between) + (1+arm|trial),family=binomial(link="log"))
    
  }

  if (measure=="MD"){
    if (ma.type=="FE")
      ma <- lm(outcome ~ factor(trial) + factor(arm)*covar)
    if (ma.type=="RE_treat")
      ma <- lmer(outcome ~ factor(trial) + factor(arm)*covar + (arm-1|trial))
    if (ma.type=="RE_uncorr")
      ma <- lmer(outcome ~ factor(arm)*covar + (1|trial) + (arm-1|trial))
    if (ma.type=="RE_corr")
      ma <- lmer(outcome ~ factor(arm)*covar + (1+arm|trial))
    if (ma.type=="RE_inter")
      ma <- lmer(outcome ~ factor(arm)*covar + (1|trial) + (arm-1|trial) + (arm:covar-arm-1|trial))
    if (ma.type=="EcoBias")
      ma <- lmer(outcome ~ factor(arm)*(covar.within+covar.between) + (1+arm|trial))
  }

  return(summary(ma))

}


##################################################################
# convert log odds/risks into standard form with CI

log_convert <- function(result){

  coeffs <- result$coefficients
  out <- data.frame(Estimate=exp(coeffs[,1]),CIlow=exp(coeffs[,1]-1.96*coeffs[,2]),CIhigh=exp(coeffs[,1]+1.96*coeffs[,2]))
  out <- round(out,3)

  return(out)
}
