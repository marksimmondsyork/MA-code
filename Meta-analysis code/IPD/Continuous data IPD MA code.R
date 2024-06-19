#########################################
# One and two-stage code for SMD/MD MA
########################################

# two-stage

two_stage_cont_MA <- function(baseline,final,study,arm,method="ancova",scale="MD"){

  require(meta)
  require(magrittr)
  require(tidyverse)

  data <- data.frame(baseline,final,study,arm)
  data$change <- data$final - data$baseline
  data <- subset(data,!is.na(baseline) & !is.na(final))
  
  tnames <- unique(data$study)
  Nt <- length(tnames)

  # data for model-based analyses 
  if (!method %in% c("final","change")){
    
    # scaling by sd (testing)
    studysize <- data %$% tapply(study,study,length)
    Nt <- length(studysize)
    sd.base <- data %$% tapply(baseline,factor(study),sd,na.rm=T)
    scale.factor <- rep(sd.base,studysize)
    data.smd <- data %>%
      mutate(baseline=baseline/scale.factor,final=final/scale.factor,change=change/scale.factor)
  }

  # final score MA (standard using metacont)
  if (method=="final"){
    n.gp <- data %$% t(matrix( tapply(baseline, factor(study):factor(arm), length),2,Nt))
    mean.final.gp <- data %$% t(matrix(tapply(final, factor(study):factor(arm), mean,na.rm=T),2,Nt))
    sd.final.gp <- data %$% t(matrix(tapply(final, factor(study):factor(arm), sd,na.rm=T),2,Nt))
    
      ma <- metacont(n.gp[,2],mean.final.gp[,2],sd.final.gp[,2],n.gp[,1],mean.final.gp[,1],sd.final.gp[,1],sm=scale,studlab=tnames)
  }
  
  # change from baseline MA (standard using metacont)
  if (method=="change"){
    n.gp <- data %$% t( matrix( tapply(baseline, factor(study):factor(arm), length),2,Nt))
    mean.change.gp <- data %$% t(matrix(tapply(change,factor(study):factor(arm),mean,na.rm=T),2,Nt))
    sd.change.gp <- data %$% t(matrix(tapply(change,factor(study):factor(arm),sd,na.rm=T),2,Nt))
  
    ma <- metacont(n.gp[,2],mean.change.gp[,2],sd.change.gp[,2],n.gp[,1],mean.change.gp[,1],sd.change.gp[,1],sm=scale,studlab=tnames)
      
  }
  
  # within-study ANCOVA
  if (method=="ancova"){
    ancova_fun <- function(data){
      summary(lm(final~baseline+factor(arm),data=data))$coefficients
    }
    n.gp <- data %$% t( matrix( tapply(baseline, factor(study):factor(arm), length),2,Nt))
    
    if (scale=="MD"){
      ancova.md <- do.call(rbind, by(data,data$study,ancova_fun))
      ancova.md <- ancova.md[grepl("arm",rownames(ancova.md)),]
      ma <- metagen(ancova.md[,1],ancova.md[,2],studlab=tnames,n.e=n.gp[,2],n.c=n.gp[,1])
    }
    if (scale=="SMD"){
      ancova.smd <- do.call(rbind, by(data.smd,data.smd$study,ancova_fun))
      ancova.smd <- ancova.smd[grepl("arm",rownames(ancova.smd)),]
      ma <- metagen(ancova.smd[,1],ancova.smd[,2],studlab=tnames,n.e=n.gp[,2],n.c=n.gp[,1])
    }
  }
  
  # change from baseline, modelled
  if (method=="change_model"){
    change_fun <- function(data){
      summary(lm(change~factor(arm),data=data))$coefficients
    }
    if (scale=="MD"){
      change.md <- do.call(rbind, by(data,data$study,change_fun))
      change.md <- change.md[grepl("arm",rownames(change.md)),]
      ma <- metagen(change.md[,1],change.md[,2],studlab=tnames)
    }
    if (scale=="SMD"){
      change.smd <- do.call(rbind, by(data.smd,data.smd$study,change_fun))
      change.smd <- change.smd[grepl("arm",rownames(change.smd)),]
      ma <- metagen(change.smd[,1],change.smd[,2],studlab=tnames)
    }
  }

  # final score, modelled
  if (method=="final_model"){
    final_fun <- function(data){
      summary(lm(final~factor(arm),data=data))$coefficients
    }
    if (scale=="MD"){
      final.md <- do.call(rbind, by(data,data$study,final_fun))
      final.md <- final.md[grepl("arm",rownames(final.md)),]
      ma <- metagen(final.md[,1],final.md[,2],studlab=tnames)
    }
    if (scale=="SMD"){
      final.smd <- do.call(rbind, by(data.smd,data.smd$study,final_fun))
      final.smd <- final.smd[grepl("arm",rownames(final.smd)),]
      ma <- metagen(final.smd[,1],final.smd[,2],studlab=tnames)
    }
  }

  return(ma)

}


#########################################################
# all time points, mixed model

two_stage_repmeasures_MA <- function(outcome,baseline,time,study,arm){
  
  require(lme4)
  require(meta)
  
  # functions to fit mixed effect models
  
  mix_fun <- function(data.one){
    res <-lmer(outcome ~ baseline + arm*factor(time) -arm + (1+arm|time), data=data.one)
    res2 <- summary(res)$coefficients[,1:2] 
    var1 <- rownames(res2)
    res.out <- data.frame(var1,res2)
    return(res.out)
  }
  
  single_fun <- function(data.one){
    res <- lm(outcome ~ baseline + arm, data=data.one)
    res2 <- summary(res)$coefficients[,1:2]
    res.out <- data.frame(var1=as.character(unique(data.one$time)),res2)
    return(res.out)
  }
  
  rep_MA_one <- function(data.one){
    ld <- length(unique(data.one$time))
    rep.meas <- if(ld>1)
      mix_fun(data.one) else
        single_fun(data.one)
    return(rep.meas)
  }
  
  # apply to each study
  data1 <- data.frame(outcome,baseline,time,study,arm)
  data1 <- subset(data1,!is.na(baseline) & !is.na(outcome))
  
  # construct results
  rep.meas.all <- by(data1,data1$study,rep_MA_one)
  lres <- lapply(rep.meas.all,function(x) length(x[,1]))
  names2 <- do.call(c, map(rep.meas.all,~rownames(.)))
  
  rep.meas2 <- do.call(rbind, rep.meas.all) %>%
    as.data.frame() %>%
    mutate(Study=rep(names(rep.meas.all),lres)) %>%
    mutate(variable=names2) %>%
    filter(grepl("arm",variable)) %>%
    mutate(Year=str_sub(var1,-1))
  
  # meta-analysis by year
  ma_fun <- function(res.data){
    ma <- metagen(Estimate,Std..Error,studlab=Study,data=res.data)
  }
  
  ma.res <- by(rep.meas2,rep.meas2$Year,ma_fun)
  
  return(ma.res)
  
}


##################################################
##################################################
# one stage analysis

# single time point

one_stage_cont_MA <- function(baseline,final,study,arm,method="ancova",scale="MD",ma.type="RE"){
  
  require(lme4)
  change <- final - baseline
  
  # convert data
  if (scale=="MD"){
    data.ma <- data.frame(study,arm,baseline,final,change)
  }
  
  if (scale=="SMD"){
    # scaling by sd (testing)
    studysize <- tapply(study,study,length)
    Nt <- length(studysize)
    sd.base <- tapply(baseline,factor(study),sd,na.rm=T)
    scale.factor <- rep(sd.base,studysize)
    data.ma <- data.frame(study,arm,baseline=baseline/scale.factor,final=final/scale.factor,change=change/scale.factor)
  }
  
  data.ma <- subset(data.ma,!is.na(baseline) & !is.na(final))
  
# analysis models
  if (method=="change"){
    if (ma.type=="RE")
      ma <- lmer(change ~ factor(study) + factor(arm) + (arm-1|study), data=data.ma)
    if (ma.type=="RE_both")
      ma <- lmer(change ~ factor(arm) + (1+arm|study), data=data.ma)
  }
  
  if (method=="final"){
    if (ma.type=="RE")
      ma <- lmer(final ~ factor(study)+ factor(arm) + (arm-1|study), data=data.ma)
    if (ma.type=="RE_both")
      ma <- lmer(final ~ factor(arm) + (1+arm|study), data=data.ma)
  }
  
  if (method=="ancova"){
    if (ma.type=="RE")
      ma <- lmer(final ~ factor(study) + baseline + factor(arm) + (arm+baseline-1|study), data=data.ma)
    if (ma.type=="RE_both")
      ma <- lmer(final ~ baseline + factor(arm) + (1+baseline+arm|study), data=data.ma)
  }
  
  return(ma)
 
}
  

###############################################################
# repeated measures mixed model

one_stage_repmeasures_MA <- function(childID,outcome,baseline,time,study,arm){
  
  require(lme4)
  
  data.ma <- data.frame(childID,outcome,baseline,time,study,arm)
  data.ma <- subset(data.ma,!is.na(baseline) & !is.na(outcome))
  
  ma <-lmer(outcome ~ baseline + arm*factor(time) - arm + (1|childID/study) + (baseline+arm-1|study), data=data.ma)
  return(ma)

  }
