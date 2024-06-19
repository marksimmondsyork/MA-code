##################################################
# basic code to calculate sens/spec etc.
# corrected using meta
####################################################

diag_data_2x2 <- function(Study,TP,FP,TN,FN){

  require(meta)
  require(purrr)
  
  expit <- function(x) {exp(x)/(1+exp(x))}

  # calculate sensitivity
  sens.ma  <- metaprop(TP, (TP+FN), sm="PLOGIT")
  sens.data <- data.frame(sens=expit(sens.ma$TE),sens.low=sens.ma$lower,sens.high=sens.ma$upper)
  
  # specificity
  spec.ma <- metaprop(TN, (TN+FP), sm="PLOGIT")
  spec.data <- data.frame(spec=expit(spec.ma$TE),spec.low=spec.ma$lower,spec.high=spec.ma$upper) 
  
  # positive rate
  #PR.ma  <- metaprop((TP+FN), (TP+FN+TN+FP), sm="PLOGIT")
  #PR.data <- data.frame(PR=expit(PR.ma$TE),PR.low=PR.ma$lower,PR.high=PR.ma$upper)

  sensspec <- 100 * cbind(sens.data,spec.data)

  # PPV/NPV
  PPV.ma  <- metaprop(TP, (TP+FP), sm="PLOGIT")
  PPV.data <- data.frame(PPV=expit(PPV.ma$TE),PPV.low=PPV.ma$lower,PPV.high=PPV.ma$upper)
  
  NPV.ma  <- metaprop(TN, (TN+FN), sm="PLOGIT")
  NPV.data <- data.frame(NPV=expit(NPV.ma$TE),NPV.low=NPV.ma$lower,NPV.high=NPV.ma$upper)

  PPVNPV <- 100 * cbind(PPV.data,NPV.data)

  # Diagnostic odds ratio
  
  one_DOR <- function(i,data){
    reg.dor <- summary(glm(cbind(positive,n-positive) ~ factor(affected),data=subset(data,test==i),family=binomial(link="logit")))
    out <- data.frame(lDOR=reg.dor$coefficients[2,1], selDOR=reg.dor$coefficients[2,2])
    return(out)
  }
  
  ntests <- length(FP)
  data2 <- data.frame(n=as.numeric(t(cbind(TP+FN,TN+FP))), 
                      positive=as.numeric(t(cbind(TP,FP))), 
                      test=rep(1:ntests,each=2), 
                      affected=rep(c(1,0),ntests))
  
  all.DORS <- map_dfr(1:ntests, ~one_DOR(.,data2))

  dor.data <- data.frame(lDOR=all.DORS$lDOR,
                         selDOR=all.DORS$selDOR,
                         DOR=exp(all.DORS$lDOR),
                         DOR.low=exp(all.DORS$lDOR - 1.96*all.DORS$selDOR),
                         DOR.high=exp(all.DORS$lDOR + 1.96*all.DORS$selDOR))

  # all data
  # sensspec: sensitivity/specificity
  #PPVNPV: PPV and NPV
  #dor.data: diagnostic odds ratios
  out <- data.frame(Study=as.character(Study),N=TP+TN+FP+FN,sensspec,PPVNPV,dor.data)
  return(out)

}


########################################################
# calculate 2x2 data from summary results

summary_to_2x2 <- function(N,Sens,Sens.low,Sens.high,Spec,tol=2){
  
  Sens.se <- (Sens.high/100-Sens.low/100) / (2*1.96)
  
  # calculate most plausible Npos and Nneg values
  Npos <- round((Sens/100*(1-Sens/100) / Sens.se^2), 0)
  Npos <- (Npos-tol):(Npos+tol)
  Nneg <- N - Npos
  
  # calculate 2x2 given this
  TP <- round(Sens/100 * Npos, 0)
  TN <- round(Spec/100 * Nneg, 0)
  FP <- Nneg - TN
  FN <- Npos - TP
  
  # recalculate Sens/Spec
  new.res <- diag_data_2x2(as.character(Npos),TP,TN,FP,FN)[,3:10]
  out <- data.frame(TP,TN,FP,FN,new.res)
  return(out)
  
}

