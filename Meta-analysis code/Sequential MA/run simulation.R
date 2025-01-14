##############################################
# simulate meta-analyses using conventional cumulative MA
############################################

library(tidyverse)
library(meta)
library(ldbounds)
library(furrr)

source("ma simulator.R")
source("all sequential code.R")
source("universal analysis code.R")

# run one simulation
simul_many_cases_2 <- function(x){
  datai <- ma_simul_gen(N[x],N.trials[x],true.eff[x],I2[x],SD[x],tsize="different")
  out <- seq_ma_analysis(datai$Eff.est, datai$SE, alpha=0.05, beta=0.1, delta=NULL)$stop.res
  return(out)
}

inf.size(alpha=0.05,beta=0.1,delta=0.1,sigma=1,H=0.5)

####################################################
# simulate data, various cases

# parameter set-up

nsims <- 5
N.1 <- 9000
N.trials.1 <- c(5,10,20,50)
true.eff.1 <- c(0,0.1)
I2.1 <- c(0,0.25,0.5,0.75,0.9)
SD.1 <- 1
tsize.1 <- "same"

all.params <- expand.grid(N=N.1,N.trials=N.trials.1,true.eff=true.eff.1,
                          I2=I2.1,SD=SD.1,tsize=tsize.1) %>%
  mutate(Sim.case=1:n(),
         nsims=nsims) %>%
  uncount(nsims) %>%
  mutate(MA.case=1:n())

n.mas <- length(all.params$N.trials)
N <- all.params$N
N.trials <- all.params$N.trials
true.eff <- all.params$true.eff
I2 <- all.params$I2
SD <- all.params$SD
tsize <- all.params$tsize

# run simulation

plan(multisession,workers=5)

MA.sim.res.1 <- future_map_dfr(1:n.mas, simul_many_cases_2)

plan(sequential)


write_csv(MA.sim.res.1,"Sequential simulation results April 22 B.csv")
write_csv(all.params,"Parameter set April 22 B.csv")




