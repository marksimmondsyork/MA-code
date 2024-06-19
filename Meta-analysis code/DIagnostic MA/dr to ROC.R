# generate symmetric ROC from DR and FPR

source("C:/mark/r/logit.r")

drtoROC = function(dr,fpr){
  
  dor = (dr*(1-fpr)) / ((1-dr)*fpr)
  
  fpr.set = seq(0,1,0.01)
  dr = expit(logit(fpr.set)+log(dor))
  out = data.frame(fpr=fpr.set,dr)
  return(out)
  
}

roc1 = drtoROC(0.32,0.25)
roc2 = drtoROC(0.96,0.13)

library(ggplot2)
source("C:/mark/r/ggplot mytheme.r")

pdata = data.frame(FPR=100*c(roc1$fpr,roc2$fpr),DR=100*c(roc1$dr,roc2$dr),Method=rep(c("Plasma amyloid beta","Age"),each=length(roc1$dr)))


# useless test line
useless = geom_line(data=data.frame(x=0:100,y=0:100),aes(x=x,y=y),linetype="dashed",colour="grey",size=0.9)
useless.text <- geom_text(aes(label="Useless test",x=13,y=11,hjust=0,yjust=0),data=data.frame(),colour="black",size=6)

xaxis <- scale_x_continuous("False-positive rate (%)",limits=c(0,100),breaks=seq(0,100,20))
yaxis <- scale_y_continuous("Detection rate (%)",limits=c(0,100),breaks=seq(0,100,20))

g1 = ggplot(pdata,aes(x=FPR,y=DR,colour=Method))
g1 + geom_line(size=2) + xaxis + yaxis + useless + useless.text + opts(legend.position="bottom",legend.direction="horizontal")