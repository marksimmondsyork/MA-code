# convert log effect to effect and  95% CI

logconv = function(TE,SE){
  
  out= data.frame(Effect = exp(TE), CIlow = exp(TE-1.96*SE), CIhigh = exp(TE+1.96*SE))
  return(out)
}