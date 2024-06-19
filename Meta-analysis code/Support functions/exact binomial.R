# exact binomial conf interval

exact.binom = function(events,total,ci=0.95){
  
  p.event = events / total
  prob.set = seq(0,1,0.005)  
  co.low = (1-ci)/2
  co.high = 1-co.low
  q.low = qbinom(co.low,total,prob.set)
  q.high = qbinom(co.high,total,prob.set)
  prob.low = rev(prob.set[q.high==events])[1]
  prob.high = rev(prob.set[q.low==events])[1]  
  out = data.frame(events,total,p.event,prob.low,prob.high)
  return(out)
}


exact.binom2 = function(events,total,ci=0.95){
  
  p.event = events / total
  prob.set = pbinom(1:events,total,p.event)
  co.low = (1-ci)/2
  nmin = length(prob.set[prob.set<=co.low])
  co.high = 1-co.low
  nmax = total - length(prob.set[prob.set>=co.high])
  out = data.frame(events,nmin,nmax,total,p.event,p.low=nmin/total,p.high=nmax/total)
  return(out)
}