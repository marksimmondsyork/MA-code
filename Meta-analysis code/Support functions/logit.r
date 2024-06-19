# the logit function

logit <-

  function(x){

    return( log(x / (1-x)) )

  }

expit <-

  function(x){

    return( exp(x)/(1+exp(x)) )

  }
