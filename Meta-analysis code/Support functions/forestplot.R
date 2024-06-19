##########################################
#
# forest.plot function for meta-analyis with forest plots
#
# Julian Higgins, 1999
# revised for R and extensively re-written, Mark Simmonds 2007
#
##########################################

forest.plot <-

  function(data, method = "random", pooling = "everything", spread.type = "var", title = " ", title1 = " ", title2 = "", xlab = "Effect size", my.pooled.result = c(0, 0, 0), my.pooled.subgroup.results = 0, xlow = 0, xhigh = 0, xinc = 0, label.data = 0, labels = c("Study"), labels1 = c("Study"), labels2 = c(""), label.pos = -10, label.just = 0, label.cex = 1, font.size = rep(1,15), axes.plot = rep(1,4), left.arrow = "Treatment better", right.arrow = "Treatment worse", interval.type = "confidence", pooled.block.size = 1, transform = 1, subgroups = 0, subgroup.title.pos = 999, to.plot = 1, study.block.size = 1, fill.colour = NA, study.blocks = "within", pooled.label.text = c("Pooled"), subgroup.pooled.label.text = c("Pooled"), studies.text = "Effect (95% CI)", show.ests=T, axis.type = "", tick.labels.pos = 999, picture.shift = 0, pooled.result.just = 0.05, confidence = 95, out.file = "", device="")

  {

       
#read in data - y = treatment effect

    data <- as.matrix(data)
    y <- c(data[, 1]) * transform
    Nstudies <- length(y)
	
#construct cis from input data
    
    if(spread.type == "ci" || spread.type == "c.i.") {
      my.ci.low <- c(data[, 2]) * transform
      my.ci.high <- c(data[, 3]) * transform
      confval <- 2 * qnorm(1 - ( (1 - (confidence/100)) / 2 ) )
      spread <- (my.ci.high - my.ci.low) / confval
      w <- 1/(transform * transform * spread * spread)
      s.e <- spread * transform
    }
    
    else if(spread.type == "ciw" || spread.type == "c.i.w") {
      my.ci.low <- c(data[, 2]) * transform
      my.ci.high <- c(data[, 3]) * transform
      spread <- data[, 4]
      w <- spread/(transform * transform)
      s.e <- transform/sqrt(spread) 
    }
  
# or read in ses weights or vars

    else{

      spread <- c(data[, 2])

      if(spread.type == "w") {
        w <- spread/(transform * transform)
        s.e <- ifelse(spread == 0, 0, transform/sqrt(spread))
      }

      else if(spread.type == "var") {
        s.e <- transform * sqrt(spread)
        w <- 1/(transform * transform * spread)
      }

      else if(spread.type == "s.e." || spread.type == "se") {
        w <- 1/(transform * transform * spread * spread)
        s.e <- spread * transform
      }
    }

# user-defined pooling

    if(method == "mine") {
      if(my.pooled.result[2] == 0 && my.pooled.result[3] == 0 &&
         length(my.pooled.subgroup.results) == 1)
        pooling <- "none"
      else if(my.pooled.result[2] == 0 && my.pooled.result[3] == 0 && length(my.pooled.subgroup.results) > 1)
        pooling <- "subgroups"
      else if((my.pooled.result[2] == 0 || my.pooled.result[3] == 0) && length(my.pooled.subgroup.results) == 1)
        pooling <- "overall"
    }

# count null trials   

    Nzerozero <- 0
    for(i in 1:Nstudies) {
      if(y[i] == 0 && spread[i] == 0) {
        Nzerozero <- Nzerozero + 1
        w[i] <- 0
      }
    }

    # confidence intervals
    
    conf.mult <- abs(qnorm(0.5 * (1 - confidence/100)))
    lower.ci <- upper.ci <- NULL
    
    for(i in 1:Nstudies) {
      lower.ci[i] <- min(y[i] - (conf.mult * s.e[i]), y[i] + (conf.mult * s.e[i]))
      upper.ci[i] <- max(y[i] + (conf.mult * s.e[i]), y[i] - (conf.mult * s.e[i]))

      if(spread.type == "c.i." || spread.type == "ci" || spread.type == "ciw" || spread.type == "c.i.w") {
        lower.ci[i] <- my.ci.low[i]
        upper.ci[i] <- my.ci.high[i]
      }
    }


          
# perform meta-analysis
# use external function "meta.analysis" 

    meta.res <- meta.analysis(y=y, w=w, spread=spread, Nstudies=Nstudies, Nzerozero=Nzerozero, method=method, pooling=pooling, my.pooled.result=my.pooled.result, confidence=confidence, subgroups=subgroups, study.blocks=study.blocks, my.pooled.subgroup.results=my.pooled.subgroup.results)

    pooled.est <- meta.res$pooled.est
    pooled.s.e <- meta.res$pooled.s.e
    pooled.lower.ci <- meta.res$pooled.lower.ci
    pooled.upper.ci <- meta.res$pooled.upper.ci

    sg.pooled.est <- meta.res$sg.pooled.est
    sg.pooled.s.e <- meta.res$sg.pooled.s.e
    sg.pooled.lower.ci <- meta.res$sg.pooled.lower.ci
    sg.pooled.upper.ci <- meta.res$sg.pooled.upper.ci

    wstar <- meta.res$wstar
    sum.wstar <- meta.res$sum.wstar
    sg.wstar <- meta.res$sg.wstar
    sg.sum.wstar <- meta.res$sg.sum.wstar

    
####################################
#    
# start plotting set-up 
#
####################################

    
# plot parameters - xpd=NA allows text overflow

    #par(fig = c(0.01, 0.99, 0.01, 0.99), plt = c(0.35 + picture.shift, 0.8 + picture.shift, 0.15, 0.85), las=1, xpd=NA)

    par(fig = c(0.01, 0.99, 0.01, 0.99), plt = c(0.4 + picture.shift, 0.75 + picture.shift, 0.15, 0.85), las=1, xpd=NA)
    
    if(title != " ")
      title1 <- title

# fonts

    roman <- 1
    bold <- 2
    italic <- 3

## font adjusting and number of trials to plot

    if(length(to.plot) != Nstudies) {
      to.plot <- rep(1, Nstudies)
    }
    if(length(font.size) < 15) {
      font.size <- c(font.size, rep(1, 15 - length(font.size)))
    }
    # supress study effect estimates
    if(show.ests==F){
      font.size[12] <- 0
    }

## labels set-up

    labels.text <- as.matrix(format(label.data))
    
    label.font.size <- 1
    if(Nstudies > 20)
      label.font.size <- 1.4 - 0.02 * Nstudies
    
    if(length(labels1) != 1 || labels1[1] != "Study")
      labels <- labels1

    Nlabels <- min(length(labels), length(labels.text[1,]))
    
    labels.plot <- 0
    if ( length(labels.text)>1 )
      labels.plot <- 1

    if(length(labels2) < Nlabels) {
      lenlab2 <- length(labels2)
      labels2[(lenlab2 + 1):Nlabels] <- rep("", (Nlabels - length(labels2)))
    }
    else if(length(labels2) > Nlabels)
      labels2 <- rep("", Nlabels)

# label location

    if(length(label.just) != Nlabels)
      #label.just <- c(0, rep(0.5, Nlabels - 1))
      label.just <- c(0, rep(0, Nlabels - 1))
    if(length(label.pos) != Nlabels && Nlabels == 1)
      label.pos <- c(-8)
    else if(length(label.pos) != Nlabels && Nlabels == 2)
      label.pos <- c(-9, -2)
    else if(length(label.pos) != Nlabels && Nlabels == 3)
      label.pos <- c(-14, -8, -2)
    else if(length(label.pos) != Nlabels && Nlabels == 4)          
      label.pos <- c(-14, -8, -2, 30)

# pooled result and subgroup text

    if(length(pooled.label.text) != Nlabels && Nlabels > 1)
      pooled.label.text[2:Nlabels] <- rep("", (Nlabels - 1))
    else if(length(pooled.label.text) != Nlabels && Nlabels == 1)
      pooled.label.text <- pooled.label.text[1]
        
    if(subgroup.title.pos == 999)
      subgroup.title.pos <- label.pos[1]
    
    if(length(pooled.label.text) == 1 && pooled.label.text[1] == "Pooled" && length(subgroup.pooled.label.text) == 1 & subgroup.pooled.label.text[1] == "Pooled" && pooling == "everything" && length(subgroups) != 1) {
      pooled.label.text <- c("Overall pooled")
    }

    
##  set up plot coords and plot size

    base.y.coord <- 0
    sg.title.y.coord <- NULL
    sg.pooled.y.coord <- NULL
    addpooled <- rep(0, Nstudies)
    addsubgpooled <- rep(0, Nstudies)
    addtitles <- rep(0, Nstudies)
        
    if(pooling == "overall" || pooling == "everything") {
      addpooled <- rep(2, Nstudies)
    }

# set up plot size and labels for subgroups

    Nsubgps <- 1

    if(length(subgroups) != 1) {
      Nsubgps <- length(subgroups)/2
      Ninsubgp <- rep(0, Nsubgps)
      for(s in 1:Nsubgps) {
        Ninsubgp[s] <- as.integer(subgroups[2 * s])
      }
      if(length(subgroup.pooled.label.text) == 1) {
        subgroup.pooled.label.text <- rep("Pooled", Nsubgps)
      }  

      Nsubgp.pooled.labels <- as.integer(length(subgroup.pooled.label.text)/Nsubgps)
      subgroup.pooled.labels <- matrix("o", Nsubgps,Nsubgp.pooled.labels)

      for(s in 1:Nsubgps) {
        subgroup.pooled.labels[s,  ] <- subgroup.pooled.label.text[(1 + Nsubgp.pooled.labels * (s - 1)):(Nsubgp.pooled.labels + Nsubgp.pooled.labels * (s - 1))]
      }

      if(pooling == "subgroups" || pooling == "everything") {
        addsubgpooled <- rep(Nsubgps:1, times = Ninsubgp)
      }
      addtitles <- 1.3 * (rep(Nsubgps:1, times = Ninsubgp) - 1)
    }


# y axis and coordinates 

    y.coord <- c(Nstudies:1) + addpooled + addtitles + addsubgpooled
    top.y.coord <- y.coord[1] + 0.15
    if(length(subgroups) != 1) {
      for(s in (1:Nsubgps)) {
        sg.pooled.y.coord[s] <- y.coord[sum(as.integer(Ninsubgp[1:s]))] - 1
        sg.title.y.coord[s] <- y.coord[sum(as.integer(Ninsubgp[1:s]))] + as.integer(Ninsubgp[s])
      }
      top.y.coord <- sg.title.y.coord[1]
    } 

    if(pooling == "overall" || pooling == "everything") {
      pooled.y.coord <- 1.2
    }

        
# define horizontal axes 

    if(xlow == xhigh) {
      if(axis.type == "log") {
        xlow <- if(exp(min(lower.ci)) > 0.5) 0.5 else 0.25
        xhigh <- ceiling(exp(max(upper.ci)))
        xinc <- 0.25
      }
      else {
        xlow <- floor(min(lower.ci))
        xhigh <- ceiling(max(upper.ci))
        if(xlow > 0)
          xlow <- 0
        if(xhigh < 0)
          xhigh <- 0
        xdist <- xhigh - xlow
        xinc <- {
          if(xdist == 1)
            0.2
          else if(xdist == 2 | xdist == 3)
            0.5
          else if(xdist >= 4 & xdist <= 6)
            1
          else if(xdist >= 7 & xdist <= 15)
            2
          else 5
        }
      }
    }

    
# results text for standard plots

    result.text <- paste(format(round(c(pooled.est, 1.111), 2))[1], "  (", format(round(c(pooled.lower.ci, 1.111), 2))[1], ",", format(round(c(pooled.upper.ci, 1.111), 2))[1], ")")

    all.studies.text <- NULL
    for (t in 1:Nstudies){
      all.studies.text[t] <- paste(format(round(c(y[t], 1.111), 2))[1], "  (", format(round(c(lower.ci[t], 1.111), 2))[1], ",", format(round(c(upper.ci[t], 1.111), 2))[1], ")")
     }

    sg.result.text <- NULL
    if(length(subgroups) != 1) {
      for(s in 1:Nsubgps) {
        sg.result.text[s] <- paste(format(round(c(sg.pooled.est[s], 1.111), 2))[1], "  (",format(round(c(sg.pooled.lower.ci[s], 1.111),2))[1], ",", format(round(c(sg.pooled.upper.ci[s], 1.111), 2))[1], ")")
      }
    }
        
# adjust m-a results for logarithmic plots
        
    vline <- 0
    axistype <- "" 

    if(axis.type == "log") {
      
      axistype <- "x"
      y <- exp(y)
      lower.ci <- exp(lower.ci)
      upper.ci <- exp(upper.ci)
      pooled.est <- exp(pooled.est)
      pooled.lower.ci <- exp(pooled.lower.ci)
      pooled.upper.ci <- exp(pooled.upper.ci)
      vline <- 1
      result.text <- paste(format(round(c(pooled.est, 1.111), 2))[1], "  (", format(round(c(pooled.lower.ci, 1.111), 2))[1],",", format(round(c(pooled.upper.ci, 1.111), 2))[1],")")

      for (t in 1:Nstudies){
        all.studies.text[t] <- paste(format(round(c(y[t], 1.111), 2))[1], "  (", format(round(c(lower.ci[t], 1.111), 2))[1], ",", format(round(c(upper.ci[t], 1.111), 2))[1], ")")
      }      

      if(length(subgroups) != 1) {
        for(s in 1:Nsubgps) {
          sg.pooled.est[s] <- exp(sg.pooled.est[s])
          sg.pooled.lower.ci[s] <- exp(sg.pooled.lower.ci[s])
          sg.pooled.upper.ci[s] <- exp(sg.pooled.upper.ci[s])
          sg.result.text[s] <- paste(format(round(c(sg.pooled.est[s], 1.111), 2))[1], "  (",format(round(c(sg.pooled.lower.ci[s], 1.111),2))[1], ",", format(round(c(sg.pooled.upper.ci[s], 1.111), 2))[1], ")")
        }
      }
    }


    

#####################################
#
# start forest plot
#
#####################################
        

## ticks and titles

    frame()
    if(tick.labels.pos[1] == 999)
      labels.at <- seq(xlow, xhigh, by = xinc)
    else if(tick.labels.pos[1] != 999)
      labels.at <- tick.labels.pos
    titlelabs.extra <- 1.3
    if(length(labels2) != 1 | labels2[1] != "") {
      titlelabs.extra <- 2.3
	}

## define plot area - plot axes
       
    plot(1, 1, type = "n", axes = F, log = axistype, ylab = "", xlab = "", xlim = c(xlow, xhigh), ylim = c(base.y.coord, top.y.coord + titlelabs.extra), col = 1)

    if(axes.plot[1] != 0)
      axis(2, pos = vline, labels = F, tick = F, col = 1)
    if(axes.plot[2] != 0) {
      axis(1, cex.axis = 1 * axes.plot[3], at = labels.at, font = roman, col = 1, pos = 0)
    }

    if(axes.plot[1] != 0) {
      lines(x = c(vline, vline), y = c(base.y.coord, top.y.coord))
    }

## study blocks and ci bars
        
    for(i in 1:Nstudies) {

      if(study.blocks == "equal") {
        wstar[i] <- sum.wstar[i]/5
      }

      if(to.plot[i] != 0) {
        if(spread[i] != 0) {

          if(y[i] >= xlow && y[i] <= xhigh)
            points(y[i], y.coord[i], pch = 15, cex = 4 * study.block.size * sqrt(wstar[i]/sum.wstar[i]), col = 1)

          if(lower.ci[i] != "NA") {
            if(upper.ci[i] <= xhigh && lower.ci[i] >= xlow) {
              lines(c(upper.ci[i], y[i], lower.ci[i]), c(y.coord[i], y.coord[i], y.coord[i]), col = 1)
            }

            if(upper.ci[i] > xhigh && lower.ci[i] >= xlow) {
              lines(c(xhigh, lower.ci[i]), c(y.coord[i], y.coord[i]), col = 1)
              arrows(xhigh - 0.05 * (xhigh - xlow), y.coord[i], xhigh, y.coord[i], length=0.1, col = 1)
            }

            if(lower.ci[i] < xlow && upper.ci[i] <= xhigh) {
              lines(c(upper.ci[i], xlow), c(y.coord[i], y.coord[i]), col = 1)
              arrows(xlow + 0.05 * (xhigh - xlow), y.coord[i], xlow, y.coord[i], length=0.1, col = 1)
            }

            if(lower.ci[i] < xlow && upper.ci[i] > xhigh) {

              lines(c(xhigh, xlow), c(y.coord[i], y.coord[i]), col = 1)
              arrows(xhigh - 0.05 * (xhigh - xlow),y.coord[i], xhigh, y.coord[i], length=0.1, col = 1)
              arrows(xlow + 0.05 * (xhigh - xlow), y.coord[i], xlow, y.coord[i], length=0.1, col = 1)

            }
          }
        }
      }  

      if(axis.type == "log" && y[i] == 1 && spread[i] == 0 && to.plot[i] == 1) {
        points(1, y.coord[i], pch = 4, cex = 1, col = 1)
      }
      else if(y[i] == 0 && spread[i] == 0 && to.plot[i] == 1) {
        points(0, y.coord[i], pch = 4, cex = 1, col = 1)
      }
      
    }

        
## pooled result diamond and text 

    if(pooling == "overall" | pooling == "everything") {

      if(Nstudies != Nzerozero) {    

        polygon(x=c(pooled.upper.ci, pooled.est, pooled.lower.ci, pooled.est, pooled.upper.ci), y=c(pooled.y.coord, pooled.y.coord - 0.5 * pooled.block.size, pooled.y.coord, pooled.y.coord + 0.5 * pooled.block.size, pooled.y.coord), col = fill.colour)    

        if(font.size[11] != 0) {
          text(xhigh, pooled.y.coord, result.text, cex = 1 * font.size[11], adj = pooled.result.just, font = bold, col = 1, las=1)             
        }

      }

      else {
        if(font.size[11] != 0) {
          text(xhigh, pooled.y.coord, "Not estimable",cex = 1 * font.size[11], adj = pooled.result.just, font = roman, col = 1, las=1)
        }
      }
    }

# dashed line at pooled estimate

    if (axes.plot[4] != 0){
      lines(x=rep(pooled.est,2),y=c(base.y.coord,top.y.coord),lty="dashed")
    }

# subgroup results diamonds and text 

    if(length(subgroups) != 1) {
      
      if(pooling == "subgroups" || pooling == "everything") {
        
        for(s in 1:Nsubgps) {
          
          if(Ninsubgp[s] != meta.res$sg.Nzerozero[s]) {
            #lines(c(sg.pooled.upper.ci[s], sg.pooled.est[s], sg.pooled.lower.ci[s], sg.pooled.est[s],sg.pooled.upper.ci[s]), c(sg.pooled.y.coord[s], sg.pooled.y.coord[s] - 0.5 * pooled.block.size, sg.pooled.y.coord[s], sg.pooled.y.coord[s] + 0.5 * pooled.block.size, sg.pooled.y.coord[s]), col = 1)  
            polygon(x=c(sg.pooled.upper.ci[s], sg.pooled.est[s], sg.pooled.lower.ci[s], sg.pooled.est[s],sg.pooled.upper.ci[s]), y=c(sg.pooled.y.coord[s], sg.pooled.y.coord[s] - 0.5 * pooled.block.size, sg.pooled.y.coord[s], sg.pooled.y.coord[s] + 0.5 * pooled.block.size, sg.pooled.y.coord[s]), col = fill.colour)  
            
            if(font.size[15] != 0) {
              text(xhigh, sg.pooled.y.coord[s], sg.result.text[s], cex = 1 * font.size[15], adj = pooled.result.just, font = roman, col = 1, las=1)
            }
          }

          else {
            if(font.size[15] != 0) {
              text(xhigh, sg.pooled.y.coord[s], "Not estimable", cex = 1 * font.size[15], adj = pooled.result.just, font = roman, col = 1, las=1)
            }
          }

        # subgroup pooled results labels

          if(font.size[14] != 0) {
            for(i in 1:Nsubgp.pooled.labels) {
              mtext(subgroup.pooled.labels[s, i], side = 2, at = sg.pooled.y.coord[s], line =  - label.pos[i], cex = label.cex * font.size[14], adj = label.just[i], font = italic, col = 1, las=1)
            }
          }
        }
      }
    }
    
# text for each study result

    if(font.size[12] != 0){
      for (t in 1:Nstudies){
        text(xhigh, y.coord[t], all.studies.text[t], cex = 0.9 * font.size[12], adj = pooled.result.just, font = roman, col = 1, las=1)
      }
    }

## label text 

    for(i in 1:Nlabels) { 

      first.coord <- top.y.coord + 1.3 + 0.04 * (top.y.coord - base.y.coord)

      second.coord <- max(top.y.coord + 1.3 - 0.03 * (top.y.coord - base.y.coord), top.y.coord + 1.3 + 0.04 * (top.y.coord - base.y.coord) - 1.1)

# label headings  

      if(font.size[4] != 0 & labels.plot==1) {
        if(labels2[i] == "")
          mtext(labels[i], side = 2, at = (first.coord + second.coord)/2, line =  - label.pos[i], cex= 1 * font.size[4], adj = label.just[i],font = bold, col = 1, las=1)

        else if(length(labels2) != 1 || labels2[i] != "") {
          mtext(labels[i], side = 2, at = first.coord, line =  - label.pos[i], cex = label.cex * font.size[4], adj = label.just[i], font = bold, col = 1, las=1)
          mtext(labels2[i], side = 2, at = second.coord, line =  - label.pos[i], cex = label.cex * font.size[4], adj = label.just[i], font = bold, col = 1, las=1)
        }
      }

# study labels

      if(font.size[6] != 0 & labels.plot==1) {
        mtext(labels.text[, i], side = 2, at = y.coord, line =  -label.pos[i], cex = label.cex * font.size[6] * label.font.size, adj = label.just[i], font = roman, col = 1, las=1)
      }

# pooled result text 

      if((pooling == "overall" | pooling == "everything") && font.size[8] != 0)
        mtext(pooled.label.text[i], side = 2, at = pooled.y.coord, line =  - label.pos[i], cex = label.cex * font.size[8], adj = label.just[i], font = italic, col = 1, las=1)
      
    }
    
# subgroup headings
        
    if(length(subgroups) != 1 && font.size[13] != 0) {
      for(s in 1:Nsubgps) {
        mtext(subgroups[s * 2 - 1], side = 2, at = sg.title.y.coord[s], line =  - subgroup.title.pos, cex = label.cex[1] * font.size[13], adj = 0, font = italic, col= 1, las=1)
      }
    }

# study estimates label

    if(font.size[5] != 0 & font.size[12] != 0){
      text(xhigh, (first.coord + second.coord)/2, studies.text, cex = label.cex * font.size[5], adj = pooled.result.just, font = bold, col = 1, las=1)
    }

# end labels

## titles 

    title1.line <- 3
    title2.line <- 2
    if(title2 == "" || title2 == " ") {
      title1.line <- 2
    }

## titles for log axes
    
    if(axis.type == "log") {
      centre <- exp(log(xlow) - ((0.35 + picture.shift) * (log(xhigh) - log(xlow)))/0.45 + (log(xhigh) - log(xlow))/0.9)

      if(font.size[1] != 0) {
        mtext(title1, side = 3, line = title1.line, at = centre, cex = 1.5 * font.size[1], font = bold, col = 1, las=1)
      }

      if(font.size[2] != 0) {
        mtext(title2, side = 3, line = title2.line, at = centre,cex = 1.5 * font.size[2], font = bold, col = 1, las=1)
      }

      # effect arrows and text
      
      if(font.size[9] != 0 & font.size[10] != 0) {
        arrows(exp(-0.02 * (log(xhigh) - log(xlow))), base.y.coord - 0.25 * (top.y.coord - base.y.coord), exp(-0.1 * (log(xhigh) - log(xlow))), base.y.coord - 0.25 * (top.y.coord - base.y.coord),length=0.1, col = 1)
      }
      if (font.size[9] != 0){
        text(exp(-0.12 * (log(xhigh) - log(xlow))), base.y.coord - 0.25 * (top.y.coord - base.y.coord), left.arrow, cex = 1.05 * font.size[9], adj = 1, font = roman, col = 1, las=1)
      }

      if(font.size[9] != 0 & font.size[10] != 0) {
        arrows(exp(0.02 * (log(xhigh) - log(xlow))), base.y.coord - 0.25 * (top.y.coord - base.y.coord), exp(0.1 * (log(xhigh) - log(xlow))), base.y.coord - 0.25 * (top.y.coord - base.y.coord),length=0.1, col = 1)
      }
      if (font.size[9] != 0){
        text(exp(0.12 * (log(xhigh) - log(xlow))), base.y.coord - 0.25 * (top.y.coord - base.y.coord), right.arrow, cex = 1.05 * font.size[9], adj= 0, font = roman, col = 1, las=1)
      }
    }

## titles for standard axes 

    else if(axis.type != "log") {
      centre <- xlow - ((0.35 + picture.shift) * (xhigh - xlow))/0.45 +(xhigh - xlow)/0.9

      if(font.size[1] != 0) {
        mtext(title1, side = 3, line = title1.line, at = centre,cex = 1.5 * font.size[1], font = bold, col = 1, las=1)
      }

      if(font.size[2] != 0 & title2 != "") {
        mtext(title2, side = 3, line = title2.line, at = centre,cex = 1.2 * font.size[2], font = bold, col= 1, las=1)
      }

# effect arrows and text

      if(font.size[9] != 0 & font.size[10] != 0) {
        arrows(vline - 0.02 * (xhigh - xlow), base.y.coord - 0.22 * (top.y.coord - base.y.coord), vline - 0.1 * (xhigh - xlow), base.y.coord - 0.22 * (top.y.coord - base.y.coord),length=0.1, col = 1)
      }
      if(font.size[9] != 0){
        text(vline - 0.125 * (xhigh - xlow), base.y.coord - 0.22 *(top.y.coord - base.y.coord), left.arrow, cex = 1 * font.size[9], adj = 1, font = italic, col = 1, las=1)
      }
    
      if(font.size[9] != 0 & font.size[10] != 0) {
        arrows(vline + 0.02 * (xhigh - xlow), base.y.coord - 0.22 * (top.y.coord - base.y.coord), vline + 0.1 * (xhigh - xlow), base.y.coord - 0.22 * (top.y.coord - base.y.coord),length=0.1, col = 1)
      }
      if(font.size[9] != 0){
        text(vline + 0.125 * (xhigh - xlow), base.y.coord - 0.22 *(top.y.coord - base.y.coord), right.arrow, cex = 1 * font.size[9], adj = 0, font = italic, col = 1, las=1)       
      }
    }
        

## remaining text

    if(font.size[3] != 0) {
      mtext(paste("Estimates with ", as.character(confidence), "% ", interval.type, " intervals", sep = ""), side = 3, line = 0.1, cex = 1.1 * font.size[3], font = italic, col = 1, las=1)
    }

# x axis label
    
    if(font.size[7] != 0) {
      text(vline, base.y.coord - 0.15 * (top.y.coord - base.y.coord), xlab, cex = 1.1 * font.size[7], adj = 0.5, font = italic, col = 1, las=1)
    }

##  output
    
    if(out.file!="" & device=="ps")
		dev.print(postscript, file = out.file, horizo = F, onefile = F, print = F, paper = "a4", width = 12, height = 15)

    if(out.file!="" & device=="pdf")
		dev.print(pdf, file = out.file, onefile = F, paper = "a4", width = 12, height = 15)

    if(out.file!="" & device=="jpeg")
		dev.print(jpeg, file = out.file, quality=100, width=2100, height=2970)

    if(out.file!="" & device=="bmp")
		dev.print(bmp, file = out.file, width=2100, height=2970)

    if(out.file!="" & device=="png")
		dev.print(png, file = out.file, width=2100, height=2970, bg="white")

    if(out.file!="" & device=="tiff")
		dev.print(tiff, file = out.file, res=400, bg="white")


    
# end forest.plot

  }










########################################################
#
# Function to perform the actual meta-analysis when producing a forest plot
#
#######################################################


meta.analysis <-

  function(y,w,spread,Nstudies,Nzerozero=0,method="random",pooling="everything",my.pooled.result=c(0,0,0),confidence=95,subgroups=0,study.blocks="across",my.pooled.subgroup.results=0)
  {

    wstar <- NULL
    conf.mult <- abs(qnorm(0.5 * (1 - confidence/100)))
    muhat <- sum(y * w)/sum(w)
    Q <- sum(w * ((y - muhat)^2))
    tausqhat1 <- (Q - Nstudies + Nzerozero + 1)/(sum(w) - (sum(w^2))/sum(w))
    tausq.hat <- ifelse(Nstudies - Nzerozero <= 1, 0, ifelse(tausqhat1 <= 0,0, tausqhat1))

    if(method != "random")
          wstar <- w
    if(method == "random") {
      for(i in 1:Nstudies) {
        wstar[i] <- 1/((1/w[i]) + tausq.hat)
        if(y[i] == 0 && spread[i] == 0)
          wstar[i] <- 0
      }

    }

    pooled.est <- sum(y * wstar)/sum(wstar)
    pooled.s.e <- sqrt(1/sum(wstar))
    pooled.lower.ci <- min(pooled.est - (conf.mult * pooled.s.e), pooled.est + (conf.mult * pooled.s.e))
    pooled.upper.ci <- max(pooled.est + (conf.mult * pooled.s.e), pooled.est - (conf.mult * pooled.s.e))

    if(method == "mine" && (pooling == "everything" || pooling == "overall")) {
      pooled.est <- my.pooled.result[1]
      pooled.lower.ci <- my.pooled.result[2]
      pooled.upper.ci <- my.pooled.result[3]
    }


    
# analysis of subgroups

    Nsubgps <- 1
    if(length(subgroups) != 1) {
          
          Nsubgps <- length(subgroups)/2
          if(length(my.pooled.subgroup.results) == 1 && method == "mine") {
            my.pooled.subgroup.results <- rep(0, 3 * Nsubgps)
          }
          Ninsubgp <- rep(0, Nsubgps)
          for(s in 1:Nsubgps) {
            Ninsubgp[s] <- as.integer(subgroups[2 * s])
          }
        }

    sum.wstar <- rep(sum(wstar), Nstudies)
    
      start <- 0
      sg.pooled.est <- rep(0, Nsubgps)
      sg.pooled.s.e <- rep(0, Nsubgps)
      sg.pooled.lower.ci <- rep(0, Nsubgps)
      sg.pooled.upper.ci <- rep(0, Nsubgps)
      sg.sum.wstar <- rep(0, Nsubgps)
      sg.wstar <- rep(0, Nsubgps)
      sg.Nzerozero <- rep(0, Nsubgps)

    if(length(subgroups) != 1) {
      
      for(s in 1:Nsubgps) {

        sg.y <- rep(0, Ninsubgp[s])
        sg.w <- rep(0, Ninsubgp[s])
        sg.wstar <- rep(0, Ninsubgp[s])
        sg.spread <- rep(0, Ninsubgp[s])

     for(i in 1:Ninsubgp[s]) {
       sg.y[i] <- y[start + i]
       sg.w[i] <- w[start + i]
       sg.spread[i] <- spread[start + i]
       if(sg.y[i] == 0 && sg.spread[i] == 0) {
         sg.Nzerozero[s] <- sg.Nzerozero[s] + 1
       }
     }

        sg.muhat <- sum(sg.y * sg.w)/sum(sg.w)
        sg.Q <- sum(sg.w * ((sg.y - sg.muhat)^2))
        sg.tausqhat1 <- (sg.Q - Ninsubgp[s] + sg.Nzerozero[s] + 1)/(sum(sg.w) - (sum(sg.w^2))/sum(sg.w))
        sg.tausq.hat <- ifelse(Ninsubgp[s] - sg.Nzerozero[s] <= 1, 0, ifelse(sg.tausqhat1 <= 0, 0, sg.tausqhat1))


        if(method != "random" || Ninsubgp[s] - sg.Nzerozero[s] == 1)
          sg.wstar <- sg.w
        if(method == "random" && Ninsubgp[s] - sg.Nzerozero[s] > 1) {
          for(i in 1:Ninsubgp[s]) {
            sg.wstar[i] <- 1/((1/sg.w[i]) + sg.tausq.hat)
            if(sg.y[i] == 0 && sg.spread[i] == 0)
              sg.wstar[i] <- 0
          }
        }


        sg.pooled.est[s] <- sum(sg.y * sg.wstar)/sum(sg.wstar)
        sg.pooled.s.e[s] <- sqrt(1/sum(sg.wstar))
        sg.pooled.lower.ci[s] <- min(sg.pooled.est[s] - (conf.mult * sg.pooled.s.e[s]), sg.pooled.est[s] +(conf.mult * sg.pooled.s.e[s]))
        sg.pooled.upper.ci[s] <- max(sg.pooled.est[s] + (conf.mult * sg.pooled.s.e[s]), sg.pooled.est[s] -(conf.mult * sg.pooled.s.e[s]))

        if(method == "mine") {
          sg.pooled.est[s] <- my.pooled.subgroup.results[((3 * s) - 2)]
          sg.pooled.lower.ci[s] <- my.pooled.subgroup.results[(3 * s) - 1]
          sg.pooled.upper.ci[s] <- my.pooled.subgroup.results[(3 * s)]
        }

        sg.sum.wstar[s] <- sum(sg.wstar)

        
        if((pooling == "subgroups" || pooling == "none") && study.blocks != "across") {

          for(i in 1:Ninsubgp[s]) {
            wstar[start + i] <- sg.wstar[i]
          }
        }


        start <- start + Ninsubgp[s]

      }


      if(pooling == "subgroups" && study.blocks != "across") {
            sum.wstar <- rep(sg.sum.wstar[1:Nsubgps], times = Ninsubgp)
          }
    }

    
# end subgroups section
        

        
# return results

    out <- list(pooled.est=pooled.est, pooled.s.e=pooled.s.e, pooled.lower.ci=pooled.lower.ci, pooled.upper.ci=pooled.upper.ci, sg.pooled.est=sg.pooled.est, sg.pooled.s.e=sg.pooled.s.e, sg.pooled.lower.ci=sg.pooled.lower.ci, sg.pooled.upper.ci=sg.pooled.upper.ci, wstar=wstar, sum.wstar=sum.wstar, sg.wstar=sg.wstar, sg.sum.wstar=sg.sum.wstar,sg.Nzerozero=sg.Nzerozero)

    return(out)
    
  }
