############################
###_METABOLIC RATES FIGS_###
############################

fig2  <-  function(data, nlModel, combs=FALSE) {
 
  xnA    <-  nlModel['A','mean']
  xnEr   <-  nlModel['Er','mean']
  xnboTc <-  nlModel['boTc','mean']
  xnEi   <-  nlModel['Ei','mean']
  xnTopt <-  nlModel['Topt','mean']
  nA     <-  xnA    + nlModel[paste0('r[',data$fam,',1]'),'mean']
  nEr    <-  nlModel['logit.Er','mean'] + nlModel[paste0('r[',data$fam,',2]'),'mean']
  nEr    <-  xnEi/(1 + exp((-1)*nEr))
  nboTc  <-  xnboTc + nlModel[paste0('r[',data$fam,',3]'),'mean']
  nEi    <-  rep(xnEi,   nrow(data))
  nTopt  <-  xnTopt + nlModel[paste0('r[',data$fam,',4]'),'mean']
  
  #corrects metabolic rates for temperature first and then size scaling
  btNlmer  <-  data$rate - nEr/8.62e-5*(1/293 - 1/data$temp) + log(1 + exp(nEi/8.62e-5*(1/nTopt - 1/data$temp))*nEr/(nEi-nEr))
  bwNlmer  <-  data$rate - nA*data$weight
  
  if(!combs) {
    par(mfrow=c(1, 2), family="Times", mai=c(1.02,0.82,0.82,0.60), omi=c(0.5,1,0.5,1))
    ylab1  <-  expression(paste("ln(Rate @ 20"*degree,"C) (g C d"^{-1}, ")"),sep="")
    ylab2  <-  expression(paste("ln(Rate @ 1g) (g C d"^{-1},")"),sep="")
    leg1   <-  '(a)'
    leg2   <-  '(b)'
  } else {
    ylab1  <-  ylab2  <-  ''
    leg1   <-  '(c)'
    leg2   <-  '(d)'
  }
  
  chosen1  <-  data$reef=="yes"
  chosen2  <-  data$reef=="no"
  
  #(a)
  plot(btNlmer[chosen2] ~ data$weight[chosen2], pch=21, col="grey70", bg=make.transparent("grey70", .5),  xlab="ln(Mass) (g)", ylab=ylab1, las=1, xlim=c(-5,10), ylim=c(-10,8), cex.lab=1.3, xpd=NA)
  
  for(j in unique(data$fam)) {
    xpoints  <-  range(data$weight[data$fam==j])
    expr     <-  (xnboTc + nlModel[paste0('r[',j,',3]'),'mean']) + 
                 (xnA    + nlModel[paste0('r[',j,',1]'),'mean'])*xpoints
    points(xpoints, expr, type="l", lty=2, col="grey30", lwd=1)
  }
  points(c(-8, 12), xnboTc + xnA*c(-8, 12), type="l", lty=1, col="black", lwd=2.5)
  points(btNlmer[chosen1] ~ data$weight[chosen1], pch=21, col="dodgerblue4", bg=make.transparent("dodgerblue4", .5))
  #fig position label
  label(px=.03, py=.97, leg1, font=2, cex=1.4)
  #trends labels
  label(px=c(.85,.95), py=c(.95,.95), text=FALSE, type="l", lty=1, lwd=2.5, col="black")
  label(px=c(.85,.95), py=c(.87,.87), text=FALSE, type="l", lty=2, lwd=1, col="grey30")
  label(px=.85, py=.95, lab=substitute("mean trend: "~italic(y) == boTc +  a*italic(x), list(boTc=round(xnboTc,2), a=round(xnA,2))), pos=2, cex=0.9)
  label(px=.85, py=.87, lab="family-level variation", pos=2, cex=0.9)
  
  if(!combs) {
  #color labels
  label(px=.03, py=1.15, text=FALSE, type="o", pch=21, col="dodgerblue4", bg=make.transparent("dodgerblue4", .5), xpd=NA, cex=1.2)
  label(px=.03, py=1.05, text=FALSE, type="o", pch=21, col="grey70", bg=make.transparent("grey70", .5), xpd=NA, cex=1.2)

  #text for color labels
    label(px=.03, py=1.145, lab=paste0("reef fishes; n = ",length(data$reef[chosen1])), pos=4, cex=0.9, xpd=NA)
    label(px=.03, py=1.045, lab=paste0("other fishes; n = ",length(data$reef[chosen2])), pos=4, cex=0.9, xpd=NA)
  }

  #(b)
  plot(1/8.62E-5*(1/293.15 - 1/data$temp)[chosen2], bwNlmer[chosen2], pch=21, col="grey70", bg=make.transparent("grey70", .5),  xlab=expression(paste("Inverse Temperature, 1/", italic(kT[c]), " - 1/", italic(kT), " (1/eV)", sep="")), ylab=ylab2, las=1, ylim=c(-12,2), cex.lab=1.3, xpd=NA)
  
  for(j in unique(data$fam)) {
    xpoints  <-  range(data$temp[data$fam==j])
    xpoints  <-  seq(xpoints[1], xpoints[2], length.out=50)
    nresEr   <-  nlModel['logit.Er','mean'] + nlModel[paste0('r[', j, ',2]'),'mean']
    nresEr   <-  xnEi/(1 + exp((-1)*nresEr))
    nresboTc <-  xnboTc + nlModel[paste0('r[',j,',3]'),'mean']
    nresTopt <-  xnTopt + nlModel[paste0('r[',j,',4]'),'mean']
    expr     <-  nresboTc + nresEr/8.62e-5*(1/293.15-1/xpoints) - log(1 + exp(xnEi/8.62e-5*(1/nresTopt - 1/xpoints))*nresEr/(xnEi-nresEr))
    points(1/8.62e-5*(1/293.15-1/xpoints), expr, type="l", lty=2, col="grey30", lwd=1)
  }
  
  mean.nx          <-  1/8.62e-5*(1/293.15-1/(263:323))
  mean.expr.nlmer  <-  xnboTc + xnEr/8.62e-5*(1/293.15-1/(263:323)) - log(1 + exp(xnEi/8.62e-5*(1/xnTopt - 1/(263:323)))*xnEr/(xnEi-xnEr))
  points(mean.nx, mean.expr.nlmer, type="l", lty=1, col="black", lwd=2.5)
  points(1/8.62E-5*(1/293.15 - 1/data$temp)[chosen1], bwNlmer[chosen1], pch=21, col="dodgerblue4", bg=make.transparent("dodgerblue4", .5))
  label(px=c(.94,.99), py=c(.79,.79), text=FALSE, type="l", lty=1, lwd=2.5, col="black")
  
  label(px=0.1, py=.92, 
    lab=
       substitute(italic(y) == xnboTc + xnEr*italic(x), 
         list(xnboTc=round(xnboTc,2), 
           xnEr=round(xnEr,2) 
         )
       ), 
    pos=4, cex=0.9
  )
 
  label(px=0.1, py=.77, 
    lab=
      substitute(-ln*bgroup('(', 1 + finexpr*'exp'*bgroup('[', expr, ']'), ')'), 
        list(xnEr=round(xnEr,2), 
            xnEi=round(xnEi,2), 
            finexpr=round(xnEr/(xnEi-xnEr),2), 
            expr=substitute(xnEi*bgroup('(',frac(1,italic(k)*xnTopt) - frac(1,italic(k)*italic(T[c])) + italic(x),')'), 
                list(xnEi=round(xnEi,2), 
                  xnTopt=round(xnTopt)
                )
            )
          )
      ),
    pos=4, cex=0.9
  )
  if(!combs) {
    axis(side=3, at=1/8.62e-5*(1/293.15-1/(273.15+seq(0,40,by=5))), labels=seq(0,40,by=5))
    label(px=.5, py=1.25, expression(paste("Temperature ("*degree,"C)",sep="")), xpd=NA, cex=1.3, adj=c(0.5,0.5))
  }
  label(px=.03, py=.97, leg2, font=2, cex=1.4)
}

figS1  <-  function(data, JAGSmodel, nlmerModel) {

  par(mfcol=c(9,5), family="Times", oma=c(4,4,2,2))
  n  <-  0 #start counting
  xnA    <-  JAGSmodel['A','mean']
  xnEr   <-  JAGSmodel['logit.Er','mean']
  xnboTc <-  JAGSmodel['boTc','mean']
  xnEi   <-  JAGSmodel['Ei','mean']
  xnTopt <-  JAGSmodel['Topt','mean']

  nres  <-  coef(nlmerModel)$family

  for(i in sort(unique(data$fam))) {
    n  <-  n + 1
    
    if(n %in% 1:9){par(mar=c(0.1,1.1,0.1,0.1)); yaxt='s'}else{par(mar=c(0.1,1.1,0.1,0.1)); yaxt='n'}
    if(n %in% c(9,18,27,36,45)){xaxt='s'}else{xaxt='n'}
    
    nA     <-  xnA    + JAGSmodel[paste0('r[',i,',1]'),'mean']
    nEr    <-  xnEr   + JAGSmodel[paste0('r[',i,',2]'),'mean']
    nboTc  <-  xnboTc + JAGSmodel[paste0('r[',i,',3]'),'mean']
    nTopt  <-  xnTopt + JAGSmodel[paste0('r[',i,',4]'),'mean']

    x.emp  <- data$temp[data$fam==i]
    x.fit  <- seq(min(x.emp),max(x.emp),length=50)
    y.emp  <- data$rate[data$fam==i] - nA*data$weight[data$fam==i] 
    
    
    JAGSfit <- c(FullMinusEa(
      nboTc,
      nEr,
      xnEi,
      nTopt,
      nA,
      1,
      x.fit))
 
    nlmerFit <- c(FullMinusEa(
      nres[unique(data$family[data$fam==i]),'ln.boTc'],
      nres[unique(data$family[data$fam==i]),'logit.Er'],
      nres[unique(data$family[data$fam==i]),'Ei'],
      nres[unique(data$family[data$fam==i]),'Topt'],
      nres[unique(data$family[data$fam==i]),'a'],
      1,
      x.fit))

    
    plot(c(x.emp, x.fit), c(y.emp, JAGSfit), xlab='', ylab='', type='n', las=1, xlim=c(270,312), ylim=c(-11.5,-2), xaxt=xaxt, yaxt=yaxt,cex.axis=1.1)
    points(x.emp, y.emp, pch=21, col='grey50', bg=make.transparent('grey50', .5), cex=1.5)
    lines(x.fit, JAGSfit, lwd=3)
    lines(x.fit, nlmerFit, lty=2, col='darkorange', lwd=2)
    label(-0.04, 0.88, unique(data$family[data$fam==i]), pos=4, cex=1.2)
    
  }
  mtext("Temperature (K)", side=1, line=2.5, outer=TRUE, cex=1.2)
  mtext(expression(paste("ln(Rate @ 1g) (g C", " d"^{-1}, ")", sep="")), side=2, line=1.7, outer=TRUE, cex=1.4)
}

figS2  <-  function(JAGSmatrix) {
  par(mfrow=c(5,5), family="Times", omi=rep(1,4), mai=rep(0,4), cex=1)
  subIter  <-  JAGSmatrix[,c('boTc', 'A', 'Er', 'Ei', 'Topt')]
  ylabs    <-  list(substitute(a *~(b), list(a=substitute(ln(bar(italic(b[o](T[c]))))), b=substitute(fir^{-alpha} * sec^{-n2}, list(fir="g C g", sec=" d", n2=1)))),
                    substitute(italic(alpha)),
                    substitute(italic(E[r])*~(eV)),
                    substitute(italic(E[i])*~(eV)),
                    substitute(italic(T[opt])*~(K)))
  val  <-  1:ncol(subIter)
  for(j in val) {
    for(i in val) {
      if(i==j) {
        par(mai=c(0.2,0.35,0.3,0.2))
        a <- hist(subIter[,i], main='', xlab='', ylab='', xaxt='n', yaxt='n', pch=21, col='dodgerblue', breaks=13)
        label(c(0, 1), c(0, 0), text=FALSE, type='l', xpd=NA)
        label(c(0, 0.5, 1), rep(-.1, 3), a$breaks[c(1,7,13)], xpd=NA, cex=0.8, adj=c(.5,.5))
        label(c(0, 0), c(0, 1), text=FALSE, type='l', xpd=NA)
        b  <-  range(a$density)
        b  <-  seq(b[1], b[2], length.out=3)
        label(rep(-.2, 3), c(0, 0.5, 1), round(b,2), xpd=NA, cex=0.8, adj=c(.5,.5))
        box("figure")
        label(0.5, 1.2, ylabs[[i]], cex=0.7, adj=c(0.5,0.5), xpd=NA)
      } else {
        if(j!=1 & i==1){yas='s'}else{yas='n'}
        if(j==5){xas='s'}else{xas='n'}

        par(mai=rep(0,4))
        rangeX  <-  max(subIter[,j]) - min(subIter[,j])
        rangeX  <-  c(min(subIter[,j]) - rangeX*0.1, max(subIter[,j]) + rangeX*0.1)
        rangeY  <-  max(subIter[,i]) - min(subIter[,i])
        rangeY  <-  c(min(subIter[,i]) - rangeY*0.1, max(subIter[,i]) + rangeY*0.1)
        plot(subIter[,i], subIter[,j], xlab='', ylab='', pch=21, col='dodgerblue', bg=make.transparent('dodgerblue', 0.5), xlim=rangeY, ylim=rangeX, axes=FALSE)
        box()
        if(yas=='s')axis(2, at=round(seq(min(subIter[,j]), max(subIter[,j]), length.out=4), 2), las=1, cex.axis=0.8)
        if(xas=='s')axis(1, at=round(seq(min(subIter[,i]), max(subIter[,i]), length.out=4), 2), cex.axis=0.8)
        if(j==1 & i==5){
          axis(3, at=round(seq(min(subIter[,i]), max(subIter[,i]), length.out=4), 2), cex.axis=0.8, xpd=NA)
          axis(4, at=round(seq(min(subIter[,j]), max(subIter[,j]), length.out=4), 2), las=1, cex.axis=0.8, xpd=NA)
        }
      }
    }
  }
}

figS3  <-  function(JAGSmatrix) {
  par(mfrow=c(3,2), family="Times", omi=c(1,1,0.5,1), cex=1)
  
  subIter  <-  JAGSmatrix[,c('boTc', 'A', 'logit.Er', 'Er', 'Ei', 'Topt')]
  ylabs    <-  list(substitute(a *~(b), list(a=substitute(ln(bar(italic(b[o](T[c]))))), b=substitute(fir^{-alpha} * sec^{-n2}, list(fir="g C g", sec=" d", n2=1)))),
                    substitute(italic(alpha)),
                    substitute(italic(E[r]), list(E="E'")),
                    substitute(italic(E[r])*~(eV)),
                    substitute(italic(E[i])*~(eV)),
                    substitute(italic(T[opt])*~(K)))
  chains   <-  nrow(subIter) / 3
  for(i in 1:ncol(subIter)) {
    if(i %in% c(5,6)){xas='s'}else{xas='n'}
    if(i %in% c(1,3,5)){par(mai=c(0.1, 0.2, 0, 0.45))}else{par(mai=c(0.1, 0.5, 0, 0.15))}
    rangeVals  <-  max(subIter[,i]) - min(subIter[,i])
    rangeVals  <-  c(min(subIter[,i]) - rangeVals*0.05, max(subIter[,i]) + rangeVals*0.15)
    plot(1:chains, subIter[1:chains,i], xlab='', ylab='', xlim=c(1,max(chains[[1]])), ylim=rangeVals, cex.axis=1, cex.lab=1.6, las=1, type='l', lwd=1, col='dodgerblue', xpd=NA, xaxt=xas)
    label(0.02, 0.9, ylabs[[i]], cex=1, pos=4)
    points(1:chains, subIter[(chains+1):(chains*2),i], type='l', lwd=1, col='firebrick')
    points(1:chains, subIter[((chains*2)+1):(chains*3),i], type='l', lwd=1, col='darkolivegreen')
  }
  mtext('Steps', outer=TRUE, line=-36, cex=1.6)
}

figS4  <-  function(data, nlModel, JAGSModel){
  #creates fam-specific parameters of interest
  par(mfrow=c(2, 2), family="Times", mai=c(0,0.82,1.84,0.60), omi=c(1,1,0,1), cex=1)

  bnot      <-  JAGSModel['boTc','mean']
  alph      <-  JAGSModel['A','mean']
  acen      <-  JAGSModel['Er','mean']
  afam      <-  alph + JAGSModel[paste0('r[',data$fam,',1]'),'mean']
  efam      <-  acen + JAGSModel[paste0('r[',data$fam,',2]'),'mean']
  ifam      <-  bnot + JAGSModel[paste0('r[',data$fam,',3]'),'mean']
  bt        <-  data$rate - efam/8.62E-5*(1/293.15 - 1/data$temp) #corrects ln.b for family-specific E (efam) from the Bayesian output
  bw        <-  data$rate - (afam*data$weight) #Temperature effect on for size-corrected values
  
  
  chosen1  <-  as.character(data$reef)=="yes"
  chosen2  <-  as.character(data$reef)=="no"
  
  #(a)
  plot(bt[chosen2] ~ data$weight[chosen2], pch=21, col="grey70", bg=make.transparent("grey70", .5),  xlab="", ylab="", las=1, xlim=c(-5,10), ylim=c(-10,8), cex.lab=1.3, xaxt='n')
  
  for(j in unique(data$fam)) {
    xpoints  <-  range(data$weight[data$fam==j])
    expr     <-  (bnot + JAGSModel[paste0('r[',j,',3]'),'mean']) + 
                 (alph + JAGSModel[paste0('r[',j,',1]'),'mean'])*xpoints
    points(xpoints, expr, type="l", lty=2, col="grey30", lwd=1)
  }

  points(c(-8, 12), bnot + alph*c(-8, 12), type="l", lty=1, col="black", lwd=2.5)
  points(bt[chosen1] ~ data$weight[chosen1], pch=21, col="dodgerblue4", bg=make.transparent("dodgerblue4", .5))
  #fig position label
  label(px=.03, py=.97, expression(paste(bold("(a)"),sep="")), cex=1.4)
  #trends labels
  label(px=c(.85,.95), py=c(.95,.95), text=FALSE, type="l", lty=1, lwd=2.5, col="black")
  label(px=c(.85,.95), py=c(.87,.87), text=FALSE, type="l", lty=2, lwd=1, col="grey30")
  label(px=.85, py=.95, lab=substitute('mean trend: '~italic(y) == A*italic(x) - B, list(A=round(alph,2), B=-1*round(bnot,2))), pos=2, cex=0.9)
  label(px=.85, py=.87, lab="family-level variation", pos=2, cex=0.9)
  #color labels
  label(px=.03, py=1.15, text=FALSE, type="o", pch=21, col="dodgerblue4", bg=make.transparent("dodgerblue4", .5), xpd=NA, cex=1.2)
  label(px=.03, py=1.05, text=FALSE, type="o", pch=21, col="grey70", bg=make.transparent("grey70", .5), xpd=NA, cex=1.2)
  #text for color labels
  label(px=.03, py=1.145, lab=paste0("reef fishes; n = ",length(data$reef[chosen1])), pos=4, cex=0.9, xpd=NA)
  label(px=.03, py=1.045, lab=paste0("other fishes; n = ",length(data$reef[chosen2])), pos=4, cex=0.9, xpd=NA)
  label(px=-0.25, py=0.1, lab=expression(paste("ln(Rate @ 20"*degree,"C) (g C d"^{-1}, ")"),sep=""), cex=1.3, srt=90, adj=c(0.5,0.5), xpd=NA)

  #(b)
  plot(1/8.62E-5*(1/293.15 - 1/data$temp)[chosen2], bw[chosen2], pch=21, col="grey70", bg=make.transparent("grey70", .5),  xlab="", ylab="", las=1, ylim=c(-12,1), cex.lab=1.3, xaxt='n')
  axis(side=3, at=1/8.62e-5*(1/293.15-1/(273.15+seq(0,40,by=5))), labels=seq(0,40,by=5))
  for(j in unique(data$fam)) {
    xpoints  <-  range(data$temp[data$fam==j])
    xpoints  <-  seq(xpoints[1], xpoints[2], length.out=50)
    nresEr   <-  acen + JAGSModel[paste0('r[',j,',2]'),'mean']
    nresboTc <-  bnot + JAGSModel[paste0('r[',j,',3]'),'mean']
    expr     <-  nresboTc + nresEr/8.62e-5*(1/293.15-1/xpoints)
    points(1/8.62e-5*(1/293.15-1/xpoints), expr, type="l", lty=2, col="grey30", lwd=1)
  }
  points(c(-4, 13), bnot + acen*c(-4, 13), type="l", lty=1, col="black", lwd=2.5)
  points(1/8.62E-5*(1/293.15 - 1/data$temp)[chosen1], bw[chosen1], pch=21, col="dodgerblue4", bg=make.transparent("dodgerblue4", .5))
  label(px=c(.85,.95), py=c(.95,.95), text=FALSE, type="l", lty=1, lwd=2.5, col="black")
  label(px=.85, py=.95, lab=substitute('mean trend: '~italic(y) == E*italic(x) - B, list(E=round(acen,2), B=-1*round(bnot,2))), pos=2, cex=0.9)
  label(px=.03, py=.97, expression(paste(bold("(b)"),sep="")), cex=1.4)
  label(px=.5, py=1.25, expression(paste("Temperature ("*degree,"C)",sep="")), xpd=NA, cex=1.3, adj=c(0.5,0.5))
  label(px=-0.25, py=0.1, lab=expression(paste("ln(Rate @ 1g) (g C d"^{-1},")"),sep=""), cex=1.3, srt=90, adj=c(0.5,0.5), xpd=NA)

  par(mai=c(1.7,0.82,0.14,0.60))
  fig2(data=data, nlModel=nlModel, combs=TRUE)
}

figS5  <-  function(data, nlmerModel) {
  #creates fam-specific parameters for NLMER results
  
  nres   <-  coef(nlmerModel)$family
  nboTc  <-  nres$ln.boTc[match(data$family,rownames(nres))]
  nA     <-  nres$a[match(data$family,rownames(nres))]
  nEi    <-  nres$Ei[match(data$family,rownames(nres))]
  nEr    <-  nres$Er[match(data$family,rownames(nres))]
  nTopt  <-  nres$Topt[match(data$family,rownames(nres))]
  xnboTc <-  fixef(nlmerModel)["ln.boTc"]
  xnEi   <-  fixef(nlmerModel)["Ei"]
  xnEr   <-  fixef(nlmerModel)["Er"]
  xnA    <-  fixef(nlmerModel)["a"]
  xnTopt <-  fixef(nlmerModel)["Topt"]
  
  #corrects metabolic rates for temperature first (both lmer and nlmer) and then size scaling
  btNlmer  <-  data$rate - nEr/8.62e-5*(1/293 - 1/data$temp) + log(1 + exp(nEi/8.62e-5*(1/nTopt - 1/data$temp))*nEr/(nEi-nEr))
  bwNlmer  <-  data$rate - nA*data$weight
  
  par(mfrow=c(1, 2), family="Times", mar=c(5.1,4.1,4.1,3), omi=c(0.5,1,0.5,1))
  chosen1  <-  data$reef=="yes"
  chosen2  <-  data$reef=="no"
  
  #(a)
  plot(btNlmer[chosen2] ~ data$weight[chosen2], pch=21, col="grey70", bg=make.transparent("grey70", .5),  xlab="ln(Mass) (g)", ylab=expression(paste("ln(Rate @ 20"*degree,"C) (g C d"^{-1}, ")"),sep=""), las=1, xlim=c(-5,10), ylim=c(-10,8), cex.lab=1.3, xpd=NA)
  
  for(j in 1:nrow(nres)) {
    xpoints  <-  range(data$weight[as.character(data$family)==rownames(nres)[j]])
    expr     <-  nres$ln.boTc[j] + nres$a[j]*xpoints
    points(xpoints, expr, type="l", lty=2, col="grey30", lwd=1)
  }
  points(c(-8, 12), xnboTc + xnA*c(-8, 12), type="l", lty=1, col="black", lwd=2.5)
  points(btNlmer[chosen1] ~ data$weight[chosen1], pch=21, col="dodgerblue4", bg=make.transparent("dodgerblue4", .5))
  #fig position label
  label(px=.03, py=.97, expression(paste(bold("(a)"),sep="")), cex=1.4)
  #trends labels
  label(px=c(.85,.95), py=c(.95,.95), text=FALSE, type="l", lty=1, lwd=2.5, col="black")
  label(px=c(.85,.95), py=c(.87,.87), text=FALSE, type="l", lty=2, lwd=1, col="grey30")
  label(px=.85, py=.95, lab=substitute("mean trend: "~italic(y)== B + A*italic(x), list(B=substr(xnboTc,1,5), A=round(xnA,2))), pos=2, cex=0.9)
  label(px=.85, py=.87, lab="family-level variation", pos=2, cex=0.9)
  
  #color labels
  label(px=.03, py=1.15, text=FALSE, type="o", pch=21, col="dodgerblue4", bg=make.transparent("dodgerblue4", .5), xpd=NA, cex=1.2)
  label(px=.03, py=1.05, text=FALSE, type="o", pch=21, col="grey70", bg=make.transparent("grey70", .5), xpd=NA, cex=1.2)
  #text for color labels
  label(px=.03, py=1.145, lab=paste0("reef fishes; n = ",length(data$reef[chosen1])), pos=4, cex=0.9, xpd=NA)
  label(px=.03, py=1.045, lab=paste0("other fishes; n = ",length(data$reef[chosen2])), pos=4, cex=0.9, xpd=NA)
  
  #(b)
  plot(1/8.62E-5*(1/293.15 - 1/data$temp)[chosen2], bwNlmer[chosen2], pch=21, col="grey70", bg=make.transparent("grey70", .5),  xlab=expression(paste("Inverse Temperature, 1/", italic(kT[c]), " - 1/", italic(kT), " (1/eV)", sep="")), ylab=expression(paste("ln(Rate @ 1g) (g C d"^{-1},")"),sep=""), las=1, ylim=c(-12,2), cex.lab=1.3, xpd=NA)
  axis(side=3, at=1/8.62e-5*(1/293.15-1/(273.15+seq(0,40,by=5))), labels=seq(0,40,by=5))
  
  for(j in 1:nrow(nres)) {
    xpoints  <-  range(data$temp[data$family==rownames(nres)[j]])
    xpoints  <-  seq(xpoints[1], xpoints[2], length.out=50)
    nresEr   <-  nres$Er[j]
    expr     <-  nres$ln.boTc[j] + nresEr/8.62e-5*(1/293.15-1/xpoints) - log(1 + exp(nres$Ei[j]/8.62e-5*(1/nres$Topt[j] - 1/xpoints))*nresEr/(nres$Ei[j]-nresEr))
    points(1/8.62e-5*(1/293.15-1/xpoints), expr, type="l", lty=2, col="grey30", lwd=1)
  }
  
  mean.nx           <-  1/8.62e-5*(1/293.15-1/(263:323))
  mean.expr.nlmer  <-  xnboTc + xnEr/8.62e-5*(1/293.15-1/(263:323)) - log(1 + exp(xnEi/8.62e-5*(1/xnTopt - 1/(263:323)))*xnEr/(xnEi-xnEr))
  points(mean.nx, mean.expr.nlmer, type="l", lty=1, col="black", lwd=2.5)
  points(1/8.62E-5*(1/293.15 - 1/data$temp)[chosen1], bwNlmer[chosen1], pch=21, col="dodgerblue4", bg=make.transparent("dodgerblue4", .5))
  label(px=c(.94,.99), py=c(.79,.79), text=FALSE, type="l", lty=1, lwd=2.5, col="black")
  
  label(px=0.1, py=.92, 
    lab=
       substitute(italic(y) == xnboTc + xnEr*italic(x), 
         list(xnboTc=substr(xnboTc,1,5), 
           xnEr=round(xnEr,2) 
         )
       ), 
    pos=4, cex=0.9
  )
 
  label(px=0.1, py=.77, 
    lab=
      substitute(-ln*bgroup('(', 1 + finexpr*'exp'*bgroup('[', expr, ']'), ')'), 
        list(xnEr=round(xnEr,2), 
            xnEi=round(xnEi,2), 
            finexpr=round(xnEr/(xnEi-xnEr),2), 
            expr=substitute(xnEi*bgroup('(',frac(1,italic(k)*xnTopt) - frac(1,italic(k)*italic(T[c])) + italic(x),')'), 
                list(xnEi=round(xnEi,2), 
                  xnTopt=round(xnTopt)
                )
            )
          )
      ),
    pos=4, cex=0.9
  )
  label(px=.03, py=.97, expression(paste(bold("(b)"),sep="")), cex=1.4)
  label(px=.5, py=1.25, expression(paste("Temperature ("*degree,"C)",sep="")), xpd=NA, cex=1.3, adj=c(0.5,0.5))
}
