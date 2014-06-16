make.transparent <- function(col, opacity=0.5) {
  if (length(opacity) > 1 && any(is.na(opacity))) {
    n <- max(length(col), length(opacity))
    opacity <- rep(opacity, length.out=n)
    col <- rep(col, length.out=n)
    ok <- !is.na(opacity)
    ret <- rep(NA, length(col))
    ret[ok] <- Recall(col[ok], opacity[ok])
    ret
  } else {
    tmp <- col2rgb(col)/255
    rgb(tmp[1,], tmp[2,], tmp[3,], alpha=opacity)
  }
}

label <- function(px, py, lab, adj=c(0, 1), text=TRUE, log=FALSE, ...) {
  usr  <-  par("usr")
  x.p  <-  usr[1] + px*(usr[2] - usr[1])
  y.p  <-  usr[3] + py*(usr[4] - usr[3])
  if(log=="x"){x.p<-10^(x.p)}
  if(log=="y"){y.p<-10^(y.p)}
  if(log=="xy"){x.p<-10^(x.p);y.p<-10^(y.p)}
  if(text){
    text(x.p, y.p, lab, adj=adj, ...)
  } else {
    points(x.p, y.p, ...)
  }
}

fig.path  <-  function(output.path="output.figs", name) {
  file.path(get(output.path), name)
}

createDir  <-  function(name1, name2) {
  if(missing(name1))
    name1  <-  getwd()
  if(!(name2 %in% list.files(name1)))
    dir.create(file.path(name1,name2))
}

to.dev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  dev(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

to.pdf <- function(expr, filename, ...) {
  to.dev(expr, pdf, filename, ...)
}

## lme4 model Functions
FullModel <- function(ln.boTc, logit.Er, Ea, Ei, Topt, a, mass, avkT, temp) 
  log(exp(ln.boTc)*exp(Ea*avkT)*mass^a*
        exp((Ei/(1 + exp(-logit.Er)))/8.62e-5*(1/293.15 - 1/temp))/(1 + exp(Ei/8.62e-5*(1/Topt - 1/temp))*((Ei/(1 + exp(-logit.Er)))/(Ei-(Ei/(1 + exp(-logit.Er)))))))
FullModel <- deriv(body(FullModel), namevec = c("ln.boTc","logit.Er","Ea","Ei","Topt","a"), func = FullModel)

FullMinusEa <- function(ln.boTc, logit.Er, Ei, Topt, a, mass, temp)
  log(exp(ln.boTc)*mass^a*
        exp((Ei/(1 + exp(-logit.Er)))/8.62e-5*(1/293.15 - 1/temp))/(1 + exp(Ei/8.62e-5*(1/Topt - 1/temp))*((Ei/(1 + exp(-logit.Er)))/(Ei-(Ei/(1 + exp(-logit.Er)))))))
FullMinusEa <- deriv(body(FullMinusEa), namevec = c("ln.boTc", "logit.Er", "Ei", "Topt", "a"), func = FullMinusEa)

FullMinusIT <- function(ln.boTc, logit.Er, Ea, a, Ei, mass, avkT, temp)
  log(exp(ln.boTc)*exp(Ea*avkT)*mass^a*
        exp((Ei/(1 + exp(-logit.Er)))/8.62e-5*(1/293.15 - 1/temp)))
FullMinusIT <- deriv(body(FullMinusIT), namevec = c("ln.boTc", "logit.Er", "Ea", "a"), func = FullMinusIT)

FullMinusEaIT <- function(ln.boTc, logit.Er, a, Ei, mass, temp)
  log(exp(ln.boTc)*mass^a*
        exp((Ei/(1 + exp(-logit.Er)))/8.62e-5*(1/293.15 - 1/temp)))
FullMinusEaIT <- deriv(body(FullMinusEaIT), namevec = c("ln.boTc", "logit.Er", "a"), func = FullMinusEaIT)

St <- function(ln.boTc, Er, Ei, Topt, a, mass, temp) 
  log(exp(ln.boTc)*mass^a*
      exp(Er/8.62e-5*(1/293 - 1/temp))/(1 + exp(Ei/8.62e-5*(1/Topt - 1/temp))*(Er/(Ei-Er))))
St <- deriv(body(St), namevec = c("ln.boTc", "Er", "Ei", "Topt", "a"), func = St)

hardWireBayesian  <-  function(iterationMatrix, community_dataset, temperature_dataset, lnNPPvector) {
  if(!(all(rownames(temperature_dataset)==names(lnNPPvector)))) stop('names between inputs do not match')

  communityTable  <-  matrix(0, nrow(iterationMatrix), (length(unique(community_dataset$sites))*6*3)+49) #6 means the mean plus 5 trophic groups; 3 means ctot, flux and epsilon. and then 49 means the boltzmann values
  
    for(i in 1:nrow(iterationMatrix))  {
        params      <-  iterationMatrix[i, c('A','Er','Ei','Topt','boTc')]
        ctots       <-  daply(community_dataset, .(sites, trophnum), function(x){log(sum(x$abun*(x$biomass_g^params['A']))) - log(unique(x$site.area))})
        ctots       <-  ctots[rownames(temperature_dataset),]
        ctots       <-  cbind(ctots, log(rowSums(exp(ctots), na.rm=TRUE)))
        boltz       <-  (exp(params['Er']/8.62E-5*(1/293.15 - 1/temperature_dataset))/(1 + (params['Er']/(params['Ei']-params['Er']))*exp(params['Ei']/8.62E-5*(1/params['Topt'] - 1/temperature_dataset))))
        boltz       <-  log(rowMeans(boltz)) 
        ln.flux     <-  ctots   + boltz + params['boTc']      
        ln.epsilon  <-  ln.flux - lnNPPvector
    
        communityTable[i,]         <-  c(sapply(c(ctots, ln.flux, ln.epsilon), function(x)as.numeric(x)), boltz)
    }
  
    colnames(communityTable)  <-  c(as.character(sapply(paste0(rep(c("ln.sizeCorrectedBiomass","ln.EnergyFlux","ln.Epsilon"), each=6), c(1:5,"")), function(x)paste0(x,"[", 1:49, "]"))), paste("timeAveragedKinetics[", 1:49, "]", sep=""))
    
    communityTable
}

