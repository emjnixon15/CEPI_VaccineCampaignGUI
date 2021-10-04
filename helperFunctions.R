
runOutbreakSimulation <- function(input, R, spillover){
  
  ##SKIP IF THERE IS NO SPILLOVER
  if(spillover == 0){
    return(list())
  }else{
    # simulate network
    total_cases <- input$maxOutbreakSize
    paramSea <- estBetaParams(input$sea_mu/365*7, (input$sea_sd/365*7)^2)
    inc_shape <- input$inc_mu^2/input$inc_sd^2
    inc_rate <- input$inc_mu/input$inc_sd^2
    rec_shape <- input$rec_mu^2/input$rec_sd^2
    rec_rate <- input$rec_mu/input$rec_sd^2
    ntw <- simulateOutbreakNetwork(total.cases = total_cases, 
                                   index.cases = spillover, 
                                   R = R,
                                   R.dispersion = ifelse(input$R0_nbDist, input$R0_mean^2/(input$R0_sd^2 - input$R0_mean), NA),
                                   shape1.infec = paramSea$shape1,
                                   shape2.infec = paramSea$shape2,
                                   shape.incub = inc_shape,
                                   scale.incub = 1/inc_rate,
                                   shape.death = rec_shape,
                                   scale.death = 1/rec_rate)
    ntw.vac <- simulateVaccination(ntw,
                                   thresholdCases = input$tho_cases, 
                                   thresholdTimeFrame = input$tho_days,
                                   vacFrac_sp = input$vaccov_pop, 
                                   vacFrac_h2h = input$vaccov_pop,
                                   vacDelay = input$vacdel_dose1, 
                                   vac2Delay = input$vacdel_dose2, 
                                   vacDuration = 0,
                                   vacEfficacy = input$vaceff_dose1, 
                                   vac2Efficacy = input$vaceff_dose2, 
                                   protDelay = input$vacprot)
    if(sum(ntw.vac$incomplete)>0){
      ntw.vac$warnIncomplete <- T
    }else{
      ntw.vac$warnIncomplete <- F 
    }
    return(ntw.vac)
  }
}

evalCampainSimulation <- function(ntwLstAA){
  withProgress(message = 'Running campaign evaluation', value = 0, expr = {
  lapply(ntwLstAA, function(ntwLst){
    incProgress(amount=1/length(ntwLstAA))
    do.call(rbind, lapply(seq_along(ntwLst), function(i){
      # SKIP IF THERE IS NO SPILLOVER
      if(is.null(ntwLst[[i]]$VacCampainStarted)){
        return(c(outbreakSize=0, 
                 numSpillover=0,
                 # numH2h=NA,
                 numRing=0,
                 numAverted=NA,
                 numAvertedHCW=NA,
                 numVacOnce=NA,
                 numVacTwice=NA,
                 numProtected=NA,
                 numDeath=NA,
                 numDeathAverted=NA,
                 numDeathAvertedHCW=NA,
                 vacCampainStarted=NA))
      }else{
        if(is.null(ntwLst[[i]])) stop(sprintf("Network should not be NULL, i=%i and j=%s",i,j))
        # VACCINATION RESPONSE
        if(!is.na(ntwLst[[i]]$VacCampainStarted))
          numRing <- sum(ntwLst[[i]]$times[,"time.incub"]>=ntwLst[[i]]$VacCampainStarted)
        else
          numRing <- 0
        return(c(outbreakSize=dim(ntwLst[[i]]$edges)[1], 
                 numSpillover=sum(ntwLst[[i]]$edges[,1]==0),
                 # numH2h=sum(ntwLst[[i]]$edges[,1]!=0),
                 numRing=numRing,
                 numAverted=sum(ntwLst[[i]]$toPrun),
                 numAvertedHCW=sum(ntwLst[[i]]$toPrun & ntwLst[[i]]$edges[,1]!=0),
                 numVacOnce=sum(!is.na(ntwLst[[i]]$times[,"time.vac1"])),
                 numVacTwice=sum(!is.na(ntwLst[[i]]$times[,"time.vac2"])),
                 numProtected=sum(!is.na(ntwLst[[i]]$times[,"time.prot"])),
                 numDeath=sum(ntwLst[[i]]$death),
                 numDeathAverted=sum(ntwLst[[i]]$death & ntwLst[[i]]$toPrun),
                 numDeathAvertedHCW=sum(ntwLst[[i]]$death & ntwLst[[i]]$toPrun & ntwLst[[i]]$edges[,1]!=0),
                 vacCampainStarted=ntwLst[[i]]$VacCampainStarted))
      }
    }))
  })
  })
}

statsSummary <- function(lst){
  do.call(cbind, lapply(lst, function(ll){
    qt <- quantile(ll, probs=c(0.25, 0.75), na.rm=T)
    rbind(Min=min(ll, na.rm=T),
          '25%'=qt[1],
          Median=median(ll, na.rm=T),
          Mean=mean(ll, na.rm=T),
          '75%'=qt[2],
          Max=max(ll, na.rm=T))
  }))
}

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(shape1=alpha, shape2=beta))
}

calcBetaMuVar <- function(alpha, beta) {
  mu <- alpha / (alpha+beta)
  var <- alpha * beta / ((alpha+beta)^2 * (alpha+beta+1))
  return(list(mu=mu, var=var))
}

wrap_dzinb <- function(x, size, mu, pi){
  prob <- size/(size+mu)
  data <- dzinb(x=x, size=size, prob=prob, pi=pi)
  return(data)
}

wrap_qzinb <- function(p, size, mu, pi){
  prob <- size/(size+mu)
  data <- qzinb(p=p, size=size, prob=prob, pi=pi)
  return(data)
}

wrap_rzinb <- function(n, size, mu, pi){
  prob <- size/(size+mu)
  data <- rzinb(n=n, size=size, prob=prob, pi=pi)
  return(data)
}

ln.outbreakSizeProb = function(outbreakSize, founders, R){
  log(founders) + 
    (outbreakSize - founders - 1) * log(outbreakSize) + 
    (outbreakSize - founders) * log(R) - 
    outbreakSize * R - 
    lfactorial(outbreakSize - founders)
}