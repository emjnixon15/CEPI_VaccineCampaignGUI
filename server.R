# Load necessary packages ----------
library(shiny)
library(extraDistr)
library(igraph)

source("outbreakSizeFunction.R")
source("helperFunctions.R")

# Define server logic for random distribution application
shinyServer(function(input, output, session) {
  
  # Create the object with no values
  values <- reactiveValues(spillover=NULL, outbreaks=NULL, campaign=NULL, population=NULL, warnIncomplete=F, campaignStats=NULL)
  
  # load default parameters
  default_params <- reactive({
    tbl <- read.csv("data/parameters.csv")
    if(input$disease=="")
      idx <- grep(pattern="dxv", x=tbl$disease)
    else{
      idx <- grep(pattern=input$disease, x=tbl$disease)
    }
    return(tbl[idx,])
  })
  
  # set default parameters of choosen diseas 
  observeEvent(input$disease,{
    validate(
      need(input$disease!="", "Please select a disease.")
    )
    tbl <- default_params()
    values$campaign <- NULL
    values$outbreaks <- NULL
    values$spillover <- NULL
    values$population <- NULL
    values$warnIncomplete <- F
    values$campaignStats <- NULL

    if(!is.na(tbl[,"spo_mean"]) & !is.na(tbl[,"spo_dispersion"]) & !is.na(tbl[,"spo_z"])){
      var = tbl[,"spo_mean"]*(1+tbl[,"spo_mean"]/tbl[,"spo_dispersion"])
      updateNumericInput(session, inputId="spo_mean", value=round(tbl[,"spo_mean"],2))
      updateNumericInput(session, inputId="spo_sd", value=round(sqrt(var),2))
      updateNumericInput(session, inputId="spo_z", value=round(tbl[,"spo_z"],2))
    }
    if(!is.na(tbl[,"r0_mean"])){
      updateNumericInput(session, inputId="R0_mean", value=round(tbl[,"r0_mean"],2))
      updateCheckboxInput(session, inputId="R0_nbDist", label="Negative Binomial distribution?", value=!is.na(tbl[,"r0_dispersion"]))
      if(!is.na(tbl[,"r0_dispersion"])){
        var = tbl[,"r0_mean"]*(1+tbl[,"r0_mean"]/tbl[,"r0_dispersion"])
        updateNumericInput(session, inputId="R0_sd", value=round(sqrt(var),2))
      }else{
        updateNumericInput(session, inputId="R0_sd", value=NA)
      }
    }
    if(!is.na(tbl[,"inc_shape"]) & !is.na(tbl[,"inc_rate"])){
      updateNumericInput(session, inputId="inc_mu", value=round(tbl[,"inc_shape"]/tbl[,"inc_rate"],2))
      updateNumericInput(session, inputId="inc_sd", value=round(sqrt(tbl[,"inc_shape"]/tbl[,"inc_rate"]^2),2))
    }
    if(!is.na(tbl[,"rec_shape"]) & !is.na(tbl[,"rec_rate"])){
      updateNumericInput(session, inputId="rec_mu", value=round(tbl[,"rec_shape"]/tbl[,"rec_rate"],2)) # mu = shape/rate
      updateNumericInput(session, inputId="rec_sd", value=round(sqrt(tbl[,"rec_shape"]/tbl[,"rec_rate"]^2),2)) # var = shape/rate^2
    }
    if(!is.na(tbl[,"sea_shape1"]) & !is.na(tbl[,"sea_shape2"])){
      param <- calcBetaMuVar(tbl[,"sea_shape1"], tbl[,"sea_shape2"])
      updateNumericInput(session, inputId="sea_mu", value=round(param$mu*365/7,1))
      updateNumericInput(session, inputId="sea_sd", value=round(sqrt(param$var)*365/7,1))
    }
    if(!is.na(tbl[,"cluster_shape1"]) & !is.na(tbl[,"cluster_shape2"])){
      param <- calcBetaMuVar(tbl[,"cluster_shape1"], tbl[,"cluster_shape2"])
      updateNumericInput(session, inputId="clu_mu", value=round(param$mu,7))
      updateNumericInput(session, inputId="clu_sd", value=round(sqrt(param$var),7))
    }
    if(!is.na(tbl[,"cluster_amount"]))
      updateNumericInput(session, inputId="numCluster", value=tbl[,"cluster_amount"])
    if(!is.na(tbl[,"cluster_mean"]) & !is.na(tbl[,"cluster_dispersion"])){
      var = tbl[,"cluster_mean"]*(1+tbl[,"cluster_mean"]/tbl[,"cluster_dispersion"])
      updateNumericInput(session, inputId="cluster_mean", value=round(tbl[,"cluster_mean"],2))
      updateNumericInput(session, inputId="cluster_sd", value=round(sqrt(var),2))
    }


    updateNumericInput(session, inputId="tho_cases", value=tbl[,"tho_cases"])
    updateNumericInput(session, inputId="tho_days", value=tbl[,"tho_days"])
    updateNumericInput(session, inputId="vaccov_hcw", value=tbl[,"vaccov_hcw"])
    updateNumericInput(session, inputId="vaccov_pop", value=tbl[,"vaccov_pop"])
    updateNumericInput(session, inputId="vaceff_dose1", value=tbl[,"vaceff_dose1"])
    updateCheckboxInput(session, inputId="vaceff_2dose", value=!is.na(tbl[,"vaceff_dose2"]))
    updateNumericInput(session, inputId="vaceff_dose2", value=tbl[,"vaceff_dose2"])
    updateNumericInput(session, inputId="vacdel_dose1", value=tbl[,"vacdel_dose1"])
    updateNumericInput(session, inputId="vacdel_dose2", value=tbl[,"vacdel_dose2"])
    updateNumericInput(session, inputId="vacprot", value=tbl[,"vacprot"])
    
    updateTabsetPanel(session, inputId="paramTab", selected="epiTab")
    
  })
  
  # plot functions for parameters
  output$defaultNumCluster <- renderText({ 
    validate(
      need(input$disease!="", "Please select a disease.")
    )
    tbl <- default_params()
    paste("Number of catchment areas:", tbl[,"cluster_amount"])
  })
  
  output$defaultTypeCluster <- renderText({ 
    validate(
      need(input$disease!="", "Please select a disease.")
    )
    tbl <- default_params()
    paste("Type of catchment areas:", tbl[,"cluster_type"])
  })
  
  output$SpillOverSize = renderUI({ plotOutput(outputId="plotSpillOver", width="100%", height="200px") })
  output$plotSpillOver <- renderPlot(width="auto", height="auto", {
    tbl <- default_params()
    validate(
      need(input$disease!="" & !is.na(tbl[,"spo_mean"]), "Please select a disease.")
    )
    spo_size <- input$spo_mean^2/(input$spo_sd^2 - input$spo_mean)
    xmax <- max(5,
                wrap_qzinb(0.99, mu=input$spo_mean, size=spo_size, pi=input$spo_z),
                wrap_qzinb(0.99, mu=tbl[,"spo_mean"], size=tbl[,"spo_dispersion"], pi=tbl[,"spo_z"]),
                na.rm = T)
    xseq <- seq(0, round(xmax,0), 1)
    yseqDefault <- wrap_dzinb(x=xseq, mu=tbl[,"spo_mean"], size=tbl[,"spo_dispersion"], pi=tbl[,"spo_z"])
    yseqInput <- wrap_dzinb(x=xseq, mu=input$spo_mean, size=spo_size, pi=input$spo_z)
    ymax <- max(yseqDefault, yseqInput, na.rm = T)
    par(mar=c(c(4, 6, 0.5, 0.5) + 0.1))
    plot(x=xseq, y=yseqDefault, ylim=c(0, ymax), type="l", xlab="Size", ylab="Density\n", las=1, lwd=2)
    abline(v=tbl[,"spo_mean"], col=1, lty=2) # TODO does the plotted mean value needs adjustment?
    if(!is.na(input$spo_mean)){
      lines(x=xseq, y=yseqInput, col=2)
      abline(v=input$spo_mean, col=2, lty=2) # TODO does the plotted mean value needs adjustment?
    }
  })
  
  output$R0 = renderUI({ plotOutput(outputId="plotR0", width="100%", height="200px") })
  output$plotR0 <- renderPlot({
    tbl <- default_params()
    validate(
      need(input$disease!="" & !is.na(tbl[,"r0_mean"]), "Please select a disease.")
    )
    
    defaultSize <- tbl[,"r0_dispersion"]
    if(is.na(defaultSize))
      defaultSize <- Inf
    inputSize <- input$R0_mean^2/(input$R0_sd^2 - input$R0_mean) 
    if(is.na(inputSize) | !input$R0_nbDist)
      inputSize <- Inf
    
    xmax <- max(5,
                qnbinom(0.99, mu=input$R0_mean, size=inputSize),
                qnbinom(0.99, mu=tbl[,"r0_mean"], size=defaultSize),
                na.rm = T)
    xseq <- seq(0, round(xmax,0), 1)
    yseqDefault <- dnbinom(x=xseq, mu=tbl[,"r0_mean"], size=defaultSize)
    yseqInput <- dnbinom(x=xseq, mu=input$R0_mean, size=inputSize)
    ymax <- max(yseqDefault, yseqInput, na.rm = T)
    
    par(mar=c(c(4, 6, 0.5, 0.5) + 0.1))
    plot(x=xseq, y=yseqDefault, type="l", xlab="Offsprings", ylab="Density\n", las=1, lwd=2)
    abline(v=tbl[,"r0_mean"], col=1, lty=2)
    if(!is.na(input$R0_mean)){
      lines(x=xseq, y=yseqInput, col=2)
      abline(v=input$R0_mean, col=2, lty=2)
    }
  })
  
  output$IncubationTime = renderUI({ plotOutput(outputId="plotIncubation", width="100%", height="200px") })
  output$plotIncubation <- renderPlot({
    tbl <- default_params()
    validate(
      need(input$disease!="" & !is.na(tbl[,"inc_shape"]) & !is.na(tbl[,"inc_rate"]), "Please select a disease.")
    )
    shape <- input$inc_mu^2/input$inc_sd^2
    rate <- input$inc_mu/input$inc_sd^2
    xmax <- max(10,
                qgamma(0.99, shape=shape, rate=rate),
                qgamma(0.99, shape=tbl[,"inc_shape"], rate=tbl[,"inc_rate"]),
                na.rm = T)
    xseq <- seq(0, round(xmax,0), length.out = 50)
    yseqDefault <- dgamma(x=xseq, shape=tbl[,"inc_shape"], rate=tbl[,"inc_rate"])
    yseqInput <- dgamma(x=xseq, shape=shape, rate=rate)
    ymax <- max(yseqDefault, yseqInput, na.rm = T)

    par(mar=c(c(4, 6, 0.5, 0.5) + 0.1))
    plot(x=xseq, y=yseqDefault,
         type="l", xlab="Days", ylab="Density\n", las=1, lwd=2)
    abline(v=tbl[,"inc_shape"]/tbl[,"inc_rate"], col=1, lty=2)
    if(!is.na(input$inc_mu) & !is.na(input$inc_sd)){
      lines(x=xseq, y=yseqInput, col=2)
      abline(v=input$inc_mu, col=2, lty=2)
    }
  })

  output$RecoveryTime = renderUI({ plotOutput(outputId="plotRecovery", width="100%", height="200px") })
  output$plotRecovery <- renderPlot({
    tbl <- default_params()
    validate(
      need(input$disease!="" & !is.na(tbl[,"rec_shape"]) & !is.na(tbl[,"rec_rate"]), "Please select a disease.")
    )
    shape <- input$rec_mu^2/input$rec_sd^2
    rate <- input$rec_mu/input$rec_sd^2
    xmax <- max(15,
                qgamma(0.99, shape=shape, rate=rate),
                qgamma(0.99, shape=tbl[,"rec_shape"], rate=tbl[,"rec_rate"]),
                na.rm = T)
    xseq <- seq(0, round(xmax,0), length.out = 50)
    yseqDefault <- dgamma(x=xseq, shape=tbl[,"rec_shape"], rate=tbl[,"rec_rate"])
    yseqInput <- dgamma(x=xseq, shape=shape, rate=rate)
    ymax <- max(yseqDefault, yseqInput, na.rm = T)
    par(mar=c(c(4, 6, 0.5, 0.5) + 0.1))
    plot(x=xseq, y=yseqDefault, ylim=c(0, ymax),
         type="l", xlab="Days", ylab="Density\n", las=1, col=1, lwd=2)
    abline(v=tbl[,"rec_shape"]/tbl[,"rec_rate"], col=1, lty=2)
    if(!is.na(input$rec_mu) & !is.na(input$rec_sd)){
      lines(x=xseq, y=yseqInput, col=2)
      abline(v=input$rec_mu, col=2, lty=2)
    }
  })
  
  output$Timing = renderUI({ plotOutput(outputId="plotTiming", width="100%", height="200px") })
  output$plotTiming <- renderPlot({
    tbl <- default_params()
    validate(
      need(input$disease!="" & !is.na(tbl[,"sea_shape1"]) & !is.na(tbl[,"sea_shape2"]), "Please select a disease.")
    )
    param <- estBetaParams(input$sea_mu/365*7, (input$sea_sd/365*7)^2)
    xseq <- seq(1/365, 1-1/365, 1/365)
    yseqDefault <- dbeta(x=xseq, shape1=tbl[,"sea_shape1"], shape2=tbl[,"sea_shape2"])
    yseqInput <- dbeta(x=xseq, shape1=param$shape1, shape2=param$shape2)
    ymax <- max(yseqDefault, yseqInput, na.rm = T)

    par(mar=c(c(4, 6, 0.5, 0.5) + 0.1))
    sea_start <- default_params()[,"sea_start"]
    plot(x=seq(1, 365-1, 1), y=yseqDefault, ylim=c(0, ymax),
         type="l", xlab=sprintf("Weeks since %s 1st", sea_start), ylab="Density\n", las=1, axes=F, lwd=2)
    abline(v=365*tbl[,"sea_shape1"]/(tbl[,"sea_shape1"]+tbl[,"sea_shape2"]), col=1, lty=2)
    axis(1, at=seq(0,365,7*4), labels=seq(0,52,4), las=1); axis(2, las=1); box()
    if(!is.na(input$sea_mu) & !is.na(input$sea_sd)){
      lines(x=seq(1, 365-1, 1), y=yseqInput, col=2)
      abline(v=input$sea_mu*7, col=2, lty=2)
    }
  })
  
  output$SpillOverProb = renderUI({ plotOutput(outputId="plotSpilloverProb", width="100%", height="200px") })
  output$plotSpilloverProb <- renderPlot({
    tbl <- default_params()
    validate(
      need(input$disease!="" & !is.na(tbl[,"cluster_shape1"]) & !is.na(tbl[,"cluster_shape2"]), "Please select a disease.")
    )
    param <- estBetaParams(input$clu_mu, input$clu_sd^2)
    xseq <- seq(0.01, 0.99, 0.01)
    yseqDefault <- dbeta(x=xseq, shape1=tbl[,"cluster_shape1"], shape2=tbl[,"cluster_shape2"])
    yseqInput <- dbeta(x=xseq, shape1=param$shape1, shape2=param$shape2)
    ymax <- max(yseqDefault, yseqInput, na.rm = T)
    
    par(mar=c(c(4, 6, 0.5, 0.5) + 0.1))
    plot(x=xseq, y=yseqDefault, ylim=c(0, ymax), 
         type="l", xlab="Spillover probability", ylab="Density\n", las=1, lwd=2)
    abline(v=tbl[,"cluster_shape1"]/(tbl[,"cluster_shape1"]+tbl[,"cluster_shape2"]), col=1, lty=2)
    if(!is.na(input$clu_mu) & !is.na(input$clu_sd)){
      lines(x=xseq, y=yseqInput, col=2)
      abline(v=input$clu_mu, col=2, lty=2)
    }
  })
  
  output$Population = renderUI({ plotOutput(outputId="plotPopulation", width="100%", height="200px") })
  output$plotPopulation <- renderPlot({
    tbl <- default_params()
    validate(
      need(input$disease!="", "Please select a disease.")
    )
    cluster_size <- input$cluster_mean^2/(input$cluster_sd^2 - input$cluster_mean) 
    xmax <- max(1000,
                qnbinom(0.99, mu=input$cluster_mean, size=cluster_size),
                qnbinom(0.99, mu=tbl[,"cluster_mean"], size=tbl[,"cluster_dispersion"]),
                na.rm = T)
    xseq <- round(seq(0, xmax, length.out = 1000),0)
    yseqDefault <- dnbinom(x=xseq, mu=tbl[,"cluster_mean"], size=tbl[,"cluster_dispersion"])
    yseqInput <- dnbinom(x=xseq, mu=input$cluster_mean, size=cluster_size)
    ymax <- max(yseqDefault, yseqInput, na.rm = T)
    
    par(mar=c(c(4, 6, 0.5, 0.5) + 0.1))
    plot(x=xseq, y=yseqDefault, ylim=c(0, ymax), type="l", xlab="Size", ylab="Density\n", las=1, lwd=2)
    abline(v=tbl[,"cluster_mean"], col=1, lty=2)
    if(!is.na(input$cluster_mean)){
      lines(x=xseq, y=yseqInput, col=2)
      abline(v=input$cluster_mean, col=2, lty=2)
    }
  })
  
  output$globalPopEst <- renderText({ 
    validate(
      need(input$disease!="", "Please select a disease.")
    )
    format(input$cluster_mean * input$numCluster, big.mark   = ",", scientific=F)
  })
  
  # output$outbreakSizeProb <- renderText({ 
  #   validate(
  #     need(input$disease!="", "Please select a disease.")
  #   )
  #   ln.outbreakSizeProb(outbreakSize=input$maxOutbreakSize, founders=input$spo_mean, R=0.99)
  #   format(input$cluster_mean * input$numCluster, big.mark   = ",", scientific=F)
  # })  
  
  ##########################
  # runs outbreak simulation
  #outbreaks <- eventReactive(input$goButton, {
  observeEvent(input$goButton, {
    
    # check input variable
    if(input$disease=="") showNotification("Please select a disease and press the GO button.", type ="error")
    validate(
      need(input$disease != "", "Please select a disease")
    )
    
    updateTabsetPanel(session, inputId="paramTab", selected="resTab")
    
    tbl <- default_params()

    withProgress(message = 'Running outbreak simulation', value = 0, expr = {
      set.seed(1) # set random value generator to same starting value
      spo_size <- input$spo_mean^2/(input$spo_sd^2 - input$spo_mean)
      values$spillover <-  wrap_rzinb(n=input$replicats, mu=input$spo_mean, size=spo_size, pi=input$spo_z)
      
      paramCluster <- estBetaParams(input$clu_mu, input$clu_sd^2)
      spillover <- do.call(cbind, lapply(values$spillover, function(sp){
        rmultinom(1, sp, prob=rbeta(input$numCluster, shape1=paramCluster$shape1, shape2=paramCluster$shape2))          
      }))
      
      values$warnIncomplete <- F
      values$outbreaks$baseline <- lapply(1:input$replicats, function(i){
        # Increment the progress bar, and update the detail text.
        incProgress(amount=1/input$replicats)
        ntw <- lapply(spillover[,i], function(sp){
          runOutbreakSimulation(input, R=input$R0_mean, spillover=sp)
        })
        values$warnIncomplete <- any(values$warnIncomplete, 
                                     unlist(lapply(ntw, function(nn) nn$warnIncomplete)))
        return(ntw)
      })
    })
  })
  
  
  baseline <- reactive({
    validate(
      need(!is.null(values$outbreaks$baseline), "Please select a disease and press the GO button.")
    )
    values$outbreaks$baseline
  })
  
  # observeEvent(output$warnIncomplete, {
  #   validate(
  #     need(!is.null(output$warnIncomplete), "")
  #   )
  #   showModal(modalDialog(
  #     title = "Warning",
  #     paste0("Warning: One or more outbreaks were larger than ", input$maxOutbreakSize, ". You might need to increase the 'Maximum allowed Outbreak Size' parameter in the simulation tab.")
  #   ))
  # })
  
  output$warnIncomplete <- renderText({ 
    validate(
      need(!is.null(values$outbreaks$baseline), "Please select a disease and press the GO button.")
    )
    if(values$warnIncomplete){
      paste0("Warning: One or more outbreaks were larger than ", input$maxOutbreakSize,
         ". The estimated numbers might be underestimated. You might need to increase the 'Maximum allowed Outbreak Size' parameter in the simulation tab.")
    }else NULL
  })
  
  # summary campain outcome
  campaign <- reactive({
    ntwLst <- baseline()
    values$campaign <- evalCampainSimulation(ntwLst)
  })

  output$CampaignStats = renderUI({ tableOutput("tableCasesAverted") })
  output$tableCasesAverted <- renderTable(rownames=T, digits=1, align="c", striped=T, expr={
    withProgress(message = 'Calculate required regimes', value = 0, expr = {
      cluster_size <- input$cluster_mean^2/(input$cluster_sd^2 - input$cluster_mean) 
      lst <- do.call(rbind, lapply(campaign(), function(tt){
        population <- rnbinom(input$numCluster, mu=input$cluster_mean, size=cluster_size)
        data.frame(Cases = sum(tt[,"outbreakSize"]),
                   Spillover = sum(tt[,"numSpillover"]),
                   Averted = sum(ifelse(!is.na(tt[,"vacCampainStarted"]), tt[,"numAverted"], 0)),
                   Regimes_Contact = input$numContact*sum(tt[,"numRing"]),
                   Regimes_Population = sum(population[!is.na(tt[,"vacCampainStarted"])]),
                   Campaigns = sum(!is.na(tt[,"vacCampainStarted"]))
                   )
      }))
      stats <- statsSummary(lst)
      colnames(stats) <- gsub("_"," ",names(lst))
      incProgress(amount=1)
    })
    values$campaignStats <- stats
    return(stats)
  })
  
  output$DistSpillover = renderUI({ plotOutput(outputId="plotDistSpillover", width="100%", height="200px") })
  output$plotDistSpillover <- renderPlot({
    validate(
      need(!is.null(values$spillover), "Please select a disease and press the GO button.")
    )
    par(mar=c(4,6,0.5,0.5)+0.1)
    hist(values$spillover, breaks=50, probability=T, col='darkgray', border="white", las=1,
         xlab='Spillover cases', ylab="Density\n", main=NULL)
    xseq <- seq(0, max(values$spillover), 1)
    spo_size <- input$spo_mean^2/(input$spo_sd^2 - input$spo_mean)
    lines(x=xseq, 
          y=wrap_dzinb(x=xseq, mu=input$spo_mean, size=spo_size, pi=input$spo_z),
          col=2, lwd=2)
    box()
    legend("topright", legend = c("Simulated","Distribution"), lty=c(1,1), lwd=c(6,2), col=c("gray","red"), bty="n")
  })
  
  output$TimmingCases = renderUI({ plotOutput(outputId="plotTimmingCases", width="100%", height="200px") })
  output$plotTimmingCases <- renderPlot({
    withProgress(message = 'Running timming evaluation', value = 0, expr = {
      ntwLst <- baseline()
      sea_start <- default_params()[,"sea_start"]
  
      timeSpillover <- unlist(lapply(ntwLst, function(ntw){
        incProgress(amount=0.5/length(ntwLst))
        unlist(lapply(ntw, function(nn){
          spIdx <- nn$edges[,"Parent"]==0
          nn$times[spIdx,"time.infec"]
        }))
      }))
      timeSpillover <- floor(timeSpillover/7)
    
      timeH2H <- unlist(lapply(ntwLst, function(ntw){
        incProgress(amount=0.5/length(ntwLst))
        unlist(lapply(ntw, function(nn){
          h2hIdx <- nn$edges[,"Parent"]!=0
          nn$times[h2hIdx,"time.infec"]
        }))
      }))
      timeH2H <- floor(timeH2H/7)
 
      timeRange <- as.character(1:max(52, timeH2H, timeSpillover, na.rm=T))
      tabSpillover <- table(timeSpillover)[timeRange]/length(ntwLst)
      tabSpillover[is.na(tabSpillover)] <- 0
      tabH2H <- table(timeH2H)[timeRange]/length(ntwLst)
      tabH2H[is.na(tabH2H)] <- 0
      tab <- rbind(tabSpillover, tabH2H)
      tab[is.na(tab)] <- 0
      
      par(mar=c(4,6,0.5,0.5)+0.1)
      barplot(tab, beside = F, width=1, space=0, names.arg=rep(NA,dim(tab)[2]),
              las=1, col=c("darkgray","lightgray"), border="white",
              ylab="Cases\n", xlab=sprintf("Weeks since %s 1st", sea_start))
      axis(1, las=1); box()
      legend("topright", legend=c("Spillover","Human-to-human"), fill=c("darkgray","lightgray"), bty="n")
        
      x <- seq(0,1,length.out=52)
      paramSea <- estBetaParams(input$sea_mu/365*7, (input$sea_sd/365*7)^2)
      y <- dbeta(x, paramSea$shape1, paramSea$shape2)
      xTrans <- floor(x*(365/7))
      spillOverDist <- cbind(xTrans, y)
      scale <- mean(tabSpillover, na.rm=T)
      lines(x=spillOverDist[,1]-0.5, y=scale*spillOverDist[,2], col=2, lwd=2)
    })
  })
  
  output$DistCampaignStart = renderUI({ plotOutput(outputId="plotDistCampaignStart", width="100%", height="200px") })
  output$plotDistCampaignStart <- renderPlot({
    withProgress(message = 'Plot campain start distribution', value = 0, expr = {
      sea_start <- default_params()[,"sea_start"]
      # Campain start distribution
      startMat <- do.call(rbind, lapply(campaign(), function(tab){
        tab[,"vacCampainStarted"]
      }))
      if(all(is.na(startMat))){
        distCampainStart <- rep(0,52)
        names(distCampainStart) <- as.character(1:52)
        scale <- 1
      } else {
        distCampainStart <- table(floor(startMat/7))[as.character(1:floor(max(startMat, na.rm=T)/7))]
        # campaignProb <- rowSums(startMat>0, na.rm=T)/dim(startMat)[1]
        # campaignProb <- sum(start>0, na.rm=T)/length(start)
        xMax <- max(52, as.integer(names(distCampainStart)), na.rm=T)
        distCampainStart <- distCampainStart[as.character(1:xMax)]
        names(distCampainStart) <- as.character(1:xMax)
        distCampainStart[is.na(distCampainStart)] <- 0
        scale <- mean(distCampainStart/dim(startMat)[1], na.rm=T)
      }
      # seasonal timming
      x <- seq(0,1,length.out=100)
      paramSea <- estBetaParams(input$sea_mu/365*7, (input$sea_sd/365*7)^2)
      y <- dbeta(x, paramSea$shape1, paramSea$shape2)
      xTrans <- x*(365/7)
      spillOverDist <- cbind(xTrans, y)
      ymax <- max(distCampainStart/dim(startMat)[1], scale*spillOverDist[,2], na.rm=T)
      
      par(mar=c(4,6,0.5,0.5)+0.1)
      barplot(distCampainStart/dim(startMat)[1], width=1, space=0, 
              names.arg=rep(NA,length(distCampainStart)),
              las=1, col="darkgray", border="white", ylim=c(0,ymax),
              ylab="Average campaigns\nper week", xlab=sprintf("Weeks since %s 1st", sea_start))
      axis(1, las=1); box()
      lines(x=spillOverDist[,1], y=scale*spillOverDist[,2], col=2, lwd=2)
      legend("topright", legend = c("Campaigns","Seasonality"), lty=c(1,1), lwd=c(6,2), col=c("darkgray","red"), bty="n")
      #legend("topright", legend=sprintf("Total probability = %.2f", campaignProb), bty="n")
      incProgress(1)
    })
  })
  
  # Downloadable csv of simulation data
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("Campaign_Simulator_", input$disease, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(values$campaignStats, file, row.names=T)
    }
  )
  
})