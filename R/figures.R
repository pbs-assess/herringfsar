############################################################
# figures.R
# 
# Plotting code for making Herring FSAR figures.
# 
# Last Modified: Aug 22, 2024
# 
############################################################


# Plot stock indicators as four panels
plot_indicators <- function(dat) {
  # Catch and TAC (panel A)
  p_catch <- ggplot(data = dat, mapping = aes(x = Year, y = Catch)) +
    geom_line() +
    geom_hline(yintercept = dat$TAC, linetype = "dashed") +
    scale_y_continuous(labels = comma) +
    labs( y = "Catch (1,000 t)") +
    expand_limits(y = 0) +
    annotate(
      geom = "text", x = Inf, y = Inf, label = "(A)", vjust = 2, hjust = 2
    ) +
    theme(axis.text.x = element_blank())
  # SSB, USR, and LRP (panel B)
  p_biomass <- ggplot(data = dat, mapping = aes(x = Year, y = SSB_med)) +
    geom_ribbon(mapping = aes(ymin = SSB_min, ymax = SSB_max), fill = "grey") +
    geom_line() +
    geom_hline(yintercept = dat$LRP, linetype = "dashed", colour = "red") +
    geom_hline(yintercept = dat$USR, linetype = "dotted", lwd=1.0, colour = "darkgreen") +
    geom_hline(yintercept = dat$B0, linetype = "dashed", colour = "black") +
    scale_y_continuous(labels = comma) +
    labs(y = "Biomass (1,000 t)") +
    expand_limits(y = 0) +
    annotate(
      geom = "text", x = Inf, y = Inf, label = "(B)", vjust = 2, hjust = 2
    ) +
    theme(axis.text.x = element_blank())
  # F and M (panel C)
  p_mortality <- ggplot(data = dat, mapping = aes(x = Year, y = M_med)) +
    geom_ribbon(mapping = aes(ymin = M_min, ymax = M_max), fill = "salmon", alpha.f = 0.5) +
    geom_line(mapping = aes(y = M_med), colour = "red") +
    labs(y = "Natural Mortality (/yr)") +
    expand_limits(y = 0) +
    annotate(
      geom = "text", x = Inf, y = Inf, label = "(C)", vjust = 2, hjust = 2
    )  +
    theme(axis.text.x = element_blank())
  # F and M (panel D)
  p_hr <- ggplot(data = dat, mapping = aes(x = Year, y = U_med)) +
    geom_ribbon(mapping = aes(ymin = U_min, ymax = U_max), fill = "grey") +
    geom_line() +
    geom_hline(yintercept = dat$Uref, linetype = "dashed", colour = "red") +
    labs(y = "Harvest Rate") +
    expand_limits(y = 0) +
    annotate(
      geom = "text", x = Inf, y = Inf, label = "(D)", vjust = 2, hjust = 2
    )
  # Recruitment (panel E)
  p_recruitment <- ggplot(data = dat, mapping = aes(x = Year, y = R_med)) +
    geom_ribbon(mapping = aes(ymin = R_min, ymax = R_max), fill = "grey") +
    geom_line() +
    scale_y_continuous(labels = comma) +
    labs(y = "Recruitment (millions)") + # 1,000
    expand_limits(y = 0) +
    annotate(
      geom = "text", x = Inf, y = Inf, label = "(E)", vjust = 2, hjust = 2
    )
  # Recruitment (panel F)
  p_surplusprod <- ggplot(data = dat, mapping = aes(x = SSB_med, y = SP_med)) +
    geom_point() +
    scale_y_continuous(labels = comma) +
    labs(y = "Surplus Production (kt)", x = "Spawning Biomass (kt)") + # 1,000
    expand_limits(y = 0) +
    annotate(
      geom = "text", x = Inf, y = Inf, label = "(F)", vjust = 2, hjust = 2
    )
  # Six-panel plot
  p <- p_catch + p_biomass + p_mortality + p_hr + p_recruitment + p_surplusprod
  p
}


# baseplot_indicators()
# Base-plot version of the multi-panel
# mandatory plot for the FSAR, showing
# most recent assessments of catch, biomass
# mortality, harvest rate, recruitment,
# and surplus production.
# Inputs: dat = table of data made from
#               an MS3 sim obj (usually weighted OM)
# Outputs: NULL, but plots to active graphics device.
baseplot_indicators <- function(dat) 
{
  
  par(mfrow = c(3,2), mar = c(2.5,1,2,1), oma = c(4,5,2,5))

  yrs <- dat$Year

  # Catch and TAC (A)
  plot(x = range(yrs), y = c(0,max(dat$Catch)),
        type = "n", xlab = "", ylab = "", axes = FALSE )
    mfg <- par("mfg")
    grid()
    box()
    axis(side = 1)
    axis(side = 2, las = 1)
    mtext(side = 2, text = "Catch (1,000 t)", line = 4)
    mtext(side = 1, text = "Year", line = 3)
    legend(x = "topright", bty = "n", legend = "(A)")
    lines(x = yrs, y = dat$Catch )

  # SSB (B)
  plot(x = range(yrs), y = c(0,max(dat$SSB_max)),
        type = "n", xlab = "", ylab = "", axes = FALSE )
    mfg <- par("mfg")
    grid()
    box()
    axis(side = 1)
    axis(side = 4, las = 1)
    mtext(side = 4, text = "Spawning Biomass\n(1,000 t)", line = 4.5)
    mtext(side = 1, text = "Year", line = 3)
    legend(x = "topright", bty = "n", legend = "(B)")
    polygon(x = c(yrs,rev(yrs)), y = c(dat$SSB_min,rev(dat$SSB_max)),
            col = "grey65", border = NA)
    lines(x = yrs, y = dat$SSB_med )
    abline(h = dat$B0[1], lty = 2, lwd = 2)
    abline(h = dat$USR[1], lty = 4, lwd = 2, col = "darkgreen")
    abline(h = dat$LRP[1], lty = 2, lwd = 2, col = "red")

  # Mortality (C)
  plot( x = range(yrs), y = c(0,max(dat$M_max)),
        type = "n", xlab = "", ylab = "", axes = FALSE )
    mfg <- par("mfg")
    grid()
    box()
    axis(side = 1)
    axis(side = 2, las = 1)
    mtext(side = 2, text = "Natural Mortality (/yr)", line = 4)
    mtext(side = 1, text = "Year", line = 3)
    legend(x = "topright", bty = "n", legend = "(C)")
    polygon(x = c(yrs,rev(yrs)), y = c(dat$M_min,rev(dat$M_max)),
            col = scales::alpha("red",0.5), border = NA)
    lines(x = yrs, y = dat$M_med, col = "red" )
    
  # Harvest Rate (D)
  plot( x = range(yrs), y = c(0,max(dat$U_max)),
        type = "n", xlab = "", ylab = "", axes = FALSE )
    mfg <- par("mfg")
    grid()
    box()
    axis(side = 1)
    axis(side = 4, las = 1)
    mtext(side = 4, text = "Harvest Rate", line = 3)
    mtext(side = 1, text = "Year", line = 3)
    legend(x = "topright", bty = "n", legend = "(D)")
    polygon(x = c(yrs,rev(yrs)), y = c(dat$U_min,rev(dat$U_max)),
            col = scales::alpha("black",0.5), border = NA)
    lines(x = yrs, y = dat$U_med, col = "black" )
    segments(x0 = 1983, x1 = 2022, y0 = 0.2, lty = 2, lwd = 2)


  # Recruitment (E)
  plot( x = range(yrs), y = c(0,max(dat$R_max)),
        type = "n", xlab = "", ylab = "", axes = FALSE )
    mfg <- par("mfg")
    grid()
    box()
    axis(side = 1)
    axis(side = 2, las = 1)
    mtext(side = 2, text = "Recruitment (millions)", line = 4)
    mtext(side = 1, text = "Year", line = 3)
    legend(x = "topright", bty = "n", legend = "(E)")
    polygon(x = c(yrs,rev(yrs)), y = c(dat$R_min,rev(dat$R_max)),
            col = scales::alpha("black",0.5), border = NA)
    lines(x = yrs, y = dat$R_med, col = "black" )
    abline(h = mean(dat$R_med), lty = 2, lwd = 2)

  # Surplus Production (D)
  yrIdx <- which(yrs >= 1987 & yrs < max(yrs))
  colfunc <- colorRampPalette(c("black", "grey75"))
  colVec <- colfunc(length(yrIdx))
  ptTypeVec <- c(17,rep(16,length(yrIdx)-2),15)
  plot( x = c(0,max(dat$SSB_med[yrIdx])), y = range(dat$SP_med[yrIdx], na.rm = T),
        type = "n", xlab = "", ylab = "", axes = FALSE )
    mfg <- par("mfg")
    grid()
    box()
    axis(side = 1)
    axis(side = 4, las = 1)
    mtext(side = 4, text = "Surplus Production\n(1,000 t)", line = 4.5)
    mtext(side = 1, text = "Spawning Biomass (1,000 t)", line = 3)
    legend(x = "topright", bty = "n", legend = "(F)")
    legend( x = "topleft", 
            bty = "n", legend = c("1988",max(yrs)),
            pch = c(17,15), col = colVec[c(length(yrIdx),1)]  )
    lines(x = dat$SSB_med[yrIdx], y = dat$SP_med[yrIdx], col = "grey70", lwd = 1 )
    points( x = dat$SSB_med[yrIdx], y = dat$SP_med[yrIdx], 
            pch = ptTypeVec, col = rev(colVec) )
    abline(v = dat$USR[1], lty = 4, lwd = 2, col = "darkgreen" )
    abline(v = dat$LRP[1], lty = 2, lwd = 2, col = "red" )
    abline(v = dat$B0[1], lty = 2, lwd = 2, col = "black" )

} # END baseplot_indicators



# plotSimDataEnvelopes()
# Simulated and real catch and index data envelopes,
# meant to help diagnose exceptional circumstances.
# Inputs:
#     obj = MS3 simulation object
#     sIdx = species index for MS3 obj
#     pIdx = area index for MS3 obj
#     fIdx = fleet index for index data
#     Cdata = input real catch data series
#     Idata = input real idx data series
#     fYear = first year of the x-axis
# Returns: Nothing, but plots to active graphics device
plotSimDataEnvelopes <- function( obj = blob,
                                  sIdx = 1, pIdx =1, fIdx = 5,
                                  Cdata = NULL,
                                  Idata = NULL,
                                  fYear = 1951 )
{
  nT        <- obj$om$nT
  tMP       <- obj$om$tMP
  goodReps  <- obj$goodReps
  B_it      <- obj$om$SB_ispt[goodReps,sIdx,pIdx,]
  C_it      <- obj$om$C_ispt[goodReps,sIdx,pIdx,]
  I_it      <- obj$mp$data$I_ispft[goodReps,sIdx,pIdx,fIdx,] 

  yrs <- seq(from = fYear, by = 1, length.out = nT)

  B_qt <- apply(X = B_it, FUN = quantile, MARGIN = 2, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  C_qt <- apply(X = C_it, FUN = quantile, MARGIN = 2, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  I_qt <- apply(X = I_it, FUN = quantile, MARGIN = 2, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

  par(mfrow = c(2,1), oma = c(4,4,2,1), mar = c(.5,.5,.5,.5))

  plot(x = range(yrs), y = c(0,max(B_qt,I_qt)), las = 1,
        type = "n", axes = FALSE)
    axis(side = 2, las = 1)
    grid()
    box()
    mtext(side = 2, text = "Spawning Biomass (kt)", line = 3)
    
    polygon(x = c(yrs,rev(yrs)), y = c(B_qt[1,],rev(B_qt[3,])),
            border = NA, col = "grey75" )
    lines(x = yrs, y = B_qt[2,], lwd = 3 )

    polygon(x = c(yrs,rev(yrs)), y = c(I_qt[1,],rev(I_qt[3,])),
            border = NA, col = scales::alpha("steelblue",0.5) )

    if(!is.null(Idata))
    {
      nDat <- length(Idata)
      datCols <- rep("black",nDat)
      datCols[tMP:nDat] <- "salmon"
      
      points(x = yrs[1:nDat], y = Idata, col = datCols, pch = 16)

      legend( x = "topleft", bty = "n",
              legend = c( "Spawning Biomass",
                          paste0("Simulated spawn index range"),
                          paste0("Observed spawn index (", yrs[1]," - ", yrs[tMP-1],")"),
                          paste0("Observed spawn index (post ",yrs[tMP-1],")")),
              pch   = c(22,22,21,21),
              pt.bg = c("grey75","steelblue","black","salmon"),
              pt.lwd = 0,
              pt.cex = 1.5,
              col = c("black","black","black","salmon"),
              lwd   = c(3,NA,NA,NA) )
    }

    plot(x = range(yrs), y = c(0,max(C_qt)), las = 1,
        type = "n", axes = FALSE)
    axis(side = 2, las = 1)
    axis(side = 1 )
    grid()
    box()
    mtext(side = 2, text = "Catch (kt)", line = 3)
    polygon(x = c(yrs,rev(yrs)), y = c(C_qt[1,],rev(C_qt[3,])),
            border = NA, col = "grey75" )
    lines(x = yrs, y = C_qt[2,], lwd = 3 )

    if(!is.null(Cdata))
    {
      nDat <- length(Cdata)
      datCols <- rep("black",nDat)
      datCols[tMP:nDat] <- "salmon"
      
      points(x = yrs[1:nDat], y = Cdata, col = datCols, pch = 16)

      legend( x = "topright", bty = "n",
              legend = c( "Historical and simulated catch",
                          paste0("Observed catch (", yrs[1]," - ", yrs[tMP-1],")"),
                          paste0("Observed catch (post ",yrs[tMP-1],")")),
              pch   = c(22,21,21),
              pt.bg = c("grey75","black","salmon"),
              pt.lwd = 0,
              pt.cex = 1.5,
              col = c("black","black","salmon"),
              lwd   = c(3,NA,NA) )

    }

  mtext(side = 1, text = "Year", outer = TRUE, line = 2)
} # END plotSimDataEnvelopes()



# plotStockRecruitEnvelope()
# Beverton-Holt stock recruitment curve for an individual
# replicate with input ref pts picked out and 
# historical/projected recruitments shown. 
# inputs:
#   obj   = simulation object (blob)
#   iRep  = simulation replicate
# outputs:
#   no output object, but draws figure on plotting device
plotStockRecruitEnvelope <- function( obj = blob,
                                      SRvals = NULL,
                                      inputBmsy = 39.127,
                                      inputUSR = 66.97,
                                      inputB0 = 90.23 )
{
  goodReps  <- which(obj$goodReps)
  nReps     <- length(goodReps)
  nT        <- obj$om$nT
  tMP       <- obj$om$tMP
  
  B0_i  <- R0_i <- h_i <- rec.a_i <- rec.b_i <- array(0,nReps)

  for(i in 1:length(goodReps))
  {
    iIdx <- goodReps[i]
    
    B0_i[i]     <- obj$rp[[iIdx]]$B0_sp[,1,1]
    R0_i[i]     <- obj$rp[[iIdx]]$R0_sp[1,1]
    h_i[i]      <- obj$rp[[iIdx]]$h_sp[1,1]
    rec.a_i[i]  <- obj$rp[[iIdx]]$rec.a_sp[1,1]
    rec.b_i[i]  <- obj$rp[[iIdx]]$rec.b_sp[1,1]

  }

  if(!is.null(SRvals))
  {
    # SRvals is a dataframe with R0, B0, and h
    SRvals <- SRvals |>
                mutate( rec.a  = 4.*h*R0/(B0*(1.-h)),
                        rec.b  = (5.*h-1.)/(B0*(1.-h)))
  }

  SB_it  <- obj$om$SB_ispt[goodReps,1,1,]
  R_it   <- obj$om$R_ispt[goodReps,1,1,]

  randRep <- sample(1:nReps, 1)

  SB_qt <- apply(X = SB_it, FUN = quantile, MARGIN = 2, probs = c(0.025, 0.5, 0.975))
  R_qt  <- apply(X = R_it, FUN = quantile, MARGIN = 2, probs = c(0.025, 0.5, 0.975))

  # First, make the curve
  Bseq_k  <- seq(from =0, to = max(SB_qt), length.out = 100)
  Rseq_ik <- array(0, dim = c(nReps,100))
  
  for(i in 1:nReps)
    Rseq_ik[i,] <- rec.a_i[i] * Bseq_k/(1 + rec.b_i[i]*Bseq_k)

  Rseq_qk <- apply(X = Rseq_ik, FUN = quantile, MARGIN = 2, probs = c(0.025, 0.5, 0.975))

  Rspline <- splinefun(x = Bseq_k, y = Rseq_qk[2,])

  if(!is.null(inputBmsy))
    inputRmsy <- Rspline(inputBmsy)

  if(!is.null(inputUSR))
    inputRusr <- Rspline(inputUSR)

  if(!is.null(inputB0))
    inputR0 <- Rspline(inputB0)


  plot( x = c(0,max(1.4*SB_qt[2,])), y = c(0,1.2*max(Rseq_qk, 1.2*R_qt[2,])), type = "n", las =1,
        xaxs = "i", yaxs = "i",
        xlab = "Spawning Biomass (kt)", ylab = "Age-1 Recruitment (1e6)" )
    polygon(  x = c(Bseq_k,rev(Bseq_k)), y = c(Rseq_qk[1,],rev(Rseq_qk[3,])), 
              border = NA, col = "grey75" )
    

    if(!is.null(SRvals))
    {
      modCols <- RColorBrewer::brewer.pal(nrow(SRvals),"Dark2")
      for(k in 1:nrow(SRvals))
      {
        rec.a <- SRvals[k,"rec.a"]
        rec.b <- SRvals[k,"rec.b"]
        lines(  x = Bseq_k, y = rec.a * Bseq_k/(1 + rec.b * Bseq_k), 
                lwd = 2, col = modCols[k] )
      }

      legend( x = "topleft", bty = "n",
              legend = SRvals$Scenario, col = modCols, lwd = 2 )
    }

    lines( x = Bseq_k, y = Rseq_qk[2,], lwd = 3)
    segments( x0 = SB_qt[2,1:(tMP-2)], 
              y0 = R_qt[1,2:(tMP-1)],
              y1 = R_qt[3,2:(tMP-1)], lwd = 2, col = "grey30" )

    points( x = SB_qt[2,1:(tMP-2)], 
            y = R_qt[2,2:(tMP-1)], 
            pch = 16, col = "grey30" )
    points( x = SB_it[randRep,(tMP-1):(nT-1)], 
            y = R_it[randRep,tMP:nT], 
            pch = 16, col = "salmon" )

    
    if(!is.null(inputBmsy) & !is.null(inputB0) & !is.null(inputUSR))
    {
      segments( x0 = c(inputB0,inputBmsy,inputUSR),
                y0 = 0, 
                y1 = c(inputR0,inputRmsy,inputRusr), lty = 2, col = c("grey40","darkgreen","steelblue"), lwd = 2)
      segments( x0 = 0,
                x1 = c(inputB0,inputBmsy,inputUSR), 
                y0 = c(inputR0,inputRmsy,inputRusr), lty = 2, col = c("grey40","darkgreen","steelblue"), lwd = 2)
    }





    legend( x = "topright", bty = "n",
            legend = c( "Stock-Recruit Relationship",
                        "Unfished equilibrium",
                        "Provisional Upper Stock Reference",
                        "MSY equilibrium",
                        "Historical period recruitments",
                        "Projection period recruitments"),
            lty   = c(1,2,2,2,1,NA),
            lwd   = c(3,2,2,2,2,NA),
            pch   = c(22,NA,NA,NA,16,16),
            pt.bg = c("grey75",NA,NA,NA,NA,NA),
            pt.lwd = 0,
            pt.cex = 2, 
            col   = c("black","grey40","steelblue","darkgreen","grey30","salmon"))
} # END plotStockRecruitEnvelope()



# plotHCRules()
# Plots a generic hockey-stick
# harvest control rule, with given
# lower, upper control points, and
# low/high Fs. Stock status is calculated
# as a proportion of Bmsy, and fishing
# mortality rate as a proportion of
# Fmsy.
plotHCRules <- function(  LCP_i = c(.3,.3,.172),
                          UCP_i = c(.6,.6,.342),
                          lowF = .0,
                          highF = 1,
                          U_i = c(.14,.20,.187) ,
                          language = "English")
{
  x <- seq(0,1.3, length.out = 100 )
  y <- rep(lowF, length(x))


  par( mar = c(4,5,1,1), oma = c(2,2,1,1) )

  lineCols <- RColorBrewer::brewer.pal(3,"Dark2")

  plot( x = range(x), y = c(0,1.2*max(U_i)), type = "n",
        las = 1,
        xlab = "",
        ylab = "",
        cex.lab = 1.5,
        cex.axis = 1.5,
        axes = FALSE)
    grid()
  if(language=="English"){
    mtext(side = 2, text = "Target Harvest Rate", line = 4, cex = 1.5 )
    mtext(side = 1, text = expression(B/B[0]), line = 4, cex = 1.5 )
    box()
    axis(side=1, at=seq(0.0, max(x), by=0.2), labels = sprintf("%.1f", seq(0.0, max(x), by=0.2)))
    axis(side=2, at=seq(0.0, 1.2*max(U_i), by=0.02), labels = sprintf("%.2f", seq(0.0, 1.2*max(U_i), by=0.02)), las=1)
    # axis(side=4, at=seq(0.0, 1.2*max(U_i), by=0.02), labels = sprintf("%.2f", seq(0.0, 1.2*max(U_i), by=0.02)), las=1)
  }else{ #French
    mtext(side = 2, text = "Taux de récolte cible", line = 4, cex = 1.5 )
    mtext(side = 1, text = expression(B/B[0]), line = 4, cex = 1.5 )
    box()
    axis(side=1, at=seq(0.0, max(x), by=0.2), labels = chartr(".", ",", sprintf("%.1f", seq(0.0, max(x), by=0.2))))
    axis(side=2, at=seq(0.0, 1.2*max(U_i), by=0.02), labels = chartr(".", ",", sprintf("%.2f", seq(0.0, 1.2*max(U_i), by=0.02))), las=1)
  }
    for(i in 1:length(LCP_i))
    {
      segments( x0 = 0, x1 = LCP_i[i],
                y0 = lowF, y1 = lowF,
                col = "grey50", lwd = 1.5 )
      
      segments( x0 = LCP_i[i], x1 = UCP_i[i],
                y0 = lowF, y1 = highF*U_i[i],
                col = lineCols[i], lwd = 2 )
      segments( x0 = UCP_i[i], x1 = max(x),
                y0 = highF*U_i[i], y1 = highF*U_i[i],
                col = lineCols[i], lwd = 2 )
    }
    # abline( v = c(0.3, 0.73),
    #         col = c("red","darkgreen"),
    #         lty = 2, lwd = 2 )
    # abline( h = highF*maxSeq[1], lty = 2, col = "grey70" )

    legend( x = "topleft", bty = "n", col = lineCols,
            lwd = 2,
            legend = c( "30-60B0_maxTHR.14",
                        "30-60B0_maxTHR.20",
                        "40-80Bmsy_maxTHRUmsy"))

} # END plotHCRules()


# plotGridTulipBtCt()
# Simulation envelopes (tulips) of
# biomass, catch, and harvest rate, comparing
# all blobs in blobList.
# Inputs: blobList = list of MS3 blobs
#         labels = MP labels (if NULL uses blob)
#         dep = Plot biomass (FALSE) or depletion (TRUE)
#         yLimX = plot y limits
#         traces = number of replicates for traces
#         traceSeed = random seed for drawing reps
#         proj = plot only projection period?
#         maxProjT = maximum projection time step
# Outputs: plots to acrive graphics device
plotGridTulipBtCtUt <- function(  blobList = mpBlobList,
                                  labels = NULL,
                                  dep   = TRUE,
                                  yLimB   = c(0,2),
                                  yLimC   = c(0,15),
                                  yLimHR  = c(0,.3),
                                  traces = 3,
                                  traceSeed = 1234,
                                  proj  = TRUE,
                                  maxProjT = NULL)
{
  # Load blobs
  nSims <- length(blobList)
  simLabs <- c()
  nTs     <- c()
  for(simIdx in 1:nSims)
  {
    simLabs[simIdx] <- blobList[[simIdx]]$ctlList$ctl$mpName
    names(blobList)[simIdx] <- simLabs[simIdx]
    nTs[simIdx] <- blobList[[simIdx]]$om$nT
  }

  if(!is.null(labels))
    plotLabs <- labels
  else plotLabs <- simLabs


  # Now set up plotting area
  par(  mfcol = c(3,nSims),
        mar   = c(.1,.1,.1,.1),
        oma   = c(4,5,2,4) )

  fYear <- 1951
  nT <- blobList[[1]]$om$nT
  tMP <- blobList[[1]]$om$tMP

  yrs <- seq(from = fYear, length.out = max(nTs))
  if(is.null(maxProjT))
    maxProjT <- max(nTs) - tMP

  if(proj)
    xLim <- range(yrs[(tMP-3):(tMP + maxProjT)])
  else xLim <- range(yrs)

  goodReps <- blobList[[1]]$goodReps

  nReps <- sum(goodReps)
  nS    <- dim(blobList[[1]]$om$SB_ispt)[2]
  nP    <- dim(blobList[[1]]$om$SB_ispt)[3]
  nT    <- dim(blobList[[1]]$om$SB_ispt)[4]


  if(traces > 0)
  {
    set.seed(traceSeed)
    traceIdx <- sample(1:nReps, size = traces )
  }


  for(simIdx in 1:nSims)
  {
    blob <- blobList[[simIdx]]

    yrs <- seq(from = fYear, length.out = nTs[simIdx])
    # cat(simLabs[simIdx],"\n")

    # Get model states
    SB_ispt   <- blobList[[simIdx]]$om$SB_ispt[which(goodReps),,,,drop = FALSE]

    nF <- blob$om$nF
    C_ispt      <- blobList[[simIdx]]$om$C_ispt[which(goodReps),,,,drop = FALSE]
    uC_ispt     <- C_ispt
    
    if( nF > 5 )
    {
      uC_ispft    <- blobList[[simIdx]]$om$C_ispft
      P_ispft     <- blobList[[simIdx]]$om$P_ispft

      uC_ispft[,,,6,] <- P_ispft[,,,6,]
      
      uC_ispt <- apply(X = uC_ispft, FUN = sum, MARGIN = c(1,2,3,5), na.rm = T)
      uC_ispt <- uC_ispt[which(goodReps),,,,drop = FALSE]
    }
    
    # TAC_ispt  <- blobList[[simIdx]]$mp$hcr$TAC_ispt[1:nReps,,,,drop = FALSE]
    
    U_ispt    <- uC_ispt
    U_ispt    <- uC_ispt/(SB_ispt + uC_ispt)

    # Load ref pts
    B0_isp     <- array(NA,dim = c(nReps,nS,nP))
    Bmsy_isp   <- array(NA,dim = c(nReps,nS,nP))
    USR_isp    <- array(NA,dim = c(nReps,nS,nP))
    Bcrash_isp <- array(NA,dim = c(nReps,nS,nP))
    Ucrash_isp <- array(NA,dim = c(nReps,nS,nP))
    MSY_isp    <- array(NA,dim = c(nReps,nS,nP))
    # Umsy_isp  <- array(0,dim = c(nReps,nS,nP))
    # MSY_isp   <- array(0,dim = c(nReps,nS,nP))

    # Pull out ref pts, use ms3Rstats functions

    goodRepIdx <- which(goodReps)
    for( i in 1:nReps)
    { 
      # j <- goodRepIdx[i]
      # Pull out Dstar and Dscalar
      tdxUSR        <- blob$ctlList$opMod$histCtl$hypo$Dstar - 1950
      scalarUSR     <- blob$ctlList$opMod$histCtl$hypo$Dscalar
      # 
      # message("i = ", i, ", j = ", j)
      USR_isp[i,,]  <- apply(X = scalarUSR * SB_ispt[i,,,tdxUSR,drop = FALSE], FUN = mean, MARGIN = 2:3 )
      #
    }

    Busr_q <- quantile(USR_isp, probs = c(0.25, 0.5, 0.75), na.rm = T)
    
    # Pull B0, Bmsy, and Umsy from blob (should be saved from ctlList)
    B0_isp    <- blob$om$B0_isp[1:nReps,,,drop = FALSE]
    Bmsy_isp  <- blob$om$Bmsy_isp[1:nReps,,,drop = FALSE]
    Umsy_isp  <- blob$om$Umsy_isp[1:nReps,,,drop = FALSE]

    LRP_isp <- 0.3*B0_isp
    # USR_isp <- 0.8*Bmsy_isp
    
    # message("USR = ", Busr_q[2], ", B0 = ", B0_isp[1,1,1])

    sbylab <- "Spawning Biomass (kt)"

    if(dep)
    {
      sbylab <- expression(SB[t]/SB[0])
      for( t in 1:nTs[simIdx] )
      {
        SB_ispt[,,,t] <- SB_ispt[,,,t] / B0_isp
      }
      LRP_isp <- LRP_isp/B0_isp
      Bmsy_isp <-  Bmsy_isp/B0_isp
      USR_isp <- USR_isp/B0_isp
      B0_isp  <- B0_isp/B0_isp
      # USR_isp <- USR_isp/Bmsy_isp
      # TRP_isp <- TRP_isp/Bmsy_isp
    }
  
    B0    <- mean(B0_isp)
    LRP   <- mean(LRP_isp)
    Bmsy  <- mean(Bmsy_isp)
    USR   <- mean(USR_isp)
    Umsy  <- mean(Umsy_isp)

    SB_qspt   <- apply(X = SB_ispt, FUN = quantile, MARGIN = 2:4, na.rm = T,
                        probs = c(0.025, 0.5, .975 ) )

    C_qspt   <- apply(X = C_ispt, FUN = quantile, MARGIN = 2:4, na.rm = T,
                        probs = c(0.025, 0.5, .975 ) )

    U_qspt   <- apply(X = U_ispt, FUN = quantile, MARGIN = 2:4, na.rm = T,
                        probs = c(0.025, 0.5, .975 ) )


    plot( x = xLim, y = yLimB,
          type = "n", axes = FALSE )
      mfg <- par("mfg")      
        mtext(side = 3, text = plotLabs[simIdx], font = 2, cex = 1)

      if(mfg[2] == 1)
      {
        axis(side = 2, las = 1)
        mtext( side = 2, text = sbylab, line = 3 )
      }
      # if(mfg[2] == mfg[4]){
      #   axis(side = 4, las = 1)
      #   legend("topright", bty = "n", lty = c(2,2,5,2,5), col=c("grey30","darkgreen","purple","orange","red"), 
      #          c("B0","USR","Bmsy","LRP","Bcrash"), lwd = c(2,2,2,2,2))
      # }
      grid()
      box()

      polygon( x = c(yrs,rev(yrs)),
               y = c(SB_qspt[1,1,1,],rev(SB_qspt[3,1,1,])),
               col = "grey60", border = NA )
      lines( x = yrs, y = SB_qspt[2,1,1,], lwd = 3 )
      
      abline(v = yrs[tMP] - 0.5, lty = 2, lwd =.9 )
      abline(v = 2052, lty = 2, lwd =1 )
      
      
      abline(h = B0, lty = 2, col = "grey30" )
      abline(h = USR, lty = 2, col = "darkgreen", lwd=2)
      abline(h = Bmsy, lty = 5, col = "purple", lwd=2)
      abline(h = LRP, lty = 2, col = "red", lwd=2 )
      if(traces > 0)
        for( t in traceIdx )
          lines( x = yrs, y = SB_ispt[t,1,1,], lwd = 1, col = "black")
      if(mfg[2] == mfg[4])
      {
        axis(side = 4, las = 1)
        legend( "topright", 
                bty = "n", 
                lty = c(2,2,5,2,5), 
                col=c("grey30","darkgreen","purple","red"), 
                legend=c(bquote(SB[0]),"USR",bquote(B[MSY]),"LRP"), 
                lwd = c(1,2,2,2,2), cex=1.2)
      }

    plot( x = xLim, y = yLimC,
          type = "n", axes = FALSE )
      mfg <- par("mfg")      

      if( mfg[2] == 1 )
      {
        axis(side = 2, las = 1)
        mtext( side = 2, text = "Catch (kt)", line = 3 )
      }
      
      grid()
      box()

      polygon( x = c(yrs,rev(yrs)),
               y = c(C_qspt[1,1,1,],rev(C_qspt[3,1,1,])),
               col = "grey60", border = NA )
      lines( x = yrs, y = C_qspt[2,1,1,], lwd = 3 )      
      # abline( h = mean(MSY_isp[,1,1]), lty = 2, col = "darkgreen")
      abline(v = yrs[tMP] - 0.5, lty = 2, lwd = .9)

      if(traces > 0)
        for( t in traceIdx )
          lines( x = yrs, y = C_ispt[t,1,1,], lwd = 1, col = "black")

      if( mfg[2] == mfg[4] ){
        axis(side = 4, las = 1)
      }


    plot( x = xLim, y = yLimHR,
          type = "n", axes = FALSE )
      mfg <- par("mfg")      

      if( mfg[2] == 1 )
      {
        axis(side = 2, las = 1)
        mtext( side = 2, text = "Harvest Rate (/yr)", line = 3 )
      }

      axis(side = 1 )


      grid()
      box()

      polygon( x = c(yrs,rev(yrs)),
               y = c(U_qspt[1,1,1,],rev(U_qspt[3,1,1,])),
               col = "grey60", border = NA )
      lines( x = yrs, y = U_qspt[2,1,1,], lwd = 3 )      
      if(simLabs[simIdx] != "NoFish"){
        abline( h = blobList[[simIdx]]$ctlList$mp$hcr$Uref_p[1], lty = 2, col = "grey30")
        abline(h = Umsy, lty = 5, col = "purple", lwd=2)
      }
      abline(v = yrs[tMP] - 0.5, lty = 2, lwd =.9 )

      if(traces > 0)
        for( t in traceIdx )
          lines( x = yrs, y = U_ispt[t,1,1,], lwd = 1, col = "black")

      if( mfg[2] == mfg[4] )
      {
        axis(side = 4, las = 1)

        legend( x = "topright", 
                lty = 5,
                bty = "n", 
                col = "purple", 
                legend = bquote(U[MSY]), 
                lwd=2, 
                cex=1.2)
      }
      

  }

  mtext(side = 1, text=  "Year", outer = TRUE, line = 2)
} # END plotGridTulipBtCtUt()

# plotMulti()
# Plot a custom stack of plots that show model esitmates  of
# time-series.
# Inputs:   ts      = vector of plot names to stack
#           repList = report object output by SISCA
#           labcex  = axis label cex
#           heading = character vector of headings to use
# Outputs:  NULL
# Side-effects: plots to active graphics device.
plotMulti <- function(  ts = c("SBtIdx","Rt","Mt","Ftg"),
                        repList = reports,
                        labcex = 0.8, heading = "dimLabs" )
{

  nP <- repList$repOpt$nP
  par(mfcol = c(length(ts),nP), oma = c(3,4,3,4), mar = c(0.2,3,0.2,2),
        cex.lab = labcex )

  argList <- list(  repList = repList,  
                    noPar = TRUE, 
                    labcex = labcex )

  for( pIdx in 1:nP)  
  {
    argList$pIdx = pIdx
    for( tsName in ts )
    {
      if( tsName %in% c("SBt","SBtIdx") )
        plotArgs <- c(argList, list(heading = heading) )
      else plotArgs <- argList
      plotCall <- paste("plot",tsName,sep = "")
      do.call(  what = eval(plotCall), 
                args = plotArgs, 
                envir = .GlobalEnv )
    }
  }
} # END plotMulti()


# plotSBtIdx()
# Time-series of esimated spawning biomass with
# indices overlaid to show goodness of fit.
plotSBtIdx <- function( repList = reports,  
                        noPar = FALSE, 
                        pIdx = 1,
                        labcex = .8,
                        heading = NULL,
                        plotCt = FALSE )
{
  initYear <- repList$fYear
  report   <- repList$repOpt

  # Pull stuff from report
  SBt       <- report$SB_pt[pIdx,]
  Cgt       <- report$totC_pgt[pIdx,,]
  Igt       <- report$I_pgt[pIdx,,]
  qg        <- report$qhat_pg[pIdx,]
  rI_gt     <- repList$data$rI_pgt[pIdx,,]
  vulnBgt   <- report$vulnB_pgt[pIdx,,]
  vulnNgt   <- report$vulnN_pgt[pIdx,,]
  mat_a     <- report$mat_a

  combIdx <- FALSE

  if( any(repList$data$whichCombIdx_g > 0) )
  {
    combIdx <- TRUE
    # Replace q values with qComb_g
    qg <- report$qComb_pg[pIdx,]
  }
  
  vulnBgt[ vulnBgt == 0 ] <- NA
  vulnNgt[ vulnNgt == 0 ] <- NA



  # Pull prob positive idx.
  probPosIdx_gt <- report$probPosIdx_pgt[pIdx,,]

  # Get survey type and whether
  # idx is relative or absolute
  surveyType  <- report$survType_g
  indexType   <- report$indexType_g
  calcIndex   <- report$calcIndex_g
  surveyGears <- which(calcIndex == 1)

  # Get number of time steps/gear types
  nT    <- report$nT
  nG    <- report$nG
  nA    <- report$nA

  B0        <- signif(report$B0_p[pIdx],3)
  M0        <- signif(report$M0_p[pIdx],3)

  plotCI <- FALSE

  if(!is.null(repList$posts))
  {
    plotCI <- TRUE
    tInitModel_p <- repList$repOpt$tInitModel_p
    SB_it <- repList$posts$SB_ipt[,pIdx,]

    SB_qt <- apply(X = SB_it, FUN = quantile,
                    MARGIN = c(2), probs = c(0.025, 0.5, 0.975),
                    na.rm = TRUE )

    SBt <- apply(X = SB_it, FUN = mean, MARGIN = c(2))

    B0 <- signif(mean(repList$posts$B0_ip[,pIdx]),3)
    M0 <- signif(mean(repList$posts$M0_ip[,pIdx]),3)

    q_ig  <- repList$posts$q_ipg[,pIdx,]
    qg   <- apply( X = q_ig, FUN = mean, MARGIN = 2)

  } else {
    SB_qt <- array(0, dim = c(3,nT+1))
  }



  SBt[ SBt == 0 ] <- NA
  SB_qt[SB_qt == 0] <- NA

  # Sum catch
  Ct        <- apply(X = Cgt, FUN = sum, MARGIN = 2)

  stockLabs <- dimnames(report$SB_pt)[[1]][pIdx]

  if( !is.null(heading))
  {
    if( heading == "dimLabs")
      stockLabel <- stockLabs
    else
      stockLabel <- try(heading[stockLabs])

    if( class(stockLabel) == "try-error" ) 
      stockLabel <- stockLabs
  }

  # if( stockLabs == "Aggregate" )
  # {
  #   # We need to pull the aggregated stock
  #   # indices in at this point
  #   mItg    <- report$mI_tg
  #   mqg     <- report$qhat_g

  #   mItg[mItg < 0] <- 0

  #   for( g in surveyGears )
  #   {
  #     if( all(Itg[,g] < 0) & any(mItg[,g] > 0))
  #     {
  #       Itg[,g] <- mItg[,g]
  #       qg[g]   <- mqg[g]
  #     }
  #   }
  # }

  # browser()

  Igt[Igt<0] <- NA

  scaledIndices <- Igt

  for( g in surveyGears )
  {
    posIdx <- which(Igt[g,1:nT] > 0)

    if( surveyType[g] == 1 )
      vulnTS <- rI_gt[g,1:nT] * vulnNgt[g,1:nT]
    
    if( surveyType[g] == 0 )
      vulnTS <- rI_gt[g,1:nT] * vulnBgt[g,1:nT]

    if( surveyType[g] == 2 )
    {
      tmpB_t     <-  SBt

      vulnTS <- rI_gt[g,1:nT] * tmpB_t[1:nT]
    }
    
    # Compute scaled indices
    scaledIndices[g,posIdx] <- Igt[g,posIdx] / qg[g]  * SBt[posIdx] / vulnTS[posIdx] / probPosIdx_gt[g,posIdx] 

  }

  tInitModel <- report$tInitModel_p[pIdx] + 1

  scaledCombIdx <- repList$data$combI_pt[pIdx,]
  scaledCombIdx[scaledCombIdx < 0 ] <- NA

  if( combIdx )
  {
    
    probPosCombIdx <- report$probPosCombIdx_pt[pIdx,]
    qComb_t <- report$qComb_pt[pIdx,]

    vulnB_at <- report$vulnB_apgt[,pIdx,5,1:nT]

    for( a in 1:nA )
      vulnB_at[a,] <- vulnB_at[a,] * mat_a[a]
    tmpB_t     <-  apply(X = vulnB_at, FUN = sum, MARGIN = 2)
    tmpB_t     <-  SBt[1:nT]


    scaledCombIdx <- scaledCombIdx / qComb_t / probPosCombIdx * SBt[1:nT] / tmpB_t
  }


  # Create x axis label vector and vertical lines
  # for easy multipanel plotting
  years <- seq(from = initYear, length = nT, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  cols <- brewer.pal(nG,"Dark2")

  depRange  <- range(SBt/B0, na.rm = T)
  maxDep    <- ceiling(depRange[2])

  depAxis   <- round(seq(0,maxDep, length.out = 5),2)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  yMax <- max(SBt,scaledIndices, na.rm =T)

  # Now plot recruitments
  plot( x = range(years), y = c(0,yMax ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    mfg <- par("mfg")
    # Add SSB axis
    axis(side = 2, las =1 )
    # Add depletion axis
    axis( side = 4, las = 1, labels = depAxis, at = B0 * depAxis )

    if( mfg[1] == mfg[3] )
      axis(side = 1 )

    box()
    grid()
    # Plot data
    if( combIdx )
    {
      posIdx  <- which(scaledCombIdx > 0)
      zeroIdx <- which(scaledCombIdx == 0)
      zeroIdx <- zeroIdx[zeroIdx >= tInitModel ]
      points( x = years[posIdx], y = scaledCombIdx[posIdx],
              col = "grey50", pch = 16)

      # points( x = years[zeroIdx], y = scaledCombIdx[zeroIdx],
      #         col = "grey40", pch = 1)
      axis( side = 3, at = years[zeroIdx], labels = FALSE,
            col.ticks = "grey40", tck = .02, lwd.ticks = 3 )
    }

    # Plot Biomass
    if( plotCI )
    {
      polyCol <- scales::alpha("red", .5)
      polygon(  x = c(years[1:nT],rev(years[1:nT])),
                y = c(SB_qt[1,1:nT], rev(SB_qt[3,1:nT])),
                col = polyCol, border = NA )
    }
    abline( h = B0, lty = 3, lwd = .8, col = "red" )
    lines( x = years[1:nT], y = SBt[1:nT], lwd = 2, col = "red" )
    if(plotCt)
      rect( xleft = years[1:nT] - .3, xright = years[1:nT] + .3,
            ybottom = 0, ytop = Ct, border = NA, col = "grey60" )
    # points(x = years[(nT+1)], y = SBt[(nT+1)], col = "black", bg = "red", pch = 21)

    if( mfg[2] == 1 )
      mtext(side = 2, text = "Spawning\nBiomass (kt)", line = 3.5, cex = labcex)
    if( mfg[2] == mfg[4] )
      mtext(side = 4, text = "Depletion", line = 3, cex = labcex)
    if( mfg[1] == 1 & !is.null(heading) )
      mtext(side = 3, text = stockLabel, line = 1, font = 2, cex = labcex )

    for( g in surveyGears )
    {
      posIdx  <- which(scaledIndices[g,] > 0)
      zeroIdx <- which(scaledIndices[g,] == 0)
      points( x = years[posIdx], y = scaledIndices[g,posIdx],
            col = alpha(cols[g],.5), pch = 16 )
      points( x = years[zeroIdx], y = scaledIndices[g,zeroIdx],
            col = alpha(cols[g],.5), pch = 1 )
    }

    
    # text( x = years[nT - 7], y = c(.4,.3,.2,.1) * yMax,
    #       labels = c( paste( "B0 = ",  B0, sep = ""),
    #                   paste( "M0 = ", M0, sep = ""),
    #                   paste( "qs = ",  round(qg[4],2), sep = ""),
    #                   paste( "qd = ",  round(qg[5],2), sep = "")),

    #       cex = .9 )
} # END plotSBtIdx



# Plot fishing mortality
plotFtg <- function(  repList = reports,  
                      noPar = FALSE, 
                      pIdx = 1,
                      labcex = 1,
                      cBrewPal='Dark2' )
{
  report    <- repList$repOpt
  initYear  <- repList$fYear

  # Pull stuff from report
  Ugt     <- report$U_pgt[pIdx,,,drop = FALSE]
  gearIDs <- dimnames(report$U_pgt)[[2]]
  nG      <- report$nG
  nT      <- report$nT

  # Pull catch
  C_gt    <- report$totC_pgt[pIdx,,]
  C_t     <- apply( X = C_gt, FUN = sum, MARGIN = 2)

  SB_t    <- report$SB_pt[pIdx,1:nT]

  # Need to recalculate U based on herring U method
  U_gt    <- C_gt
  for(g in 1:nG)
    U_gt[g,] <- C_gt[g,]/(SB_t + C_t)

  U_t <- C_t/(SB_t + C_t)

  if(length(repList$posts))
  {
    tInitModel_p <- repList$repOpt$tInitModel_p

    U_igt <- repList$posts$U_ipgt[,pIdx,,]
    
    SB_it <- repList$posts$SB_ipt[,pIdx,1:nT]
    U_it <- SB_it

    for( i in 1:nrow(SB_it))
    {
      
      for( g in 1:nG )
      {
        U_igt[i,g,] <- C_gt[g,]/(SB_it[i,] + C_t)
      }
    }
    U_it <- apply(X = U_igt, FUN = sum, MARGIN = c(1,3))


    U_gt   <- apply( X = U_igt, FUN = median,
                      MARGIN = c(2,3), na.rm = T )

    U_t <- apply( X = U_it, FUN = median,
                      MARGIN = c(2), na.rm = T )

  } 
  


  commGears <- which(repList$ctlList$data$fleetType[gearIDs] != 0 )


  minTime_g <- rep(0,nG)
  maxTime_g <- rep(nT,nG)
  for( g in commGears )
  {
    minTime_g[g] <- min(which(Ugt[1,g,] > 0),na.rm = T)
    maxTime_g[g] <- max(which(Ugt[1,g,] > 0),na.rm = T)
  }
  minTime_g[!is.finite(minTime_g)] <- 1
  maxTime_g[!is.finite(maxTime_g)] <- nT


  cols    <- brewer.pal( length(commGears), cBrewPal )

  years <- seq(from = initYear, length.out = nT, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot F series
  plot( x = range(years), y = c(0,1 ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE, yaxs = "i" )
    mfg <- par("mfg")
    axis(side = 2, las =1 )
    if( mfg[1] == mfg[3] )
      axis(side = 1 )

    # Add catch axis
    maxC <- max(C_t,na.rm = T)
    maxCaxis <- round(5 * maxC / 4)
    CaxisLabs <- round(seq(from = 0, to = maxCaxis, length.out = 6),1)
    CaxisTicks <- seq(from = 0, to = 1, length.out = 6)
    axis( side = 4, at = CaxisTicks, labels = rev(CaxisLabs), las = 1 )
    box()
    # Plot recruitment
    grid()
    # Rescale catch so that it goes from 0 to .8?
    C_t <- C_t / maxC * 0.8
    rect( xleft = years - .3,
          xright = years + .3,
          ytop = 1, ybottom = 1 - C_t,
          col = "grey50", border = NA )

    for(gIdx in 1:length(commGears))
    {
      g <- commGears[gIdx]
      gYrs <-  max(minTime_g[g]-1,1):min(maxTime_g[g]+1,nT)
      lines( x = years[gYrs], y = U_gt[g,gYrs], lwd = 2, col = cols[gIdx] )
    }
    lines(x = years, y = U_t, lwd = 2, col = "grey30", lty = 2)
    
    if(mfg[2] == 1)
    {
      legend( x = "topright", col = c(cols,"grey30"),
              legend = c(gearIDs[commGears],"Total"), 
              lwd = 2, lty = c(rep(1,length(commGears)),2), 
              bty = "n")
      mtext(side = 2, text = "Harvest\nRate", line = 3.5, cex = labcex )
    }
    if( mfg[2] == mfg[4] )
      mtext( text = "Catch (kt)", side = 4, line = 3, cex = labcex)
}


# Plot natural mortality
plotMt <- function( repList = reports, 
                    noPar = FALSE, 
                    pIdx = 1,
                    labcex = .8)
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear

  # Pull stuff from report
  Mat       <- report$M_apt[,pIdx,,drop = FALSE]
  Mbar      <- report$M
  M         <- report$M_p[pIdx]
  M0        <- report$M0_p[pIdx]
  nT        <- report$nT
  Mjuve     <- report$Mjuve_p[pIdx]
  juveMage  <- report$juveMage

  # Pull info on predation mortality
  Ugt     <- report$U_pgt[pIdx,,,drop = TRUE]
  gearIDs <- dimnames(report$U_pgt)[[2]]
  nG      <- report$nG
  nT      <- report$nT

  predG   <- which(gearIDs %in% c("hakeLt50", "hakeGt50", "HS", "SSL", "HB"))

  Mat[Mat == 0] <- NA

  plotCI <- FALSE

  
  if(!is.null(repList$posts))
  {
    plotCI <- TRUE
    tInitModel_p  <- repList$repOpt$tInitModel_p
    M_iapt        <- repList$posts$M_iapt[,,pIdx,,drop = FALSE]
    Mjuve_i       <- repList$posts$M_iapt[,juveMage,pIdx,tInitModel_p[pIdx]+1]

    M_qapt <- apply(  X = M_iapt, FUN = quantile,
                      MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                      na.rm = TRUE )

    Mat <- apply(X = M_iapt, FUN = mean, MARGIN = c(2,3,4))

    Mjuve <- mean(Mjuve_i)


  } else {
    M_qapt <- 0
  }

  # if(addPredM & length(predG)>0 )
  # {

  #   predM   <- apply(Fgt[predG,], FUN=sum, MAR=c(2))
  #   basalM  <- Mat[juveMage+1,1,1:nT]
  #   Mt    <- predM + basalM

  #   YLIM <- c(0,max(Mt, na.rm =T) )      
  # }  else
  # YLIM <- c(0,max(Mat,M_qapt, na.rm =T) )
 

  years     <- seq(from = initYear, length.out = nT, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)



  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot F series
  plot( x = range(years), y = c(0,max(Mat,M_qapt, na.rm =T) ),
        type = "n", xlab = "", ylab = "", las = 1, axes = FALSE )
    mfg <- par("mfg")
    axis(side = 2, las =1 )
    if( mfg[1] == mfg[3] )
      axis(side = 1 )
    box()
    grid()
    # Plot recruitment
    polyCol <- scales::alpha("salmon", .5)
    if( plotCI )
    {
      polygon(  x = c(years[1:nT],rev(years[1:nT])),
                y = c(M_qapt[1,juveMage+1,1,1:nT], rev(M_qapt[3,juveMage+1,1,1:nT])),
                col = polyCol, border = NA )
    }
    
    # abline( h = M, lty = 2, lwd = 2, col = "salmon")
    abline( h = Mbar, lty = 2, lwd = 2, col = "grey50")
    abline( h = Mjuve, lty = 4, lwd = 2, col = "salmon" )
    lines( x = years[1:nT], y = Mat[juveMage+1,1,1:nT], col = "salmon",
            lwd = 3  )
    # points( x = years[nT+1], y = Mat[juveMage+1,1,nT+1], col = "black",
    #         pch = 21, bg = "salmon" )

    # if(addPredM & length(predG)>0 )
    # {
    #   lines( x = years[1:nT], y = predM, col = "#1b9e77", lwd = 3  )
    #   lines( x = years[1:nT], y = Mt, col = "black", lwd = 3  )

    #   legend( x = "top", 
    #           lwd = c(3,3),
    #           lty = c(1,1),
    #           cex = .8,
    #           col = c("#1b9e77","black"),
    #           legend = c("Mp","Mp + Mat M"), bty = "n")
    # }  

    if(mfg[2] == 1)
    {
      mtext(side = 2, text = "Natural\nMortality (/yr)", line = 3.5, cex = labcex)
      legend( x = "topleft", 
              lwd = c(2,2,2,NA),
              lty = c(1,4,2,NA),
              pch = c(22,NA,NA),
              pt.bg = c(polyCol,NA,NA), 
              pt.lwd = c(0,NA,NA),
              pt.cex = c(1.5,NA,NA),
              cex = .8,
              col = c("salmon","salmon","grey50"),
              legend = c("Age-2+ M","M0","Mbar"), bty = "n")
    }
}




# recruitments
plotRt <- function( repList = reports,
                    noPar = FALSE, 
                    pIdx = 1, 
                    labcex = .8,
                    plotLog = FALSE ,
                    heading = "dimLabs" )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear

  # Pull stuff from report
  Rt      <- report$R_pt[pIdx,]
  omegaRt <- report$omegaR_pt[pIdx,]
  sigmaR  <- report$sigmaR
  R0      <- report$R0_p[pIdx]

  ageData_agt <- repList$data$A_apgt[,pIdx,,]

  colRec_t <- rep("grey40", length(Rt))

  checkPos <- function( x )
  {
    if(any(x > 0))
      return(1)
    else return(0)
  }

  posAgeIdx_t <- apply( X = ageData_agt, FUN = checkPos, MARGIN = 3 )

  colRec_t[posAgeIdx_t == 0] <- "white"
  


  plotCI <- FALSE

  if(!is.null(repList$posts))
  {
    plotCI <- TRUE
    tInitModel_p <- repList$repOpt$tInitModel_p
    R_it <- repList$posts$R_ipt[,pIdx,]

    R_qt <- apply(X = R_it, FUN = quantile,
                    MARGIN = c(2), probs = c(0.025, 0.5, 0.975),
                    na.rm = TRUE )

    Rt <- apply(X = R_it, FUN = mean, MARGIN = c(2))

    R0 <- round(mean(repList$posts$R0_ip[,pIdx]),2)

  } else {
    R_qt <- 0
  }

  R0lab <- round(R0,2)

  Rt[Rt == 0] <- NA
  R_qt[R_qt == 0] <- NA

  if(plotLog)
  {
    R0 <- log(R0, base = 10)
    R_qt <- log(R_qt, base = 10)
    Rt <- log(Rt, base = 10)

  }

  # stock labels
  stockLabs <- dimnames(report$R_pt)[[1]]


  if( !is.null(heading))
  {
    if( heading == "dimLabs")
      stockLabel <- stockLabs
    else
      stockLabel <- try(heading[stockLabs])

    if( class(stockLabel) == "try-error" ) 
      stockLabel <- stockLabs
  }

  # Get number of time steps
  nT    <- report$nT
  years <- seq(from = initYear, length = nT, by = 1)

  vertLines <- seq(from = 1950, to = max(years), by = 10)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot recruitments
  plot( x = range(years), y = c(0,max(Rt,R_qt, na.rm =T) ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    mfg <- par("mfg")
    if(plotLog)
    {
      maxY <- ceiling(max(Rt,R_qt, na.rm =T))
      yTicks <- 0:maxY
      yLabs <- 10^yTicks
    }
    if( mfg[1] == mfg[3])
      axis(side = 1)
    if(!plotLog)
      axis(side = 2, las = 1  )
    if(plotLog)
      axis(side = 2, las = 1, at = yTicks, labels = yLabs )
    box()
    grid()
    # Plot recruitment
    # lines( x = years[1:nT], y = Rt[1:nT], lwd = 2, col = "grey40" )
    abline( h = R0, lty = 2, lwd = 2)
    if( plotCI )
      segments( x0 = years[1:nT],
                y0 = R_qt[1,1:nT],
                y1 = R_qt[3,1:nT],
                lwd = 1, col = "grey40" )
    points( x = years[1:nT], y = Rt[1:nT], lwd = 2, bg = "grey40",
            pch = 21, col = "grey40" )

    axis( side = 3, at = years[posAgeIdx_t == 0], labels = FALSE,
          col.ticks = "grey40", tck = .02, lwd.ticks = 3 )
    
    panLab( x = 0.8, y = 0.1, txt = paste("R0 = ", R0lab, sep = "") )
    panLab( x = 0.8, y = 0.2, txt = "|  Miss Ages")
    if( mfg[1] == mfg[3])
      mtext( side = 1, text = "Year", line = 2)
    if( mfg[2] == 1 )
      mtext(side = 2, text = "Recruits\n(millions)", line = 3.5, cex = labcex)
    if( mfg[1] == 1 & !is.null(heading) )
      mtext(side = 3, text = stockLabel[pIdx], line = 1, font = 2, cex = labcex )

}



# plotRule()
# Plot of Harvest control rule and recommended
# TAC based on SSPM estimate of biomass
# Input: repList = SISCAH reports object
# Output: NULL
# Side-effects: plots to active graphics device.
plotRule <- function( repList = reports )
{
  ctlList <- repList$ctlList
  repObj  <- repList$repOpt

  # Pull estimates
  C_t   <- repObj$C_pgt[1,1,]
  B_t   <- repObj$SB_pt[1,]
  B0    <- repObj$B0_p[1]

  # Control points
  ctlPts <- c(ctlList$mp$LCP,ctlList$mp$UCP)
  
  lowF   <- 0
  Fref   <- ctlList$mp$maxTHR
  inputF <- TRUE
  if( Fref == "est")
  {
    inputF <- FALSE
    Fref <- repObj$Fmsy
  }

  fYear <- 1951
  nT    <- repObj$nT
  year  <- fYear + nT

  Bseq <- seq(from = 0, to = B0, length.out = 100 )
  Fseq <- sapply( X = Bseq, FUN = .calcRampedHCR, 
                  LCP = ctlPts[1]*B0,
                  UCP = ctlPts[2]*B0,
                  Fref = Fref, lowFmult = lowF )

  predB <- B_t[length(B_t)]
  
  targF <- .calcRampedHCR(  B = predB/B0, 
                            LCP = ctlPts[1],
                            UCP = ctlPts[2],
                            Fref = Fref, lowFmult = lowF )
  Q     <- .calcLegalCatch( B = predB, 
                            LCP = ctlPts[1],
                            UCP = ctlPts[2],
                            Fref = Fref, lowFmult = lowF )

  Bhat    <- round( predB * 1000, digits=0 )
  estB0   <- trunc( B0*1000, digits=0 )
  targU   <- round( targF,   digits=4 )
  Q       <- trunc( Q*1000 )

  units <- " t"
  
  
  plot(x = c(0,B0), y=c(0,1.2*Fref), type = "n",
        xlab = "Estimated Biomass (kt)", ylab = "Harvest rate", las = 1)
    lines(x = Bseq, y = Fseq)
    abline( v = ctlPts * B0, lwd = 1, lty = 2 )
    points( x = predB, y = targF, pch = 16, col = "grey30" )

    txt <- bquote( hat(italic(B))[.(year)] == .(Bhat) )
    panLab( 0.75, 0.5, adj=0, cex=1, txt )
    txt <- bquote( hat(italic(B))[0] == .(estB0) )
    panLab( 0.75, 0.4, adj=0, cex=1, txt )
    txt <- bquote( italic(U)[.(year)] == .(targU) )
    panLab( 0.75, 0.3, adj=0, cex=1, txt )
    txt <- bquote( italic(Q)[.(year)] == .(Q) )
    panLab( 0.75, 0.2, adj=0, cex=1, txt )  

} # END plotRule()
