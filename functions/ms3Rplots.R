# <><><><><><><><><><><><><><><><><><><><><><><><><><><>
# ms3Rplots.R
#
# Plots for ms3R.
#
# Author: SDN Johnson
# Date: October 8, 2019
#
# Last Update: July 9, 2020
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><>


# plotStockRecruitRep()
# Beverton-Holt stock recruitment curve for an individual
# replicate with input ref pts picked out and 
# historical/projected recruitments shown. 
# inputs:
#   obj   = simulation object (blob)
#   iRep  = simulation replicate
# outputs:
#   no output object, but draws figure on plotting device
plotStockRecruitRep <- function( obj = blob, iRep = 1)
{

  # Pull B0, h, etc.
  B0    <- obj$rp[[iRep]]$B0_sp[1,1,1]
  R0    <- obj$rp[[iRep]]$R0_sp[1,1]
  h     <- obj$rp[[iRep]]$h_sp[1,1]
  rec.a <- obj$rp[[iRep]]$rec.a_sp[1,1]
  rec.b <- obj$rp[[iRep]]$rec.b_sp[1,1]

  if(!is.null(obj$ctlList$refPts))
  {
    inputB0   <- obj$ctlList$refPts$B0
    inputBmsy <- obj$ctlList$refPts$Bmsy
  } else {
    inputB0     <- B0
    inputBmsy   <- NA
  }

  SB_t  <- obj$om$SB_ispt[iRep,1,1,]
  R_t   <- obj$om$R_ispt[iRep,1,1,]

  nT  <- obj$om$nT
  tMP <- obj$om$tMP

  # First, make the curve
  Bseq <- seq(from =0, to = 2*B0, length.out = 100)
  Rseq <- rec.a * Bseq/(1 + rec.b*Bseq)

  inputRmsy   <- rec.a*inputBmsy/(1 + rec.b*inputBmsy)
  inputR0     <- rec.a*inputB0/(1 + rec.b*inputB0)
  
  plot( x = c(0,max(Bseq,1.1*SB_t)), y = c(0,max(Rseq, 1.1*R_t)), type = "n", las =1,
        xaxs = "i", yaxs = "i",
        xlab = "Spawning Biomass (kt)", ylab = "Age-1 Recruitment (1e6)" )
    lines( x = Bseq, y = Rseq, lwd = 3)
    points(x = SB_t[1:(tMP-2)], y = R_t[2:(tMP-1)], pch = 16, col = "grey30")
    points(x = SB_t[(tMP-1):(nT-1)], y = R_t[tMP:nT], pch = 16, col = "salmon")
    segments( x0 = c(inputB0,inputBmsy),
              y0 = 0, 
              y1 = c(inputR0,inputRmsy), lty = 2, col = c("grey40","darkgreen"), lwd = 2)
    segments( x0 = 0,
              x1 = c(inputB0,inputBmsy), 
              y0 = c(inputR0,inputRmsy), lty = 2, col = c("grey40","darkgreen"), lwd = 2)

    legend( x = "topright", bty = "n",
            legend = c( "Stock-Recruit Relationship",
                        "Unfished equilibrium",
                        "MSY equilibrium",
                        "Historical period",
                        "Projection period"),
            lty = c(1,2,2,NA,NA),
            lwd = c(3,2,2,NA,NA),
            pch = c(NA,NA,NA,16,16),
            col = c("black","grey40","darkgreen","grey30","salmon"))
}



refPtLHpairs <- function( refPt_sims = CC_sims )
{
  
  nReps <- length(refPt_sims$m1_i)



  colNames <- c("B0","h","m1","M0","Bmsy","Umsy","MSY")
  refPtTable <- array(NA, dim = c(nReps, length(colNames)))
  colnames(refPtTable) <- colNames
  refPtTable <- as.data.frame(refPtTable)


  refPtTable <- refPtTable

  for(i in 1:nReps)
  {
    eqbria <- solveSimEqbriaRep( empRefCurves = refPt_sims,
                                      maxXspline = 0.8,
                                      sIdx = 1, pIdx = 1,
                                      USR = 38,
                                      iRep = i )

    refPtTable[i, "B0" ]   <- eqbria$B0
    refPtTable[i, "m1" ]   <- refPt_sims$m1_i[i]
    refPtTable[i, "h" ]    <- refPt_sims$h_i[i]
    refPtTable[i, "M0" ]   <- refPt_sims$Mb_[i] + exp(-refPt_sims$m1_i[i])
    refPtTable[i, "Bmsy"]  <- eqbria$Bmsy
    refPtTable[i, "Umsy"]  <- eqbria$Umsy
    refPtTable[i, "MSY"]   <- eqbria$MSY

  }

  # Now we need to plot the pairs
  nPars <- length(colNames)
  par(mfrow = c(nPars,nPars), oma = c(5,5,3,3), mar = c(.1,.1,.1,.1))
  for( k in 1:nPars )
    for(j in 1:nPars)
    {
      if(k == j)
      {
        # plot diagonal
        parVec <- refPtTable[,k]
        histPar    <- hist(parVec, plot = FALSE)
        maxY <- 1.5*max(histPar$density)

        plot(x = range(parVec), y = c(0,maxY),
              type = "n", axes = FALSE)
          mfg <- par("mfg")
          # grid(lwd = 2, lty = 2)

          hist(parVec, add = TRUE, freq = FALSE, probability = TRUE,
                col = "grey60")

          legend(x = "top", legend = colNames[k], bty = "n")
      }

      if(k != j )
      {
        parVec1 <- refPtTable[,k] 
        parVec2 <- refPtTable[,j] 
        plot(x = range(parVec2), y = range(parVec1),
              type = "n", axes = FALSE)
          mfg <- par("mfg")
          grid(lwd = 2, lty = 2)

          points(x = parVec2, y = parVec1, pch = 16,
                  col = scales::alpha("grey40", alpha = 0.75))
      }

      if(mfg[2] == 1 & k != j )
      {
        axis(side = 2, las = 1)
        mtext(side = 2, text = colNames[k], line = 3)
      }

      if(mfg[2] == mfg[4] & k != j )
      {
        axis(side = 4, las = 1)
      }

      if(mfg[1] == 1 )
      {
        axis(side = 3 )
      }

      if(mfg[1] == mfg[3] )
      {
        axis(side = 1 )
        mtext(side = 1, text = colNames[j], line = 3)
      }

      box()
    }

}





plotHistProjProcError <- function(obj = blob,
                                iRep = 1,
                                sIdx = 1,pIdx = 1,
                                recDevOffset = 1)
{
  omegaR_t <- obj$om$errors$omegaR_ispt[iRep,sIdx,pIdx,]
  omegaM_t <- obj$om$errors$omegaM_ispt[iRep,sIdx,pIdx,]

  tMP <- obj$om$tMP
  nT  <- obj$om$nT

  fYear <- obj$ctlList$opMod$fYear

  yrs <- seq(from = fYear, by = 1, length.out = nT)

  par(  mfrow = c(2,1), 
        mar = c(.5,1,.5,1),
        oma = c(3,3,2,1))

  plot(x = range(yrs), y = range(omegaR_t),
        type = "n", xlab = "Year", ylab = "Std. Recruitment Deviations",
        axes = FALSE)
    axis(side = 2, las = 1)
    mtext(side = 2, text ="Std. R Deviations", line =2, cex = 2)
    grid()
    box()

    points(x = yrs[1:(tMP-recDevOffset)], y = omegaR_t[1:(tMP-recDevOffset)], col = "grey40", pch = 16, cex = 1.5)
    points(x = yrs[(tMP-recDevOffset+1):nT], y = omegaR_t[(tMP-recDevOffset+1):nT], col = "salmon", pch = 16, cex = 1.5)
    
    omegaR_t[1:(tMP-recDevOffset)] <- omegaR_t[1:(tMP-recDevOffset)] - 0.5*obj$om$sigmaR_isp[iRep,sIdx,pIdx]
    histMean <- mean(omegaR_t[1:(tMP-recDevOffset)])
    histSD <- sd(omegaR_t[1:(tMP-recDevOffset)])

    projMean <- mean(omegaR_t[(tMP-recDevOffset+1):nT])
    projSD <- sd(omegaR_t[(tMP-recDevOffset+1):nT])
    segments(x0 = yrs[1], x1 = yrs[tMP-recDevOffset], y0 = histMean, lwd = 2, col = "grey40")
    segments(x0 = yrs[1], x1 = yrs[tMP-recDevOffset], y0 = histMean+c(1,-1)*histSD, lwd = 2, col = "grey40", lty= 2)
    segments(x0 = yrs[tMP-recDevOffset+1], x1 = yrs[nT], y0 = projMean, lwd = 2, col = "salmon")
    segments(x0 = yrs[tMP-recDevOffset+1], x1 = yrs[nT], y0 = projMean+c(1,-1)*projSD, lwd = 2, col = "salmon", lty= 2)

  plot(x = range(yrs), y = range(omegaM_t),
        type = "n", xlab = "Year", ylab = "Std. M Deviations",
        axes = FALSE)
    axis(side = 2, las = 1)
    axis(side = 1)
    grid()
    box()
    mtext(side = 2, text ="Std. M Deviations", line =2, cex = 2)
    mtext(side = 1, text = "Year", line = 2, cex = 2)

    points(x = yrs[1:(tMP-1)], y = omegaM_t[1:(tMP-1)], col = "grey40", pch = 16, cex = 1.5)
    points(x = yrs[(tMP-1+1):nT], y = omegaM_t[(tMP-1+1):nT], col = "salmon", pch = 16, cex = 1.5)
    
    omegaM_t[1:(tMP-1)] <- omegaM_t[1:(tMP-1)] - 0.5*obj$om$sigmaM_isp[iRep,sIdx,pIdx]
    histMean <- mean(omegaM_t[1:(tMP-1)])
    histSD <- sd(omegaM_t[1:(tMP-1)])



    projMean <- mean(omegaM_t[(tMP-1+1):nT])
    projSD <- sd(omegaM_t[(tMP-1+1):nT])
    segments(x0 = yrs[1], x1 = yrs[tMP-1], y0 = histMean, lwd = 2, col = "grey40")
    segments(x0 = yrs[1], x1 = yrs[tMP-1], y0 = histMean+c(1,-1)*histSD, lwd = 2, col = "grey40", lty= 2)
    segments(x0 = yrs[tMP-1+1], x1 = yrs[nT], y0 = projMean, lwd = 2, col = "salmon")
    segments(x0 = yrs[tMP-1+1], x1 = yrs[nT], y0 = projMean+c(1,-1)*projSD, lwd = 2, col = "salmon", lty= 2)
}

# genSimInfo()
genSimInfo <- function(simFolder = here('mserproject'))
{

  # Read in sims
  sims <- list.files(file.path(simFolder))
  sims <- sims[grepl("sim",sims)]

  readInfoFile <- function( sim )
  {
    infoPath <- file.path(simFolder,sim,'infoFile.txt' ) 
    info <- lisread(infoPath)
    info.df <- as.data.frame(info)
    info.df$simLabel <- sim

    info.df
  }

  # Read in info files, sort by  scenarios
  info.df <- lapply( X = sims, FUN = readInfoFile )
  info.df <- do.call( "rbind", info.df )

}  #END genSimInfo()


plotMtBt <- function( obj = blob,
                      iRep = 1 )
{

  # Check if depM
  opMod     <- obj$ctlList$opMod
  repObj    <- obj$ctlList$opMod$histRpt
  depM      <- repObj$densityDepM
  juveMage  <- repObj$juveMage

  # Dimensions
  nS <- obj$om$nS 
  nP <- obj$om$nP
  nT <- obj$om$nT

  M_spt <- array(0,dim = c(nS,nP,nT))
  B_spt <- array(0,dim = c(nS,nP,nT))

  # Get Bt and Mt series
  M_spt[1:nS,1:nP,] <- obj$om$M_iaxspt[iRep,juveMage+1,1,,,]
  B_spt[1:nS,,]     <- obj$om$B_ispt[iRep,,,]

  # Pull m1 and Mb values
  Mb_sp             <- obj$rp[[iRep]]$Mb_sp
  m1_sp             <- obj$rp[[iRep]]$m1_sp

  # Pull totB0
  totB0_sp <- obj$rp[[iRep]]$totB0_sp

  Bseq    <- seq(from = 1, to = 1.2*max(B_spt,na.rm = T), length.out = 1000)
  Mseq_sp <- array(0,dim = c(nS,nP,1000))

  M0_sp   <- array(0,dim = c(nS,nP))

  for(s in 1:nS )
    for(p in 1:nP )
    {
      Mseq_sp[s,p,] <- Mb_sp[s,p] + exp(-m1_sp[s,p] * Bseq/totB0_sp[s,p])
      M0_sp[s,p]    <- Mb_sp[s,p] + exp(-m1_sp[s,p])
    }

  par(mfrow = c(nS,nP), 
      oma   = c(3,3,1,1),
      mar   = c(1,1,1,1) )
  for( s in 1:nS )
    for( p in 1:nP )
    {
      plot( x = c(0,1.2*max(B_spt,na.rm = T)),
            y = c(0,1.2*max(M_spt,na.rm = T)),
            axes = FALSE, type = "n", xaxs = "i", yaxs = "i")
      axis( side = 1)
      axis( side = 2, las = 1)
      grid()
      box()
      lines(x  = Bseq, y = Mseq_sp[s,p,], col = "salmon", lwd = 3)
      points( x = B_spt[s,p,], y = M_spt[s,p,], pch = 16 )
      segments( x0 = totB0_sp[s,p], y0 = 0, y1 = M0_sp[s,p], lty = 2)
      segments( x0 = 0, x1 = totB0_sp[s,p], y0 = M0_sp[s,p], lty = 2)
    }

  mtext( side = 1, outer = T, text = "Age 2+ Biomass (kt)", line = 2)
  mtext( side = 2, outer = T, text = "Natural Mortalty (/yr)", line = 2)

}



# multiSimBtIt()
# Plot BtIt for same rep across different mp/scenario combinations
multiSimBtIt <- function( iRep = 1,
                           simFolder,
                           OMs = c("rw5_loPondM, rw15_hiPondM"),
                           MPs = c("NoFish"),
                           yLimB     = c(0,22))


{


  info.df <- genSimInfo(simFolder=simFolder) %>%
              mutate_if( is.factor, as.character) %>%
              filter( scenario %in% OMs,
                      mp %in% MPs )
  
  sims    <- expand.grid(OMs, MPs)
  names(sims) <- c('OM','MP')
        

  # Set up plot environment 
  par(mfcol = c(3,nrow(info.df) ), mar = c(1.5,1.5,0.5,1), 
      oma = c( 3,4,1,0))

  # Loop over MPs
  for( i in 1:nrow(info.df) )
  {
    subDF <-  info.df %>%
              filter( mp == sims$MP[i] &
                      scenario == sims$OM[i])

    simID <- subDF[1,"simLabel"]

    simFile <- paste(simID,".RData",sep = "")
    simPath <- file.path(simFolder,simID,simFile)

    # Load blob
    load(simPath)
    
    # Plot BtIt
    plotBtIt_p(obj = blob, iRep = 1, f=5, addCatch=TRUE, 
               parArg=FALSE, YlabOn=FALSE, legdOn=FALSE)

    mtext(side = 3, text = paste(sims$OM[i],sims$MP[i]) , line = 33, cex = 1)

  }

  mtext( side =2, outer = TRUE, text = "Spawning biomass kt)", line = 1)


            
} # END multiSimBtCtRt()

plotMultiHCR <- function(LCPs, UCPs, hcrNames, xlab)
{
  # number of plots
  nI <- length(LCPs)

  par(mfrow=c(1,nI), mgp=c(1.4,0.6,0), mar=c(2,3,1,1), oma = c(2,1,0,0))
  for(i in 1:nI)
    plotHockeyStickHCR(LCP = LCPs[i], UCP = UCPs[i], 
                        hcrName=hcrNames[i])

  mtext( side = 1, text = xlab, outer = TRUE)

}

# plotHockeyStickHCR
plotHockeyStickHCR <- function( LCP = .5, UCP = .6, refHR = .1, 
                                refHRaxis = NULL, hcrName='',
                                refHRlab='', mpLab=NULL,
                                refB = 1, yAXT='n',
                                yLim = NULL,
                                xLab='', yLab='Harvest Rate')
                                # xLab = expression(paste("Stock Status (",B[t]/B[MSY],")",sep = "")),
                                # yLab = "Target Fishing Mortality (F)" )
{
  if(is.null(yLim))
    yLim <- c(0,1.5*refHR)


  plot( x = c(0,refB*1.2), y = yLim, type = "n", xlab = "", ylab = "", las = 1,
        yaxt=yAXT, main=hcrName )
    segments( x0 = 0, x1 = LCP * refB, y0 = 0, lwd = 1.5 )
    segments( x0 = LCP*refB, x1 = UCP*refB, y0 = 0, y1 = refHR, lwd = 1.5 )
    segments( x0 = UCP*refB, x1 = refB*1.2, y0 = refHR, y1 = refHR, lwd = 1.5 )
    abline( v = c(UCP*refB, LCP * refB), lwd = .8, lty = 2)
    abline( h = refHR, lty = 3, lwd = .8)
    mtext( side = 1, text = xLab, line = 2 )
    mtext( side = 2, text = yLab, line = 2 )

    if(is.null(refHRaxis))
      refHRaxis = refHR

    if(!is.null(mpLab))
      mtext( side = 4, text = mpLab, line = 0.5, cex=0.7 )


    axis(side=2, at=c(0,refHR), tick=TRUE, labels=c('0',refHRaxis), las=1)
    text(x=0.2,y=refHR+0.004, refHRlab)


} # END plotHockeyStickHCR()

# plotTACallocationError()
# A multi-panel for understanding the
# level of error in TAC allocations under 
# aggregate MPs. Not sure of the best way
# to look at this...
plotTACallocationError <- function( obj = blob,
                                    iRep = 1 )
{
  goodReps <- obj$goodReps

  # Get biomass arrays
  SB_ispt        <- obj$om$SB_ispt[goodReps,,,]
  VB_ispt        <- obj$om$vB_ispft[goodReps,,,2,]
  totB_ispt      <- obj$om$B_ispt[goodReps,,,]
  retroSB_itspt  <- obj$mp$assess$retroSB_itspt[goodReps,,,,]

  ctlList <- obj$ctlList

  retroSB_itspt[retroSB_itspt < 0] <- NA
  nReps     <- sum(goodReps)

  # Model dims
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT



}



makeAMEstDists <- function( obj = blob,
                            sIdx = 1, pIdx = 1,
                            pt = 1 )
{
  # Pull model dims
  nS    <- obj$om$nS
  nP    <- obj$om$nP
  nX    <- obj$om$nX
  nF    <- obj$om$nF
  tMP   <- obj$om$tMP
  nT    <- obj$om$nT

  # Get goodreps
  goodReps <- which(obj$goodReps)
  nReps    <- length(goodReps)

  repObj  <- obj$ctlList$opMod$histRpt
  
  # Biological pars
  B0_i      <- rep(0,nReps)
  R0_i      <- rep(0,nReps)
  totB0_i   <- rep(0,nReps)
  h_i       <- rep(0,nReps)
  M_i       <- array(0,dim = c(nReps))
  M0_i      <- array(0,dim = c(nReps))
  m1_i      <- rep(0,nReps)
  M_it      <- array(0,dim = c(nReps,tMP - 1))

  # Nuisance pars
  q_if      <- array(0,dim = c(nReps,nF))
  tauObs_if <- array(0,dim = c(nReps,nF))

  # Selectivity parameters
  SelA_ift    <- array(0,dim = c(nReps,nF,tMP - 1))
  SelB_ift    <- array(0,dim = c(nReps,nF,tMP - 1))

  estSelA_ift    <- array(0,dim = c(nReps,nF,tMP - 1))
  estSelB_ift    <- array(0,dim = c(nReps,nF,tMP - 1))
  

  
  
  B0_i     <- obj$mp$assess$retroB0_itsp[goodReps,pt,sIdx,pIdx]
  totB0_i  <- obj$mp$assess$retrototB0_itsp[goodReps,pt,sIdx,pIdx]
  R0_i     <- obj$mp$assess$retroR0_itsp[goodReps,pt,sIdx,pIdx]
  h_i      <- obj$mp$assess$retroh_itsp[goodReps,pt,sIdx,pIdx]
  M_i      <- obj$mp$assess$retroM_itsp[goodReps,pt,sIdx,pIdx]
  m1_i     <- obj$mp$assess$retrom1_itsp[goodReps,pt,sIdx,pIdx]
  M0_i     <- M_i + exp(-m1_i)
  q_if     <- obj$mp$assess$retroq_itspf[goodReps,pt,sIdx,pIdx,]
  
  # And now adjust the estimates as well
  SelA_ift[1:nReps,,]   <- obj$mp$assess$retroSelA_itspft[goodReps,pt,sIdx,pIdx,,1:(tMP-1)]
  SelB_ift[1:nReps,,]   <- obj$mp$assess$retroSelA_itspft[goodReps,pt,sIdx,pIdx,,1:(tMP-1)] + obj$mp$assess$retroSelB_itspft[goodReps,pt,sIdx,pIdx,,1:(tMP-1)]

    
  
  B0_q   <- quantile(B0_i, probs = c(0.025,0.5,0.975),na.rm = T)
  totB0_q<- quantile(totB0_i, probs = c(0.025,0.5,0.975),na.rm = T)
  R0_q   <- quantile(R0_i, probs = c(0.025,0.5,0.975),na.rm = T)
  h_q    <- quantile(h_i, probs = c(0.025,0.5,0.975), na.rm = T)
  M_q    <- quantile(M_i, probs =c(0.025,0.5,0.975), na.rm = T)
  M0_q   <- quantile(M0_i, probs =c(0.025,0.5,0.975), na.rm = T)
  m1_q   <- quantile(m1_i, probs =c(0.025,0.5,0.975), na.rm = T)
  q_qf   <- apply(X = q_if, FUN = quantile, probs =c(0.025,0.5,0.975), MARGIN = 2, na.rm = T)
  tau_qf <- apply(X = tauObs_if, FUN = quantile, probs =c(0.025,0.5,0.975), MARGIN = 2, na.rm = T)

  SelA_qft <- apply(X = SelA_ift, FUN = quantile, probs =c(0.025,0.5,0.975), MARGIN = 2:3, na.rm = T)
  SelB_qft <- apply(X = SelB_ift, FUN = quantile, probs =c(0.025,0.5,0.975), MARGIN = 2:3, na.rm = T)

  distList <- list( B0_i       = B0_i,
                    totB0_i    = totB0_i,
                    R0_i       = R0_i,
                    M0_i       = M0_i,
                    h_i        = h_i, 
                    M_i        = M_i,
                    m1_i       = m1_i,
                    q_if       = q_if, 
                    tau_if     = tauObs_if,
                    SelA_ift   = SelA_ift,
                    SelB_ift   = SelB_ift )

  quantList <- list(  B0_q       = B0_q,
                      totB0_q    = totB0_q,
                      R0_q       = R0_q, 
                      M0_q       = M0_q, 
                      h_q        = h_q, 
                      M_q        = M_q,
                      m1_q       = m1_q,
                      q_qf       = q_qf, 
                      tau_qf     = tau_qf,
                      SelA_qft   = SelA_qft,
                      SelB_qft   = SelB_qft )

  outList <- list(  distList  = distList,
                    quantList = quantList )

  outList
} # END makeAMEstDists


makeAMREs <- function(  obj = blob,
                        sIdx = 1, pIdx = 1,
                        pt = 1 )
{
  # Pull model dims
  nS    <- obj$om$nS
  nP    <- obj$om$nP
  nX    <- obj$om$nX
  nF    <- obj$om$nF
  tMP   <- obj$om$tMP
  nT    <- obj$om$nT

  # Get goodreps
  goodReps <- which(obj$goodReps)
  nReps    <- length(goodReps)

  repObj  <- obj$ctlList$opMod$histRpt
  
  # Biological pars
  B0_i      <- rep(0,nReps)
  R0_i      <- rep(0,nReps)
  totB0_i   <- rep(0,nReps)
  h_i       <- rep(0,nReps)
  M_i       <- array(0,dim = c(nReps))
  M0_i      <- array(0,dim = c(nReps))
  m1_i      <- rep(0,nReps)
  M_it      <- array(0,dim = c(nReps,tMP - 1))

  # Nuisance pars
  q_if      <- array(0,dim = c(nReps,nF))
  tauObs_if <- array(0,dim = c(nReps,nF))

  # Selectivity parameters
  SelA_ift    <- array(0,dim = c(nReps,nF,tMP - 1))
  SelB_ift    <- array(0,dim = c(nReps,nF,tMP - 1))

  estSelA_ift    <- array(0,dim = c(nReps,nF,tMP - 1))
  estSelB_ift    <- array(0,dim = c(nReps,nF,tMP - 1))
  

  
  # Pull true values
  for( idx in 1:length(goodReps) )
  {
    i <- goodReps[idx]

    B0_i[idx]       <- obj$rp[[i]]$B0_sp[sIdx,pIdx,1]
    R0_i[idx]       <- obj$rp[[i]]$R0_sp[sIdx,pIdx]
    totB0_i[idx]    <- obj$rp[[i]]$totB0_sp[sIdx,pIdx]
    h_i[idx]        <- repObj$rSteepness_p[pIdx]
    M_i[idx]        <- obj$rp[[i]]$Mb_sp[sIdx,pIdx]
    m1_i[idx]       <- obj$rp[[i]]$m1_sp[sIdx,pIdx]
    M0_i[idx]       <- M_i[idx] + exp(-m1_i[idx])
    q_if[idx,]      <- repObj$qComb_pg[pIdx,]
    tauObs_if[idx,] <- repObj$tauComb_pg[pIdx,]
    
    SelA_ift[idx,,1:(tMP-1)]   <- repObj$SelAlpha_pgt[1,,]
    SelB_ift[idx,,1:(tMP-1)]   <- SelA_ift[idx,,1:(tMP-1)] + repObj$SelBeta_pgt[1,,]
    
  }

  # And now adjust the estimates as well
  estSelA_ift[1:nReps,,]   <- obj$mp$assess$retroSelA_itspft[goodReps,pt,sIdx,pIdx,,1:(tMP-1)]
  estSelB_ift[1:nReps,,]   <- obj$mp$assess$retroSelA_itspft[goodReps,pt,sIdx,pIdx,,1:(tMP-1)] + obj$mp$assess$retroSelB_itspft[goodReps,pt,sIdx,pIdx,,1:(tMP-1)]

    
  if(!is.null(obj$mp$assess$retroM0_itsp))
    errM0_i     <- -1*(M0_i - obj$mp$assess$retroM0_itsp[goodReps,pt,sIdx,pIdx])/M0_i

  errB0_i     <- -1*(B0_i - obj$mp$assess$retroB0_itsp[goodReps,pt,sIdx,pIdx])/B0_i
  errtotB0_i  <- -1*(totB0_i - obj$mp$assess$retrototB0_itsp[goodReps,pt,sIdx,pIdx])/totB0_i
  errR0_i     <- -1*(R0_i - obj$mp$assess$retroR0_itsp[goodReps,pt,sIdx,pIdx])/R0_i
  errh_i      <- -1*(h_i - obj$mp$assess$retroh_itsp[goodReps,pt,sIdx,pIdx])/h_i
  errM_i      <- -1*(M_i - obj$mp$assess$retroM_itsp[goodReps,pt,sIdx,pIdx])/M_i
  errm1_i     <- -1*(m1_i - obj$mp$assess$retrom1_itsp[goodReps,pt,sIdx,pIdx])/m1_i
  errq_if     <- -1*(q_if - obj$mp$assess$retroq_itspf[goodReps,pt,sIdx,pIdx,])/q_if
  errtau_if   <- -1*(tauObs_if - obj$mp$assess$retrotauObs_itspf[goodReps,pt,sIdx,pIdx,])/tauObs_if

  if(is.null(obj$mp$assess$retroM0_itsp))
  {
    retroM0_itsp <- obj$mp$assess$retroM_itsp + exp(-1 * obj$mp$assess$retrom1_itsp)
    errM0_i       <- -1*(M0_i - retroM0_itsp[goodReps,pt,sIdx,pIdx])/M0_i
  }

  # Get errors for historical period only - will update to
  # future once we have tv sel in the projections. Right now,
  # we're good as all tvSel happens in the history and the
  # historical estimates for the newRV will be captured here.
  errSelA_ift   <- -1*(SelA_ift - estSelA_ift)/SelA_ift
  errSelB_ift   <- -1*(SelB_ift - estSelB_ift)/SelB_ift
  
  errB0_q   <- quantile(errB0_i, probs = c(0.025,0.5,0.975),na.rm = T)
  errtotB0_q<- quantile(errtotB0_i, probs = c(0.025,0.5,0.975),na.rm = T)
  errR0_q   <- quantile(errR0_i, probs = c(0.025,0.5,0.975),na.rm = T)
  errh_q    <- quantile(errh_i, probs = c(0.025,0.5,0.975), na.rm = T)
  errM_q    <- quantile(errM_i, probs =c(0.025,0.5,0.975), na.rm = T)
  errM0_q   <- quantile(errM0_i, probs =c(0.025,0.5,0.975), na.rm = T)
  errm1_q   <- quantile(errm1_i, probs =c(0.025,0.5,0.975), na.rm = T)
  errq_qf   <- apply(X = errq_if, FUN = quantile, probs =c(0.025,0.5,0.975), MARGIN = 2, na.rm = T)
  errtau_qf <- apply(X = errtau_if, FUN = quantile, probs =c(0.025,0.5,0.975), MARGIN = 2, na.rm = T)

  errSelA_qft <- apply(X = errSelA_ift, FUN = quantile, probs =c(0.025,0.5,0.975), MARGIN = 2:3, na.rm = T)
  errSelB_qft <- apply(X = errSelB_ift, FUN = quantile, probs =c(0.025,0.5,0.975), MARGIN = 2:3, na.rm = T)

  errList <- list(  errB0_i       = errB0_i,
                    errtotB0_i    = errtotB0_i,
                    errR0_i       = errR0_i,
                    errM0_i       = errM0_i,
                    errh_i        = errh_i, 
                    errM_i        = errM_i,
                    errm1_i       = errm1_i,
                    errq_if       = errq_if, 
                    errtau_if     = errtau_if,
                    errSelA_ift   = errSelA_ift,
                    errSelB_ift   = errSelB_ift )

  quantList <- list(  errB0_q       = errB0_q,
                      errtotB0_q    = errtotB0_q,
                      errR0_q       = errR0_q, 
                      errM0_q       = errM0_q, 
                      errh_q        = errh_q, 
                      errM_q        = errM_q,
                      errm1_q       = errm1_q,
                      errq_qf       = errq_qf, 
                      errtau_qf     = errtau_qf,
                      errSelA_qft   = errSelA_qft,
                      errSelB_qft   = errSelB_qft )

  outList <- list(  errList = errList,
                    quantList = quantList )

  outList
} # END makeAMREs


# plotAMREs()
plotAMREs <- function(  obj = blob,
                        sIdx = 1, pIdx = 1,
                        pt = 1, plotFleets = c(4:5),
                        fleetNames = c("Surface","Dive") )
{
  # Pull model dims
  nReps <- obj$nSims
  nS    <- obj$om$nS
  nP    <- obj$om$nP
  nX    <- obj$om$nX
  nF    <- obj$om$nF
  tMP   <- obj$om$tMP
  nT    <- obj$om$nT

  repObj  <- obj$ctlList$opMod$histRpt

  
  quantErrors <- makeAMREs( obj = obj, sIdx = sIdx,
                            pIdx = pIdx, pt = pt)$quantList

  errB0_q       <- quantErrors$errB0_q
  errtotB0_q    <- quantErrors$errtotB0_q
  errR0_q       <- quantErrors$errR0_q
  errM0_q       <- quantErrors$errM0_q
  errh_q        <- quantErrors$errh_q
  errM_q        <- quantErrors$errM_q
  errm1_q       <- quantErrors$errm1_q
  errq_qf       <- quantErrors$errq_qf
  errtau_qf     <- quantErrors$errtau_qf
  errSelA_qft   <- quantErrors$errSelA_qft
  errSelB_qft   <- quantErrors$errSelB_qft
  
  labels <- c("B0","totB0","R0","h","Mb","m1","M0")
  
  nErr <- length(labels) + 2*length(plotFleets)

  
  
  for(f in 1:length(plotFleets))
    labels <- c(labels,paste0("q_",fleetNames[f]),paste0("tau_",fleetNames[f]))

  plot( x = c(1,nErr), xlab = "", ylab = "",
        y = c(-1,1), axes = FALSE, type = "n" )
    axis( side = 1, at = 1:nErr, labels = labels )
    axis(side = 2, las =1)
    grid()
    box()

    # B0 error
    segments(x0 = 1, y0 = errB0_q[1], y1 = errB0_q[3], lwd = 3, col = "grey60" )
    points(x = 1, y = errB0_q[2], pch = 16, col = "grey60" )

    # totB0 error
    segments(x0 = 2, y0 = errtotB0_q[1], y1 = errtotB0_q[3], lwd = 3, col = "grey60" )
    points(x = 2, y = errtotB0_q[2], pch = 16, col = "grey60" )


    # R0 error
    segments(x0 = 3, y0 = errR0_q[1], y1 = errR0_q[3], lwd = 3, col = "grey60" )
    points(x = 3, y = errR0_q[2], pch = 16, col = "grey60" )


    # h error
    segments(x0 = 4, y0 = errh_q[1], y1 = errh_q[3], lwd = 3, col = "grey60" )
    points(x = 4, y = errh_q[2], pch = 16, col = "grey60" )    

    # Mb error
    segments(x0 = 5, y0 = errM_q[1], y1 = errM_q[3], lwd = 3, col = "grey60" )
    points(x = 5, y = errM_q[2], pch = 16, col = "grey60" )    

    # m1 error
    segments(x0 = 6, y0 = errm1_q[1], y1 = errm1_q[3], lwd = 3, col = "grey60" )
    points(x = 6, y = errm1_q[2], pch = 16, col = "grey60" )    

    # M0 error
    segments(x0 = 7, y0 = errM0_q[1], y1 = errM0_q[3], lwd = 3, col = "grey60" )
    points(x = 7, y = errM0_q[2], pch = 16, col = "grey60" )    

    nPlotFleets <- length(plotFleets)
    # Catchability error
    for(j in 1:length(plotFleets))
    {
      fIdx <- plotFleets[j]
      segments(x0 = 8 + (j-1) * 2, y0 = errq_qf[1,fIdx], y1 = errq_qf[3,fIdx], lwd = 3, col = "grey60" )
      points(x = 8 + (j-1) * 2, y = errq_qf[2,fIdx], pch = 16, col = "grey60" )    

      segments(x0 = 9 + (j-1) * 2, y0 = errtau_qf[1,fIdx], y1 = errtau_qf[3,fIdx], lwd = 3, col = "grey60" )
      points(x = 9 + (j-1) * 2, y = errtau_qf[2,fIdx], pch = 16, col = "grey60" )    
    }

    
    mtext( side = 2, text = "Relative Error in AM estimates", line = 3)
}


# plotAMREsel()
# Relative error distributions for selectivity
# parameters at a given projection time-step
plotAMREsel <- function(  obj = blob,
                          sIdx = 1, pIdx = 1,
                          pt = 1 )
{
  # Pull model dims
  nReps <- obj$nSims
  nS    <- obj$om$nS
  nP    <- obj$om$nP
  nX    <- obj$om$nX
  nF    <- obj$om$nF
  tMP   <- obj$om$tMP
  nT    <- obj$om$nT

  repObj  <- obj$ctlList$opMod$histRpt

  # Get fleet names
  fleetNames_f <- dimnames(repObj$C_pgt)[[2]]

  # Make error distributions
  quantErrors <- makeAMREs( obj = obj, sIdx = sIdx,
                            pIdx = pIdx, pt = pt)$quantList

  errList <- list()
  errList$errSelA_qft   <- quantErrors$errSelA_qft
  errList$errSelB_qft   <- quantErrors$errSelB_qft

  parNames <- c("L50A","L95A")

  fleetCols <- RColorBrewer::brewer.pal(nF,"Paired")

  if(pt > 4)
    nFleets <- nF
  else nFleets <- nF-1
  
  par(mfrow = c(4,1), oma = c(3,5,3,3), mar = c(.1,.1,.1,.1))
  for( j in 1:4)
  {
    # maxY <- max(abs(errList[[j]]),na.rm =T)
    plot( x = c(0,nFleets+1), xlab = "", ylab = "",
          y = c(-1,1), axes = FALSE, type = "n" )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3])
        axis( side = 1, at = 1:(nFleets), labels = fleetNames_f[1:(nFleets)] )
      
      abline(h = 0, lty = 2, lwd = 1)

      axis(side = 2, las =1)
      legend(x = "topright", legend = parNames[j], bty= "n")
      grid()
      box()

      # Now, loop over fleets and plot the quantiles      
      for( f in 1:(nFleets))
      {
        if( (j %in% c(3,4)) & (!f %in% c(1,7,8,9)) )
          next
        
        segments( x0 = f, y0 = errList[[j]][1,f,1], y1 = errList[[j]][3,f,1], col = fleetCols[f],
                  lwd = 3 )
        points( x = f, y = errList[[j]][2,f,1], col = fleetCols[f],
                pch = 16, cex = 2 )
      }
  }

      
  mtext( side = 2, text = "Relative Error in AM estimates", line = 3, outer = TRUE)
} # END plotAMREsel()


# plotDDmort()
# Plots density dependent mortality relationship
plotDDmort <- function( obj = blob,
                        iRep = 1)
{
  
  opMod     <- obj$ctlList$opMod
  repObj    <- obj$ctlList$opMod$histRpt
  depM      <- repObj$densityDepM
  juveMage  <- repObj$juveMage
  

  # Dimensions
  nA <- obj$om$nA
  nS <- obj$om$nS 
  nP <- obj$om$nP
  nT <- obj$om$nT
  tMP <- obj$om$tMP

  fYear <- repObj$fYear
  yrs   <- seq( from = fYear, by = 1, length.out = nT )
  
  M_aspt <- array(0,dim = c(nA,nS,nP,nT))
  M_spt <- array(0,dim = c(nS,nP,nT))
  B_spt <- array(0,dim = c(nS,nP,nT))
  
  # Get Bt and Mt series
  M_aspt[1:nA,1:nS,1:nP,] <- obj$om$M_iaxspt[iRep,,1,,,]
  M_spt[1:nS,1:nP,] <- obj$om$M_iaxspt[iRep,juveMage+1,1,,,]
  B_spt[1:nS,,]     <- obj$om$B_ispt[iRep,,,]
  
  # Pull m1 and Mb values
  Mb_sp             <- obj$rp[[iRep]]$Mb_sp
  m1_sp             <- obj$rp[[iRep]]$m1_sp
  
  # Pull totB0
  totB0_sp <- obj$rp[[iRep]]$totB0_sp
  
  M0_sp   <- array(0,dim = c(nS,nP))
  
  for(s in 1:nS )
    for(p in 1:nP )
    {
      M0_sp[s,p]    <- Mb_sp[s,p] + exp(-m1_sp[s,p])
    }
  
  
  D_spt <- B_spt
  for(s in 1:nS )
    for(p in 1:nP )
      D_spt[s,p,] <- B_spt[s,p,]/totB0_sp[s,p]
  
  
  
  
  Dseq <- seq(from = 0, to = max(D_spt), length.out = 1000 )
  
  M_spb <- array(0, dim = c(nS,nP,length(Dseq)) )
  
  for( s in 1:nS)
    for( p in 1:nP )
    {
      M_spb[s,p,] <- Mb_sp[s,p] + exp(-m1_sp[s,p] * Dseq )
      
    }
  
  colFunc <- colorRampPalette(c("red","royalblue"))
  yearCols <- colFunc(dim(D_spt)[3])
  
  
  par(mfrow = c(nP,1), mar = c(.1,1,.1,1), oma = c(3,4,1,1) )
  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = c(0,max(D_spt)), y = c(0,max(M_spb[s,p,])), type = "n",
            axes = FALSE, xaxs = "i", yaxs = "i" )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis(side = 1)
      axis( side = 2, las = 1)
      grid()
      box()
      lines( x = Dseq, y = M_spb[s,p,], lwd = 2, col = "salmon" )
      points( x = D_spt[s,p,1:(tMP-1)], y = M_aspt[juveMage+1,s,p,1:(tMP-1)], pch = 16, col = yearCols[1:(tMP-1)] )
      points( x = D_spt[s,p,tMP:nT], y = M_aspt[juveMage+1,s,p,tMP:nT], pch = 17 , col=yearCols[tMP:nT])
      segments(x0 = 1, y0 = 0, y1 = M0_sp[s,p], lty = 2 )
      segments(x0 = 0, x1 = 1, y0 = M0_sp[s,p], lty = 2 )
      
      legend( x = "topright", bty = "n",
              legend = c( "Expected depM curve",
                          paste0("M ", fYear," - ", yrs[tMP]),
                          paste0("M ", yrs[tMP + 1]," - ", yrs[nT]),
                          fYear, yrs[nT],
                          paste0("M_b = ", round(Mb_sp[s,p],2)),
                          paste0("m1 = ", round(m1_sp[s,p],2))),
              lty = c(1,NA,NA,NA,NA),
              pch = c(NA,16,17,15,15,NA,NA),
              col = c("red","grey40","grey40",yearCols[1],yearCols[nT],NA,NA))
    }
  }
  mtext( side = 1, text = "Age 2+ biomass depletion", line = 2, outer = TRUE)
  mtext( side = 2, text = "Age 2+ mortality (/yr)", line = 2, outer = TRUE)
  
} # END plotDDmort()

plotMultiRetroErrDist <- function(  blobList = list(blob),
                                    sIdx = 1, pIdx = 1,
                                    pt = 1, yrange = 1,
                                    plotFleets = c(1:3),
                                    fleetNames = c( "Reduction",
                                                    "Seine-roe",
                                                    "Gillnet-roe"),
                                    blobLabs = c( "simIdx",
                                                  "simAge",
                                                  "simAll") )
{
  retroBioDist <- lapply(X = blobList, FUN = makeAMbioDists,
                          sIdx = sIdx, pIdx = pIdx, pt = pt )

  nReps <- blobList[[1]]$nSims
  nS    <- blobList[[1]]$om$nS
  nP    <- blobList[[1]]$om$nP
  nX    <- blobList[[1]]$om$nX
  nF    <- blobList[[1]]$om$nF
  nT    <- blobList[[1]]$om$nT
  tMP   <- blobList[[1]]$om$tMP

  nRows <- length(plotFleets)+1
  nBlobs <- length(blobList)

  blobCols <- RColorBrewer::brewer.pal(nBlobs,"Dark2")


  yrs <- seq(from = 1968, by = 1, length.out = tMP+pt-2)
  nYrs <- length(yrs)

  par(  mfcol = c(nRows,nBlobs), oma = c(5,5,1,1),
        mar = c(0.1,0.1,.1,.1) )

  # Plot spawning biomass
  for( k in 1:nBlobs)
  {
    plot( x = range(yrs), 
          y = c(-yrange,yrange),
          xlab = "", ylab = "", las = 1, type = "n",
          axes = FALSE )
      mfg <- par("mfg")
      if(mfg[2] == 1)
        axis(side = 2, las =1)
      abline(h = 0, lty =2, lwd = .8)
      box()
      if(k == 1)
        legend(x = "topleft", legend = "Spawning Biomass", bty = "n")
      
      mtext(side = 3, text  = blobLabs[k], font = 2)
      
      sbErr_qt <- retroBioDist[[k]]$errSB_qt
      polygon( x = c(yrs,rev(yrs)),
               y = c(sbErr_qt[1,1:nYrs],sbErr_qt[3,nYrs:1]),
                col = scales::alpha(blobCols[k],.5), border = NA)
      lines( x = yrs, y = sbErr_qt[2,1:nYrs], col = blobCols[k],
              lty = k, lwd = 3 ) 
        # lines( x = yrs, y = sbErr_qt[1,1:nYrs], col = blobCols[k],
        #         lty = 2, lwd = 1 ) 
        # lines( x = yrs, y = sbErr_qt[3,1:nYrs], col = blobCols[k],
        #         lty = 2, lwd = 1 ) 

    for( fIdx in 1:length(plotFleets) )
    {
      f <- plotFleets[fIdx]
      plot( x = range(yrs), 
          y = c(-yrange,yrange),
          xlab = "", ylab = "", las = 1, type = "n", axes = FALSE )
      
      abline(h = 0, lty =2, lwd = .8)
      box()
      mfg <- par("mfg")
      if( mfg[1] == mfg[3])
        axis( side = 1)
      if(mfg[2] == 1)
        axis(side = 2, las =1)
      if(k == 1)
        legend( x = "topleft", bty = "n",
              legend = fleetNames[fIdx])
      
      vbErr_qft <- retroBioDist[[k]]$errVB_qft
      polygon( x = c(yrs,rev(yrs)),
               y = c(vbErr_qft[1,f,1:nYrs],vbErr_qft[3,f,nYrs:1]),
                col = scales::alpha(blobCols[k],.5), border = NA)
      lines( x = yrs, y = vbErr_qft[2,f,1:nYrs], col = blobCols[k],
              lty = k, lwd = 3 ) 
        
    }

    # if(f == max(plotFleets) & k == nBlobs)
    #   legend( x = "bottomleft", bty = "n",
    #           legend = c( "1.0x Len SD",
    #                       "0.5x Len SD",
    #                       "0.1x Len SD"),
    #         lwd = 3, col = blobCols )
  }

  mtext( side = 2, text = "Relative Error in AM estimates", line = 3, outer = TRUE)
  mtext( side = 1, text = "Year", line = 3, outer = TRUE)

}


# plotAMREs()
plotMultiAMREs <- function( blobList = list(blob1 = blob),
                            sIdx = 1, pIdx = 1,
                            pt = 1, plotFleets = c(4,5),
                            fleetNames = c("Surface","Dive") )
{
  obj <- blobList[[1]]

  # Pull model dims
  nReps <- obj$nSims
  nS    <- obj$om$nS
  nP    <- obj$om$nP
  nX    <- obj$om$nX
  nF    <- obj$om$nF

  repObj <- obj$ctlList$opMod$histRpt
  
  nPlotFleets <- length(plotFleets)
  nErr        <- 2+nX+2*nPlotFleets
  nBlobs      <- length(blobList)


  xJitter   <- seq( from = -.3, to = .3, length.out = nBlobs)
  blobCols  <- RColorBrewer::brewer.pal(nBlobs,"Dark2")

  browser()

  labels <- c("B0","totB0","R0","h","Mb","m1","M0")
  
  for(f in 1:length(plotFleets))
    labels <- c(labels,paste0("q_",fleetNames[f]),paste0("tau_",fleetNames[f]))


  plot( x = c(0.5,nErr+.5), xlab = "", ylab = "",
        y = c(-1,1), axes = FALSE, type = "n",
        xaxs = "i")
    axis( side = 1, at = 1:nErr, labels = labels )
    axis(side = 2, las =1)
    abline(v = 0.5 + 1:(nErr-1))
    abline(h = 0, lty = 2, col = "grey60")
    box()

    for( k in 1:nBlobs)
    {
      errList <- makeAMREs( blobList[[k]],
                            sIdx = sIdx, pIdx = pIdx,
                            pt = pt )$quantList
      # B0 error
      segments(x0 = xJitter[k] + 1, y0 = errList$errB0_q[1], y1 = errList$errB0_q[3], lwd = 3, col = blobCols[k] )
      points(x = xJitter[k] + 1, y = errList$errB0_q[2], pch = 16, col = blobCols[k] )

      # h error
      segments(x0 = xJitter[k] + 2, y0 = errList$errh_q[1], y1 = errList$errh_q[3], lwd = 3, col = blobCols[k] )
      points(x = xJitter[k] + 2, y = errList$errh_q[2], pch = 16, col = blobCols[k] )    

      # M_m error
      segments(x0 = xJitter[k] + 3, y0 = errList$errM_qx[1,1], y1 = errList$errM_qx[3,1], lwd = 3, col = blobCols[k] )
      points(x = xJitter[k] + 3, y = errList$errM_qx[2,1], pch = 16, col = blobCols[k] )    

      # M_f error
      segments(x0 = xJitter[k] + 4, y0 = errList$errM_qx[1,2], y1 = errList$errM_qx[3,2], lwd = 3, col = blobCols[k] )
      points(x = xJitter[k] + 4, y = errList$errM_qx[2,2], pch = 16, col = blobCols[k] )    

      # Catchability error
      for(j in 1:nPlotFleets)
      {
        fIdx <- plotFleets[j]
        segments(x0 = xJitter[k] + 5 + (j-1)*2, y0 = errList$errq_qf[1,fIdx], y1 = errList$errq_qf[3,fIdx], lwd = 3, col = blobCols[k] )
        points(x = xJitter[k] + 5 + (j-1)*2, y = errList$errq_qf[2,fIdx], pch = 16, col = blobCols[k] )    

        segments(x0 = xJitter[k] + 6 + (j-1)*2, y0 = errList$errtau_qf[1,fIdx], y1 = errList$errtau_qf[3,fIdx], lwd = 3, col = blobCols[k] )
        points(x = xJitter[k] + 6 + (j-1)*2, y = errList$errtau_qf[2,fIdx], pch = 16, col = blobCols[k] )    
      }

    }
    legend( x = "topleft", bg ="white",
            legend = c( "Sim Lengths",
                        "Sim Index",
                        "Sim All Data"),
            lwd = 3, col = blobCols,
            pch = 16, pt.cex = 2 )

    mtext( side = 2, text = "Relative Error in AM estimates", line = 3)
}

# plotAMREs()
plotMultiAMestimates <- function( blobMultiList = list(CC = list(blob1 = blob)),
                                  histPosts = blob$ctlList$opMod$posts,
                                  sIdx = 1, pIdx = 1,
                                  pt = 1, plotFleets = c(1,2,3),
                                  fleetNames = c("RedFB","Seine","Gillnet") )
{

  nLists <- length(blobMultiList)

  labels      <- c("B0","totB0","R0","h","m1","M0","q")
  fleetLabels <- c(paste0("selAlpha_",fleetNames))

  allLabels <- c(labels,fleetLabels)
  nErr    <- length(labels) + length(fleetLabels)

  par(mfcol = c(nErr,nLists), mar = c(.2,2,.2,2),
      oma = c(5,5,2,5))


  for(l in 1:nLists)
  {
    blobList  <- blobMultiList[[l]]
    obj       <- blobList[[1]]


    # Pull model dims
    nReps <- obj$nSims
    nS    <- obj$om$nS
    nP    <- obj$om$nP
    nX    <- obj$om$nX
    nF    <- obj$om$nF

    repObj <- obj$ctlList$opMod$histRpt
    
    nPlotFleets <- length(plotFleets)
    nBlobs      <- length(blobList)


    xJitter   <- seq( from = -.3, to = .3, length.out = nBlobs + 1)
    blobCols  <- RColorBrewer::brewer.pal(nBlobs+1,"Dark2")

    
    blobLabels <- c("Posterior")

    for( k in 1:nBlobs)
      blobLabels <- c(blobLabels,blobList[[k]]$ctlList$ctl$mpName)

    # Now we make a plot for each parameter

    for( j in 1:nErr)
    {
      
        
      thisLabel <- allLabels[j]
      if(thisLabel == "B0")
        post_i <- histPosts$B0_ip[,1]

      if(thisLabel == "totB0")
        post_i <- histPosts$totB0_ip[,1]

      if(thisLabel == "R0")
        post_i <- histPosts$R0_ip[,1]

      if(thisLabel == "h")
        post_i <- histPosts$h_ip[,1]

      if(thisLabel == "m1")
        post_i <- histPosts$m1_ip[,1]

      if(thisLabel == "M0")
        post_i <- histPosts$M0_ip[,1]

      if(thisLabel == "q")
        post_i <- histPosts$q_ipg[,1,4]

      if(thisLabel == "selAlpha_RedFB")
        post_i <- histPosts$selAlpha_ipg[,1,1]

      if(thisLabel == "selAlpha_Seine")
        post_i <- histPosts$selAlpha_ipg[,1,2]

      if(thisLabel == "selAlpha_Gillnet")
        post_i <- histPosts$selAlpha_ipg[,1,3]

      post_q <- quantile(post_i, probs = c(0.025, 0.5, 0.975))

      # if(thisLabel %in% fleetLabels)
      #   browser()

      plot( x = c(0.5,nBlobs + 1 +.5), xlab = "", ylab = "",
            y = c(0.95*post_q[1],1.05*post_q[3]), axes = FALSE, type = "n",
            xaxs = "i")



        mfg <- par("mfg")

        if(mfg[1] == 1)
          mtext(side = 3, text = names(blobMultiList)[l], font = 2, line = 0)

        if(mfg[1] == mfg[3])
          axis( side = 1, at = 1:(nBlobs+1), labels = blobLabels, las = 2 )
        axis( side = 2, las =1)
        if(mfg[2] == 1)
          mtext(side = 2, text = allLabels[j], line = 4)

        if(mfg[2] == mfg[4])
        {
          axis(side = 4, las = 1)
          mtext(side = 4, text = allLabels[j], line = 4)
        }

        abline(h = 0, lty = 2, col = "grey60")
        grid()
        box()
      
      abline(h = post_q[2], lwd = 2, lty = 2)

      # B0 error
      segments( x0 = 1, 
                y0 = post_q[1], 
                y1 = post_q[3], lwd = 3, 
                lty = 1,
                col = blobCols[1] )
      points(x = 1, y = post_q[2], pch = 16, 
              col = blobCols[1], cex = 2 )


      for( k in 1:nBlobs)
      {

        quantList <- makeAMEstDists(  blobList[[k]],
                                      sIdx = sIdx, pIdx = pIdx,
                                      pt = pt )$quantList

        if(thisLabel == "B0")
          quant_q <- quantList$B0_q

        if(thisLabel == "totB0")
          quant_q <- quantList$totB0_q

        if(thisLabel == "R0")
          quant_q <- quantList$R0_q

        if(thisLabel == "h")
          quant_q <- quantList$h_q

        if(thisLabel == "m1")
          quant_q <- quantList$m1_q

        if(thisLabel == "M0")
          quant_q <- quantList$M0_q

        if(thisLabel == "q")
          quant_q <- quantList$q_qf[,4]

        if(thisLabel == "selAlpha_RedFB")
          quant_q <- quantList$SelA_qft[,1,1]

        if(thisLabel == "selAlpha_Gillnet")
          quant_q <- quantList$SelA_qft[,3,1]

        if(thisLabel == "selAlpha_Seine")
          quant_q <- quantList$SelA_qft[,2,1]

        # B0 error
        segments( x0 = k + 1, 
                  y0 = quant_q[1], 
                  y1 = quant_q[3], lwd = 3, 
                  lty = 1,
                  col = blobCols[k+1] )
        points(x = k + 1, y = quant_q[2], pch = 16, 
                col = blobCols[k+1], cex = 2 )
      }
        # h or
      # if(j == 1)
      #   legend( x = "topleft", bg ="white", bty = "n",
      #           legend = c( "Posterior",
      #                       "Sim Lengths",
      #                       "Sim Index",
      #                       "Sim All Data"),
      #           lwd = 3, col = blobCols,
      #           pch = 16, pt.cex = 2 )
    }

    # mtext( side = 2, text = "Relative Error in AM estimates", line = 3)
  }
}

# plotTulipData()
# Simulation envelopes of index and catch data,
# with the observed data in the history plotted over
# the envelopes. Sometimes used to test for exceptional
# circumstances (i.e., model mis-specification error).
# inputs: 
#     obj = simulation object (blob)
# outputs: none
# Author: SDNJ, LFR
# Date: Apr 30, 2024
plotTulipData <- function(obj = blob)
{

  idxData <- read.csv("history/SOG_BlendedIdx.csv")
  catData <- read.csv("history/SOGcatchData_2023.csv")

  # Pull qs
  q_s     <- obj$ctlList$opMod$histRpt$qhat_pg[,4]

  # Make blended index
  idxData <- idxData |>
              mutate( I_t = Dive + Surface,
                      Dive = Dive/I_t,
                      Surface = Surface/I_t,
                      q_t = Dive + Surface * q_s )

  catData <- catData |>
              group_by(Year) |>
              summarise( Cat = sum(Value))

  # Pull series
  SB_it <- obj$om$SB_ispt[,1,1,]
  I_it  <- obj$mp$data$I_ispft[,1,1,5,]
  q_it  <- obj$om$q_ispft[,1,1,5,]
  C_it  <- obj$om$C_ispt[,1,1,]

  Iscale_it <- I_it/q_it

  # Now pull historical data
  nT  <- obj$om$nT
  tMP <- obj$om$tMP
  pT  <- nT - tMP + 1
  yrs <- seq(from = 1951, by = 1, length.out = nT )

  # Take quantiles
  SB_qt <- apply(X = SB_it, FUN = quantile, MARGIN = 2, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  I_qt  <- apply(X = Iscale_it, FUN = quantile, MARGIN = 2, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  C_qt  <- apply(X = C_it, FUN = quantile, MARGIN = 2, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)



  par(mfrow = c(2,1), mar = c(.25, 1, .25, 1), oma = c(3,5,2,2))


  plot(x = range(yrs), y = c(0, max(SB_qt, I_qt, na.rm = T)),
        axes = FALSE, type = "n")
    axis(side = 2, las = 1)
    box()
    grid()
    mtext(side = 2, text = "Spawning Biomass (kt)", line = 3)
    polygon(x = c(yrs[tMP:nT],yrs[nT:tMP]), y = c(I_qt[1,tMP:nT],rev(I_qt[3,tMP:nT])), 
            col = "grey75", border = NA )
    polygon(x = c(yrs,rev(yrs)), y = c(SB_qt[1,],rev(SB_qt[3,])), col = "steelblue",
              border = NA)
    
    lines(x = yrs, y = SB_qt[2,], lwd = 3 )
    points(x = idxData$Year, y = idxData$I_t/idxData$q_t, pch = 16,
            col = c(rep("grey30", tMP-1),rep("salmon",pT)) )
    legend( x = "topleft",
            legend = c( "Median SSB",
                        "SSB envelope",
                        "SI pre 2019",
                        "Observed SI 2019+",
                        "SI sim envelope"),
            pch = c(NA,15,16,16,15),
            col = c("black","steelblue","grey30","salmon","grey75"),
            lty = c(1,NA,NA,NA,NA),
            lwd = c(3,NA,NA,NA,NA),
            cex = 1.5, bty = "n")

  plot(x = range(yrs), y = c(0, max(C_it, na.rm = T)),
        axes = FALSE, type = "n")
    axis(side = 2, las =1 )
    box()
    grid()
    axis(side = 1)
    mtext(side = 2, text = "Catch (kt)", line = 3)
    mtext(side = 1, text = "Year", line = 2)
    polygon(x = c(yrs,rev(yrs)), y = c(C_qt[1,],rev(C_qt[3,])), col = "grey75", border = NA)
    lines(x = yrs, y = C_qt[2,], lwd = 3 )
    points(x = catData$Year, y=catData$Cat,
            pch = 16, col = c(rep("grey30", tMP-1),rep("salmon",pT)))

    legend( x = "topright",
            legend = c( "Median Catch",
                        "Catch sim envelope",
                        "Landings pre-2019",
                        "Landings 2019+"),
            pch = c(NA,15,16,16),
            col = c("black","grey75","grey30","salmon"),
            lty = c(1,NA,NA,NA),
            lwd = c(3,NA,NA,NA),
            cex = 1.5, bty = "n")


}

# Envelopes of simulated assessment errors
plotTulipAssError <- function(  obj = blob, 
                                pt = 1, 
                                sIdx = 1, 
                                pIdx = 1,
                                noPar = FALSE )
{
  
  goodReps <- obj$goodReps

  # Model dims
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT


  # Get biomass arrays
  SB_it        <- obj$om$SB_ispt[goodReps,sIdx,pIdx,1:(tMP-1)]
  retroSB_it   <- obj$mp$assess$retroSB_itspt[goodReps,pt,sIdx,pIdx,1:(tMP-1)]

  ctlList <- obj$ctlList

  retroSB_it[retroSB_it < 0] <- NA
  nReps     <- sum(goodReps)

  
  SB_it[SB_it == 0]     <- NA

  assErr_it <- array(NA, dim = c(nReps,tMP-1) )
  assErr_it[,1:(tMP-1)] <- (retroSB_it[,1:(tMP-1)] - SB_it[,1:(tMP-1)])/SB_it[,1:(tMP-1)]

  assErr_qt <- apply( X = assErr_it, na.rm = T,
                      FUN = quantile,
                      probs = c(0.025, 0.5, 0.975),
                      MARGIN = c(2) )

  years <- seq(from = 1951, by = 1, length.out = tMP-1 )

  tIdx <- 1:(tMP-1)

  traces <- sample( 1:nReps, 3)

  # Now plot
  plot( x = range(years[tIdx]),
        y = range(assErr_qt[,tIdx],na.rm = T ),
        type = "n", xlab = "", ylab = "", axes = FALSE )
        mfg <- par("mfg")
        axis( side = 1 )
        axis( side = 2, las = 1 )
        box()
        grid()

        # Plot baseline
        polygon(  x = c(years,rev(years)),
                  y = c(assErr_qt[1,],rev(assErr_qt[3,])),
                  border = NA, col = "grey60")
        lines( x = years, y = assErr_qt[2,],
                col = "black", lwd = 2)
        for( traceIdx in traces )
          lines(  x = years, y = assErr_it[traceIdx,],
                  col = "black", lwd = .8 )

        abline( h = 0, lty = 2, lwd = 1 )
  if(!noPar)
  {
    mtext(  side = 1, text = "Year", outer = TRUE, line = 2 )
    mtext(  side = 2, text = "Relative assessment error",
            line = 2, outer = TRUE )
  }

  
} # END plotTulipAssError

# Envelopes of simulated assessment errors
plotTulipAssErrorMort <- function(  obj = blob, 
                                    pt = 1, sIdx = 1, pIdx = 1 )
{
  
  goodReps <- obj$goodReps

  # Model dims
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT


  # Get biomass arrays
  M_it        <- obj$om$M_iaxspt[goodReps,2,1,sIdx,pIdx,1:(tMP-1)]
  retroM_it   <- obj$mp$assess$retroM_itaxspt[goodReps,pt,2,1,sIdx,pIdx,1:(tMP-1)]

  ctlList <- obj$ctlList

  retroM_it[retroM_it < 0] <- NA
  nReps     <- sum(goodReps)

  
  M_it[M_it == 0]     <- NA

  assErr_it <- array(NA, dim = c(nReps,tMP-1) )
  assErr_it[,1:(tMP-1)] <- (retroM_it[,1:(tMP-1)] - M_it[,1:(tMP-1)])/M_it[,1:(tMP-1)]

  assErr_qt <- apply( X = assErr_it, na.rm = T,
                      FUN = quantile,
                      probs = c(0.025, 0.5, 0.975),
                      MARGIN = c(2) )

  years <- seq(from = 1951, by = 1, length.out = tMP-1 )

  tIdx <- 1:(tMP-1)

  traces <- sample( 1:nReps, 3)

  # Now plot
  plot( x = range(years[tIdx]),
        y = range(assErr_qt[,tIdx],na.rm = T ),
        type = "n", xlab = "", ylab = "", axes = FALSE )
        mfg <- par("mfg")
        axis( side = 1 )
        axis( side = 2, las = 1 )
        box()
        grid()

        # Plot baseline
        polygon(  x = c(years,rev(years)),
                  y = c(assErr_qt[1,],rev(assErr_qt[3,])),
                  border = NA, col = "grey60")
        lines( x = years, y = assErr_qt[2,],
                col = "black", lwd = 2)
        for( traceIdx in traces )
          lines(  x = years, y = assErr_it[traceIdx,],
                  col = "black", lwd = .8 )

        abline( h = 0, lty = 2, lwd = 1 )
  mtext(  side = 1, text = "Year", outer = TRUE, line = 2 )
  mtext(  side = 2, text = "Relative assessment error",
          line = 2, outer = TRUE )

  
} # END plotTulipAssError



# plotTulipBtCtBaseSim()
# Overlays the biomass can catch tulips
# from the baseline omniscient sim (black) and 
# stochastic simulation (red) so we can see
# how the different results are made
plotTulipBtCtBaseSim <- function( sim = 1, 
                                  groupFolder = "diffCV_fixedF_longGrid",
                                  var         = "SB_ispt",
                                  save        = FALSE,
                                  fYear       = 1956,
                                  proj        = TRUE )
{
  # Load loss reports
  objLoss <- .loadLoss(sim = sim, groupFolder)

  # Load blob, save reference points
  .loadSim(sim = sim, groupFolder)

  rp <- blob$rp[[1]]

  scenName  <- blob$ctlList$ctl$scenarioName
  mpName    <- blob$ctlList$ctl$mpName

  stamp <- paste(scenName,":",mpName, sep = "")

  EmsyMSRefPts <- rp$EmsyMSRefPts
  FmsyRefPts   <- rp$FmsyRefPts


  simFolder <- objLoss$simFolder

  speciesNames    <- objLoss$speciesNames
  stockNames      <- objLoss$stockNames

  baseState_ispt  <- objLoss$baseStates[[var]]
  simState_ispt   <- objLoss$simStates[[var]]

  B0_isp          <- objLoss$baseStates[["SB_ispt"]][,,,1]


  # Pull model dimensions
  tMP <- objLoss$tMP
  nT  <- objLoss$nT
  nF  <- objLoss$nF
  nS  <- objLoss$nS
  nP  <- objLoss$nP
  pT  <- objLoss$pT

  if( proj )
    tIdx <- (tMP-1):nT
  else tIdx <- 1:nT

  years <- seq( from = fYear, by = 1, length.out = nT )

  nReps <- dim(simState_ispt)[1]
  
  traces <- sample(1:nReps, size = 3)

  # BeqFmsySS_sp <- array(NA, dim =c(nS+1, nP+1))
  # BeqEmsyMS_sp <- array(NA, dim =c(nS+1, nP+1))

  # BeqFmsySS_sp[1:nS,1:nP] <- FmsyRefPts$BeqFmsy_sp
  # BeqFmsySSS_sp[1:nS,p+1] <- apply( X = FmsyRefPts$BeqFmsy_sp, 
  #                                   FUN =  )
  
  if( save )
  {
    graphics.off()
    filename <- paste("baseSimTulipOverlay_",var,".pdf")
    outFile <- file.path(simFolder,filename)

    pdf( outFile, width = 11, height = 8 )
  }

  # put on depletion scale
  if( var == "SB_ispt")
  {
    maxY_sp <- matrix(1, nrow = nS + 1, ncol = nP +1 )
    yLab = expression(paste("Spawning biomass relative to unfished (", B[t]/B[0], ")",sep = ""))
    for( p in 1:(nP+1) )  
    {
      
      for( s in 1:(nS+1) )
        for( i in 1:nReps)
        {
          baseState_ispt[i,s,p,]  <- baseState_ispt[i,s,p,] / B0_isp[i,s,p]
          simState_ispt[i,s,p,]   <- simState_ispt[i,s,p,] / B0_isp[i,s,p]
        }

      maxY_sp[,p] <- max( baseState_ispt[,,p,tIdx], simState_ispt[,,p,tIdx])
    }

    
  } else {
    yLab    <- "Catch (kt)"
    maxY_sp <- matrix(NA, nrow = nS+1, ncol = nP+1 )

    for( s in 1:(nS+1))
      for( p in 1:(nP+1))
        maxY_sp[s,p] <- max( baseState_ispt[,s,p,tIdx], simState_ispt[,s,p,tIdx], na.rm = T )
  }

  baseState_qspt <- apply( X = baseState_ispt, FUN = quantile,
                            probs = c(0.025, 0.5, 0.975 ),
                            MARGIN = c(2,3,4) )

  simState_qspt <- apply( X = simState_ispt, FUN = quantile,
                            probs = c(0.025, 0.5, 0.975 ),
                            MARGIN = c(2,3,4) )

  par( mfcol = c((nP+1), nS+1),
        mar = c(0.1,0.1,0.1,0.1),
        oma = c(4,4.5,3,3) )

  if( var == "C_ispt" )
    par( mar = c(.1, 1.5, .1, 1.5) )


  basePolyCol <- "grey50"
  baseLineCol <- "black"
  simPolyCol <- scales::alpha("red",.3)
  simLineCol <- "red"

  for( s in 1:(nS+1) )
    for( p in 1:(nP+1) )
    {
      plot( x = range(years[tIdx]),
            y = c(0, maxY_sp[s,p] ),
            type = "n", xlab = "", ylab = "", axes = FALSE )
        mfg <- par("mfg")
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        if( mfg[2] == 1 )
          axis( side = 2, las = 1 )
        if( var == "C_ispt" & mfg[2] != 1 )
          axis( side = 2, las = 1)
        if( mfg[2] == mfg[4] )
          mtext( side = 4, text = stockNames[p], line = 0.5)
        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s])
        box()
        grid()

        # Plot baseline
        polygon(  x = c(years,rev(years)),
                  y = c(baseState_qspt[1,s,p,],rev(baseState_qspt[3,s,p,])),
                  border = NA, col = basePolyCol)
        lines( x = years, y = baseState_qspt[2,s,p,],
                col = baseLineCol, lwd = 2)
        for( traceIdx in traces )
          lines(  x = years, y = baseState_ispt[traceIdx,s,p,],
                  col = baseLineCol, lwd = .8 )

        # Plot stochastic sim
        polygon(  x = c(years,rev(years)),
                  y = c(simState_qspt[1,s,p,],rev(simState_qspt[3,s,p,])),
                  border = NA, col = simPolyCol)
        lines( x = years, y = simState_qspt[2,s,p,],
                col = simLineCol, lwd = 2)
        for( traceIdx in traces )
          lines(  x = years, y = simState_ispt[traceIdx,s,p,],
                  col = simLineCol, lwd = .8 )

        # Plot some lines
        abline( v = years[tMP], lty = 2, lwd = 1 )
        if( s <= nS & p <= nP & var == "SB_ispt")
        {
          abline( h = FmsyRefPts$BeqFmsy_sp[s,p]/B0_isp[1,s,p], 
                  lty = 3, lwd = 1 )
          abline( h = EmsyMSRefPts$BeqEmsy_sp[s,p]/B0_isp[1,s,p], 
                  lty = 4, lwd = 1 )
        }


    }
  mtext( side = 1, text = "Year", outer = TRUE, line = 2 )
  mtext( side = 2, text = yLab, outer = TRUE, line = 2 ) 

  mtext( side = 1, text = stamp, outer = TRUE, line = 2.5,
          adj = .85, cex = .6, col = "grey75")

  if(save)
    dev.off()

}



# plotBatchConvergenceRate()
# Abandoned performance metric for simulations,
# basically measuring the usefulness of a model
# under lower data quality conditions - can be
# used to set a threshold for including
# results in the paper...
plotBatchConvergenceRate <- function( groupFolder = "diffCV_newObsCV_short",
                                      prefix = "diffCV",
                                      AMlabs = c( singleStock = "singleStock",
                                                  hierMultiStock = "hierMultiStock",
                                                  dataPooled = "dataPooled",
                                                  coastwide = "coastwide",
                                                  totalAgg = "totalAgg" ) )
{
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs(here::here("Outputs",groupFolder))
  simFolderList <- simFolderList[grepl(prefix, simFolderList)]


  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
                filter( grepl( prefix, simLabel ))

  # Break up MP names into factor levels
  # MP labels are AM_Fsrce_eqbm
  splitMP <- function(  mpLab, 
                        breaks = "_",
                        factorNames = c("AM","Fsrce","eqbm") )
  {
    splitMP <- stringr::str_split(mpLab, breaks)[[1]]

    outList <- vector(  mode = "list", 
                        length = length(splitMP) )
    names(outList) <- factorNames
    for( k in 1:length(splitMP))
      outList[[k]] <- splitMP[k]

    outList
  }

  MPfactors <- lapply( X = info.df$mp, FUN = splitMP )
  MPfactors <- do.call(rbind, MPfactors)

  # Read in the performance tables
  mpTable <- cbind( info.df, MPfactors ) %>%
              mutate( perfPath = here::here("Outputs",groupFolder,simLabel,"simPerfStats.csv") )

  
  perfTableList <- lapply(  X = mpTable$perfPath, 
                            FUN = read.csv,
                            header = TRUE,
                            stringsAsFactors = FALSE )

  names( perfTableList ) <- mpTable$simLabel

  # Summarise the perf tables into convergence rates
  summariseTable <- function( table )
  {
    table <- table %>%
              group_by(simLabel) %>%
              summarise(  obsCVmult = mean(projObsErrMult),
                          minConv   = min(pGoodReps),
                          maxConv   = max(pGoodReps),
                          meanConv  = mean(pGoodReps) )

    table
  }

  # Summarise and make table
  perfTableSummList <- lapply(  X = perfTableList, 
                                FUN = summariseTable )

  convRateTable <- do.call(rbind, perfTableSummList )


  mpTable <- mpTable %>% left_join( convRateTable, 
                                    by = "simLabel" )

  CVmults <- unique(mpTable$obsCVmult)
  CVmults <- CVmults[order(CVmults)]

  AMs       <- names(AMlabs)
  AMcols    <- RColorBrewer::brewer.pal(length(AMs), "Set2")
  AMjitter  <- seq( from = -.3, to = .3, length.out = length(AMs) )

  AMwidth <- AMjitter[2] - AMjitter[1]

  # Now set up plotting area
  plot( x = c(0,length(CVmults) + 1), y = c(0,1),
        xlab = "", ylab = "", axes=FALSE,
        type = "n" )
    axis( side = 1, at = 1:length(CVmults),
          labels = CVmults )
    axis( side = 2, las = 1  )
    box()
    for( aIdx in 1:length(AMs) )
    {
      amLab <- AMs[aIdx]
      SStable <- mpTable %>% filter( AM == amLab,
                                      grepl("SS", eqbm) ) %>%
                              dplyr::select( obsCVmult,
                                              meanConv )

      SStable <- SStable[order(SStable$obsCVmult),]

      MStable <- mpTable %>% filter( AM == amLab,
                                      grepl("MS", eqbm) ) %>%
                              dplyr::select( obsCVmult,
                                              meanConv )
      MStable <- MStable[order(MStable$obsCVmult),]

      rect( xleft = 1:3 + AMjitter[aIdx] - AMwidth/2,
            xright = 1:3 + AMjitter[aIdx],
            ybottom = 0,
            ytop = SStable$meanConv,
            border = "black",
            col = AMcols[aIdx] )

      rect( xleft = 1:3 + AMjitter[aIdx],
            xright = 1:3 + AMjitter[aIdx] + AMwidth/2,
            ybottom = 0,
            ytop = MStable$meanConv,
            border = "black",
            col = AMcols[aIdx] )
      

      
    }



} # END plotBatchConvergenceRate()


# plotBatchLossDists_CV()
# Multi-panel plot of relative and absolute
# loss under a batch of MPs for all stock/species
# combinations, including data-pooled and 
# coast-wide aggregations.
plotBatchLossDists_CV <- function(  groupFolder = "diffCV_fixedF_longGrid",
                                    prefix = "sim_parBatfixedF",
                                    lossType = "abs",
                                    var = "C_ispt",
                                    period = 62:72,
                                    lossList = NULL,
                                    dim1 = 3,   # species (D,E,R, DP)
                                    dim2 = 1,   # stocks (H,Q,W, CW)
                                    qProbs = c(.05,.5,.95),
                                    refPts = "MSrefPts" )
{
  # First, read info files from the relevant
  # sims
  
  simFolderList <- list.dirs(here::here("Outputs",groupFolder),recursive = FALSE)
  simFolderList <- simFolderList[grepl(prefix, simFolderList)]

  infoFiles     <- lapply(  X = file.path(simFolderList,"infoFile.txt"),
                            FUN = lisread )


  if( lossType == "rel" )
  {
    lossArrayName <- "totRelLoss_isp"
    yLab <- "Total Relative Loss"
  }

  if( lossType == "abs" )
  {
    lossArrayName <- "totAbsLoss_isp"
    yLab <- "Total Absolute Loss (kt)"
  }

  # Break up MP names into factor levels
  # MP labels are AM_Fsrce_eqbm
  splitMP <- function(  mpLab, 
                        breaks = "_",
                        factorNames = c("AM","Fsrce","eqbm") )
  {
    splitMP <- stringr::str_split(mpLab, breaks)[[1]]

    outList <- vector(  mode = "list", 
                        length = length(splitMP) )
    names(outList) <- factorNames
    for( k in 1:length(splitMP))
      outList[[k]] <- splitMP[k]

    outList
  }



  for( lIdx in 1:length(infoFiles) )
  {
    infoFiles[[lIdx]] <- c(infoFiles[[lIdx]], splitMP(infoFiles[[lIdx]]$mp))
    infoFiles[[lIdx]] <- as.data.frame(infoFiles[[lIdx]])
  }

  info.df <- do.call(rbind, infoFiles) %>% 
              mutate_if(is.factor, as.character)

  if( !is.null(refPts) )
    info.df <- info.df %>% filter( eqbm == refPts )




  splitCVmult <- function( scenName )
  {
    splitUnderscore <- stringr::str_split( scenName, "_" )[[1]][2]

    splitObsErr     <- stringr::str_split( splitUnderscore, "obsErr")[[1]][1]

    CV <- as.numeric(splitObsErr)

    CV
  }

  info.df<- info.df %>%
              mutate( CVmult = sapply( X = scenario, FUN = splitCVmult) )


  # Load loss files, place in a list
  simFolderNames  <- info.df$simLabel
  if( is.null(lossList) )
  {
    lossList        <- lapply(  X = simFolderNames,
                                FUN = .loadLoss,
                                folder = groupFolder )
    names(lossList) <- simFolderNames
  }
  

  # Calculate total loss for variable/period
  totLossList  <- lapply( X = lossList,
                          FUN = calcTotalLossPeriod,
                          var = var, period = period )
  names(totLossList) <- names(lossList)


  nS  <- lossList[[1]]$nS
  nP  <- lossList[[1]]$nP
  nT  <- lossList[[1]]$nT
  tMP <- lossList[[1]]$tMP
  pT  <- lossList[[1]]$pT

  # Get number of MPs
  nSims <- length(lossList)


  AMs       <- unique(info.df$AM)
  Fsources  <- unique(info.df$Fsrce)
  eqbm      <- unique(info.df$eqbm)
  MPs       <- unique(info.df$mp)
  CVs       <- unique(info.df$CVmult)
  CVs       <- CVs[order(CVs)]

  nAM       <- length(AMs)
  nSrce     <- length(Fsources)
  nEqbm     <- length(eqbm)
  nMP       <- length(MPs)



  # Some axis labels and translations from
  # simple axis label to info file AM name
  AMcodes <- c( SS = 1,
                MS = 2,
                DP = 3,
                CW = 4,
                TA = 5 )

  AMdictionary <- c(  singleStock     = "SS",
                      hierMultiStock  = "MS",
                      dataPooled      = "DP",
                      coastwide       = "CW",
                      totalAgg        = "TA" )

  AMcols          <- RColorBrewer::brewer.pal(nAM, "Dark2")
  names(AMcols)   <- names(AMdictionary)
  eqbmPCH  <- 15 + 1:nEqbm
  names(eqbmPCH) <- eqbm
  FsrceLTY   <- 1:nSrce
  names(FsrceLTY)  <- Fsources


  MPgrid <- list( AMs = names(AMdictionary),
                  Fsrce = Fsources,
                  eqbm = eqbm )

  MPgrid <- expand.grid(MPgrid)


  xJitter <- seq(from = -.4, to = .4, length.out = nMP )


  nReps         <- dim(lossList[[1]]$lossRel[[var]])[1]
  speciesNames  <- lossList[[1]]$speciesNames
  stockNames    <- lossList[[1]]$stockNames

  # Make an array to hold loss function values
  totLossArray_misp <- array(NA,  dim = c(nSims, nReps, nS+1, nP+1 ),
                                  dimnames = list(  simID = info.df$simLabel,
                                                    rep = 1:nReps,
                                                    species = speciesNames,
                                                    stock = stockNames ) )

  for( k in 1:nSims )
  {
    simID <- info.df$simLabel[k]
    totLossArray_misp[k,,,] <- totLossList[[simID]][[lossArrayName]] 
  }

  totLossQuantiles_qmsp <- apply( X = totLossArray_misp,
                                  FUN = quantile,
                                  probs = qProbs,
                                  MARGIN = c(1,3,4) )

  # Start plotting windows
  nSpec   <- length(dim1)
  nStock  <- length(dim2) 

  par(  mfrow = c(nSpec,nStock),
        mar = c(0.2,1.5,.2,1),
        oma = c(3,3,3,3) )

  for( sIdx in 1:nSpec )
    for( pIdx in 1:nStock )
    {
      s <- dim1[sIdx]
      p <- dim2[pIdx]

      lossRange <- max(abs(range(totLossQuantiles_qmsp[,,s,p])))

      plot( y = c(0,lossRange), x = c(0.5,3.5),
            type = "n", axes = FALSE  )
        
        mfg <- par("mfg")
        if( mfg[1] == mfg[3])
          axis( side = 1, at = 1:length(CVs), labels = CVs )

        axis( side = 2, las = 1 )

        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2, line = 1 )

        if( mfg[2] == mfg[4] )
          mtext( side = 4, text = stockNames[p], font = 2, line = 1 )

        grid()
        box()
        
        jitIdx <- 0
        for( amIdx in 1:nAM )
          for( eqIdx in 1:nEqbm )
          {
            jitIdx <- jitIdx + 1
            AMid    <- names(AMdictionary)[amIdx]
            eqbmID  <- eqbm[eqIdx]
            
            subInfo <- info.df %>%
                        filter( AM    == AMid,
                                eqbm  == eqbmID ) %>%
                        arrange( CVmult )
            
            simLabels <- subInfo$simLabel
            
            FsrceID <- subInfo$Fsrce[1]
            AMcode  <- AMdictionary[amIdx]

            MPlab <- paste( AMid, FsrceID, eqbmID, sep = "_")

            points( x = 1:3 + xJitter[jitIdx],
                    y = totLossQuantiles_qmsp[2,simLabels,s,p],
                    pch = eqbmPCH[eqbmID],
                    col = AMcols[AMid], cex = 1.5 )
            # lines(  x = 1:3 + xJitter[jitIdx],
            #         y = totLossQuantiles_qmsp[2,simLabels,s,p],
            #         lty = 1, lwd = .8, col = AMcols[AMid] )

            segments( y0 = totLossQuantiles_qmsp[1,simLabels,s,p],
                      y1 = totLossQuantiles_qmsp[3,simLabels,s,p],
                      x0 = 1:3 + xJitter[jitIdx],
                      lty = FsrceLTY[FsrceID],
                      col = AMcols[AMid], lwd = 2  )

            abline( v = 1:2 + .5,
                    col = "black", lty = 2, lwd = 3 )


          }

      }


  legend( x       = "bottomright",
          bty     = "n",
          legend  = c(AMdictionary,eqbm),
          pch     = c(rep(NA,5),eqbmPCH),
          lty     = c(rep(1,5),NA,NA),
          col     = c(AMcols,"black","black"),
          lwd     = 2 )


  mtext( side = 2, outer = TRUE, text = yLab, line = 2 )
  mtext( side = 1, outer = TRUE, text = "Observation Error SD multiplier ", line = 2 )

  
}


# plotBatchLossDists()
# Multi-panel plot of relative and absolute
# loss under a batch of MPs for all stock/species
# combinations, including data-pooled and 
# coast-wide aggregations.
plotBatchLossDists <- function( groupFolder = "diffCV_fixedF_longGrid",
                                prefix = "sim_parBatfixedF",
                                lossType = "rel",
                                var = "C_ispt",
                                period = 62:72,
                                lossList = NULL )
{
  # First, read info files from the relevant
  # sims
  
  simFolderList <- list.dirs(here::here("Outputs",groupFolder))
  simFolderList <- simFolderList[grepl(prefix, simFolderList)]

  infoFiles     <- lapply(  X = file.path(simFolderList,"infoFile.txt"),
                            FUN = lisread )

  if( lossType == "rel" )
  {
    lossArrayName <- "totRelLoss_isp"
    yLab <- "Total Relative Loss"
  }

  if( lossType == "abs" )
  {
    lossArrayName <- "totAbsLoss_isp"
    yLab <- "Total Absolute Loss (kt)"
  }

  # Break up MP names into factor levels
  # MP labels are AM_Fsrce_eqbm
  splitMP <- function(  mpLab, 
                        breaks = "_",
                        factorNames = c("AM","Fsrce","eqbm") )
  {
    splitMP <- stringr::str_split(mpLab, breaks)[[1]]

    outList <- vector(  mode = "list", 
                        length = length(splitMP) )
    names(outList) <- factorNames
    for( k in 1:length(splitMP))
      outList[[k]] <- splitMP[k]

    outList
  }



  for( lIdx in 1:length(infoFiles) )
  {
    infoFiles[[lIdx]] <- c(infoFiles[[lIdx]], splitMP(infoFiles[[lIdx]]$mp))
    infoFiles[[lIdx]] <- as.data.frame(infoFiles[[lIdx]])
  }

  info.df <- do.call(rbind, infoFiles) %>% 
              mutate_if(is.factor, as.character)


  # Load loss files, place in a list
  simFolderNames  <- info.df$simLabel
  if( is.null(lossList) )
  {
    lossList        <- lapply(  X = simFolderNames,
                                FUN = .loadLoss,
                                folder = groupFolder )
    names(lossList) <- simFolderNames
  }
  

  # Calculate total loss for variable/period
  totLossList  <- lapply( X = lossList,
                          FUN = calcTotalLossPeriod,
                          var = var, period = period )
  names(totLossList) <- names(lossList)


  nS  <- lossList[[1]]$nS
  nP  <- lossList[[1]]$nP
  nT  <- lossList[[1]]$nT
  tMP <- lossList[[1]]$tMP
  pT  <- lossList[[1]]$pT

  # Get number of MPs
  nSims <- length(lossList)


  AMs       <- unique(info.df$AM)
  Fsources  <- unique(info.df$Fsrce)
  eqbm      <- unique(info.df$eqbm)

  nAM       <- length(AMs)
  nSrce     <- length(Fsources)
  nEqbm     <- length(eqbm)

  # Some axis labels and translations from
  # simple axis label to info file AM name
  AMcodes <- c( SS = 1,
                MS = 2,
                DP = 3,
                CW = 4,
                TA = 5 )

  AMdictionary <- c(  singleStock     = "SS",
                      hierMultiStock  = "MS",
                      dataPooled      = "DP",
                      coastwide       = "CW",
                      totalAgg        = "TA" )

  AMcols          <- RColorBrewer::brewer.pal(nAM, "Dark2")
  names(AMcols)   <- names(AMdictionary)
  FsrcePCH  <- 15 + 1:nSrce
  names(FsrcePCH) <- Fsources
  eqbmLTY   <- 1:nEqbm
  names(eqbmLTY)  <- eqbm

  nReps         <- dim(lossList[[1]]$lossRel[[var]])[1]
  speciesNames  <- lossList[[1]]$speciesNames
  stockNames    <- lossList[[1]]$stockNames

  # Make an array to hold loss function values
  totLossArray_misp <- array(NA,  dim = c(nSims, nReps, nS+1, nP+1 ),
                                  dimnames = list(  simID = info.df$simLabel,
                                                    rep = 1:nReps,
                                                    species = speciesNames,
                                                    stock = stockNames ) )

  for( k in 1:nSims )
  {
    simID <- info.df$simLabel[k]
    totLossArray_misp[k,,,] <- totLossList[[simID]][[lossArrayName]] 
  }

  totLossQuantiles_qmsp <- apply( X = totLossArray_misp,
                                  FUN = quantile,
                                  probs = c(0.25, 0.5, 0.75),
                                  MARGIN = c(1,3,4) )

  

  # Loop and plot - start by plotting ALL
  # combos
  yJitter <- seq(from = -.3, to = .3, length = nSrce * nEqbm )
  names(yJitter) <- apply(expand.grid(Fsources, eqbm), 1, paste, collapse="_")

  # Plot total loss for given period
  # on grid of Species/stocks
  par(  mfcol = c(nS+1,nP+1), 
        mar = c(.2,1.5,.2,1),
        oma = c(3,3.5,3,3) )

  for( s in 1:(nS+1) )
    for( p in 1:(nP+1) )
    {
      lossRange <- max(abs(range(totLossQuantiles_qmsp[,,s,p])))

      plot( y = c(0,lossRange), x = c(0,6),
            type = "n", axes = FALSE  )
        
        mfg <- par("mfg")
        if( mfg[1] == mfg[3])
          axis( side = 1, at = AMcodes, labels = names(AMcodes) )

        axis( side = 2, las = 1 )

        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2, line = 1 )

        if( mfg[2] == mfg[4] )
          mtext( side = 4, text = stockNames[p], font = 2, line = 1 )

        grid()
        box()

        for( k in 1:nSims )
        {
          simID   <- info.df$simLabel[k]
          AMid    <- info.df$AM[k]
          FsrceID <- info.df$Fsrce[k]
          eqbmID  <- info.df$eqbm[k]
          AMcode  <- AMdictionary[AMid]

          jitY    <- yJitter[paste(FsrceID,eqbmID,sep = "_")]

          points( y   = totLossQuantiles_qmsp[2,simID,s,p],
                  x   = AMcodes[AMcode] + jitY,
                  pch = FsrcePCH[FsrceID],
                  col = AMcols[AMid], cex = 1.5 )
          segments( y0 = totLossQuantiles_qmsp[1,simID,s,p],
                    y1 = totLossQuantiles_qmsp[3,simID,s,p],
                    x0 = AMcodes[AMcode] + jitY,
                    lty = eqbmLTY[eqbmID],
                    col = AMcols[AMid], lwd = 2  )
        }
    } # END p loop
    # END s loop

  legend( x       = "bottomright",
          bty     = "n",
          legend  = c(Fsources,eqbm),
          pch     = c(FsrcePCH,NA,NA),
          lty     = c(NA,NA,eqbmLTY),
          lwd     = 2 )


  mtext( side = 2, outer = TRUE, text = yLab, line = 2 )
  mtext( side = 1, outer = TRUE, text = "AM Configuration", line = 2 )


} # END plotBatchLossDists()

# plotTotLossDists()
# Plots relative and absolute loss
# for a given simulation. Requires
# loss for a given baseline to have
# been calculated first, and saved
# into the sim folder
plotTotLossDists <- function( sim = 1, 
                              groupFolder = "shortGrid",
                              lossType    = "rel",
                              var         = "SB_ispt",
                              save        = FALSE )
{
  # Load loss reports
  objLoss <- .loadLoss(sim = sim, groupFolder)
  simFolder <- objLoss$simFolder

  # Species/stock names and model dims
  speciesNames    <- objLoss$speciesNames
  stockNames      <- objLoss$stockNames
  tMP             <- objLoss$tMP
  nT              <- objLoss$nT 
  nF              <- objLoss$nF 
  nS              <- objLoss$nS 
  nP              <- objLoss$nP 
  pT              <- objLoss$pT
  simFolder       <- objLoss$simFolder
  
  loss_shortTerm  <- calcTotalLossPeriod(objLoss,var, period = tMP:(tMP+10))
  loss_medTerm    <- calcTotalLossPeriod(objLoss,var, period = tMP:(tMP+20))
  loss_longTerm   <- calcTotalLossPeriod(objLoss,var, period = tMP:nT)

  if( lossType == "rel" )
  {
    lossListName  <- "totRelLoss_isp"  
    yLab          <- "Total relative loss (unitless)"
  }

  if( lossType == "abs" )
  {
    lossListName  <- "totAbsLoss_isp"    
    yLab          <- "Total loss (kt)"
  }

  if( save )
    pdf( file = file.path(simFolder,paste(lossType,var,"LossCleveland.pdf",sep = "_")),
          width = 11, height = 8.5 )

  # Plot window - just for loss right now
  par(  mfcol = c(nS + 1, nP + 1),
        mar = c(.5,.5,.5,.5),
        oma = c(3,3,3,3) )

  for( s in 1:(nS+1) )
    for( p in 1:(nP + 1) )
    {
      lossShort_i <- loss_shortTerm[[lossListName]][,s,p]
      lossMed_i   <- loss_medTerm[[lossListName]][,s,p]
      lossLong_i  <- loss_longTerm[[lossListName]][,s,p]

      lossShort_q <- quantile(lossShort_i, probs = c(0.025, 0.5, 0.975) )
      lossMed_q   <- quantile(lossMed_i, probs = c(0.025, 0.5, 0.975) )
      lossLong_q  <- quantile(lossLong_i, probs = c(0.025, 0.5, 0.975) )

      maxLoss     <- max( abs(lossShort_q[is.finite(lossShort_q)]),
                          abs(lossMed_q[is.finite(lossMed_q)]),
                          abs(lossLong_q[is.finite(lossLong_q)]), na.rm = T )
      
      plot( x = c(0,4), y = c(0,maxLoss), type = "n",
            axes = FALSE )
        mfg <- par("mfg")
        if( mfg[1] == mfg[3] )
          axis( side = 1, at = 1:3, labels = c("S", "M", "L") )
        if( mfg[2] == 1 )
          axis( side = 2, las = 1)
        if( mfg[1] == 1 )
        {
          mtext( side = 3, text = speciesNames[s] )
        }
        if( mfg[2] == mfg[4] )
        {
          mtext( side = 4, text = stockNames[p] )
        }
        grid()
        box()

        segments( x0 = 1, y0 = lossShort_q[1], y1 = lossShort_q[3]  )
        points( x = 1, y = lossShort_q[2]  )
        segments( x0 = 2, y0 = lossMed_q[1], y1 = lossMed_q[3]  )
        points( x = 2, y = lossMed_q[2]  )
        segments( x0 = 3, y0 = lossLong_q[1], y1 = lossLong_q[3]  )
        points( x = 3, y = lossLong_q[2]  )

    }

  mtext( side = 1, text = "Time Period", outer = T, line = 2)
  mtext( side = 2, text = yLab, outer = T, line = 2)

  if( save )
    dev.off()
} # END plotLossDists()

# plotLossTulips()
# Plots relative and absolute loss
# for a given simulation. Requires
# loss for a given baseline to have
# been calculated first, and saved
# into the sim folder
plotLossTulip <- function(  sim = 1, 
                            groupFolder = "shortGrid",
                            lossType    = "rel",
                            var         = "SB_ispt",
                            save        = FALSE )
{
  # Load loss reports
  objLoss <- .loadLoss(sim = sim, groupFolder)
  simFolder <- objLoss$simFolder

  if(lossType == "rel" )
  {
    loss_ispt <- objLoss$lossRel[[var]]
    yLab      <- "Relative Loss"
  }

  if(lossType == "abs" )
  {
    loss_ispt <- objLoss$lossRaw[[var]]
    yLab      <- "Absolute loss (kt)"
  }

  # Pull model dimensions
  tMP <- objLoss$tMP
  nT  <- objLoss$nT
  nF  <- objLoss$nF
  nS  <- objLoss$nS
  nP  <- objLoss$nP
  pT  <- objLoss$pT


  nReps <- dim(loss_ispt)[1]
  
  traces <- sample(1:nReps, size = 3)

  # Time steps
  fYear   <- objLoss$fYear
  tdxPlot <- tMP:nT
  years   <- seq(from = fYear, by = 1, length.out = nT)

  if( save )
  { 
    graphics.off()
    pdf( file = file.path(simFolder,paste(lossType,var,"LossEnvelopes.pdf",sep = "_")),
          width = 11, height = 8.5 )
  }

  # Plot window - just for loss right now
  par(  mfcol = c(nS + 1, nP + 1),
        mar = c(.5,.5,.5,.5),
        oma = c(3,3,3,3) )

  speciesNames  <- c(objLoss$speciesNames,"data-pooled")
  stockNames    <- c(objLoss$stockNames,"coast-wide")

  loss_qspt <- apply( X = loss_ispt, FUN = quantile,
                      probs = c(0.025, 0.5, 0.975),
                      MARGIN = c(2,3,4) )

  plotYrs <- years[tdxPlot]


  # Loop and plot
  for( s in 1:(nS+1) )
    for( p in 1:(nP+1) )
    {
      maxLoss <- 1
      

      plot( x = range(years[tdxPlot]), y = c(-maxLoss,maxLoss),
            type = "n", axes = FALSE  )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[2] == 1 )
        axis( side = 2, las = 1)
      if( mfg[1] == 1 )
      {
        mtext( side = 3, text = speciesNames[s] )
      }
      if( mfg[2] == mfg[4] )
      {
        mtext( side = 4, text = stockNames[p] )
      }
      grid()
      box()

      abline( v = years[tMP], lty = 2, lwd = 3 )

      polygon(  x = c(years,rev(years)),
                y = c(loss_qspt[1,s,p,],rev(loss_qspt[3,s,p,])),
                border = NA, col = scales::alpha("grey10",.4) )
      lines(  x = years, y = loss_qspt[2,s,p,],
              col = "black", lwd = 3 )
      for( traceIdx in traces )
        lines( x = years, y = loss_ispt[traceIdx,s,p,],
                col = "black", lwd = 1 )
      grid()

      abline( h = 0, lty = 2, lwd = 1 )


    }
    mtext( side = 1, text = "Year", line = 2, outer = TRUE )
    mtext( side = 2, text = yLab, line = 2, outer = TRUE )

  if( save )
    dev.off()

} # END plotLoss()


# compareSimCondIdxData()
# Compares simulated index data to the data 
# from the conditioning assessment.
compareSimCondIdxData <- function(  obj = blob, sIdx = 1, 
                                    pIdx = 1, iRep = 1 )
{
  # get model dims and repObj
  nT <- obj$om$nT
  tMP <- obj$om$tMP
  nF <- obj$om$nF

  repObj <- obj$ctlList$opMod$histRpt

  # First, pull out the real data
  condI_ft              <- repObj$I_pgt[pIdx,,]
  condI_ft[5,1:(tMP-1)] <- repObj$combI_pt[pIdx,]
  simI_ft               <- obj$mp$data$I_ispft[iRep,sIdx,pIdx,,1:(tMP-1)]

  
  # Check which flets have data
  dataFleets <- c()
  for( f in 1:nF )
    if(any(condI_ft[f,] > 0))
      dataFleets <- c(dataFleets,f)

  condI_ft[condI_ft < 0] <- NA
  simI_ft[simI_ft < 0] <- NA

  fYear <- obj$ctlList$opMod$fYear
  yrs   <- seq( from = fYear, length.out = tMP + 3, by = 1) 

  # Multi-panel plot overlaying
  #   1. Data
  #   2. Residuals - annually with distro at RHS
  #   3. Ratios of designs - if blended

  # Get standardised resids (obs errors)
  condResid_ft <- array(NA, dim = c(nF,tMP-1))
  condResid_ft[dataFleets,] <- repObj$zComb_pt
  simResid_ft               <- obj$om$errors$delta_ispft[iRep,sIdx,pIdx,,1:(tMP-1)]

  # get ratios
  condRatioI_ft <- repObj$rI_pgt[pIdx,,]
  simRatioI_ft  <- obj$om$rI_ispft[iRep,sIdx,pIdx,,]

  par(  mfrow = c(3,length(dataFleets)),
        oma   = c(3,3,2,2),
        mar   = c(.5,2,.5,2) )
  for( f in dataFleets )
  {
    # First plot data
    plot(x = range(yrs), y = c(0,max(condI_ft[f,], simI_ft[f,],na.rm = T)),
          type = "n", axes = FALSE )
      axis( side = 2, las = 1)
      mtext(side = 2, text= "Spawn Index (kt)", line = 2)
      grid()
      box()
      points( x = yrs[1:(tMP-1)], y = condI_ft[f,],  pch = 16, col = "grey40" )
      points( x = yrs[1:(tMP-1)], y = simI_ft[f,],   pch = 21, col = "salmon", bg = NA )

    # Now plot residuals
    plot(x = range(yrs), y = range(condResid_ft[f,], simResid_ft[f,],na.rm = T),
          type = "n", axes = FALSE )
      abline(h = 0, lty = 2)
      mtext(side = 2, text= "Standardised residual", line = 2)
      axis( side = 2, las = 1)
      grid()
      box()
      points( x = yrs[1:(tMP-1)], y = condResid_ft[f,], pch = 16, col = "grey40" )
      segments( x0 = yrs[1], x1 = yrs[tMP-1],
                y0 = mean(condResid_ft[f,],na.rm = T), lty = 2, col = "grey40" )
      segments( x0 = yrs[1], x1 = yrs[tMP-1], lty = 2, col = "salmon",
                y0 = mean(simResid_ft[f,],na.rm = T) )
      points( x = yrs[1:(tMP-1)], y = simResid_ft[f,], pch = 21, col = "salmon", bg = NA )

    # need to add violin plots here

    # Now ratios - although they are the same
    plot(x = range(yrs), y = c(0,1),
          type = "n", axes = FALSE )
      axis( side = 2, las = 1)
      axis( side = 1 )
      mtext(side = 2, text= "Design ratio", line = 2)
      grid()
      box()
      rect( xleft = yrs[1:(tMP-1)]-0.3,
            xright =yrs[1:(tMP-1)], 
            ybottom = 0,
            ytop = condRatioI_ft[f,], col = "grey40" )
      rect( xleft = yrs[1:(tMP-1)],
            xright =yrs[1:(tMP-1)]+0.3, 
            ybottom = 0,
            ytop = simRatioI_ft[f,], col = "salmon" )

  }
  mtext(side = 1, text= "Year", outer = TRUE, line = 2)
} # END compareSimCondIdxData()

# compareSimCondAgeData()
# Compares simulated age data to the data used
# to fit the conditioning assessment.
compareSimCondAgeData <- function(  obj = blob, sIdx = 1, 
                                    pIdx = 1, iRep = 1,
                                    fIdx = 1 )
{
  # get model dims and repObj
  tMP <- obj$om$tMP
  nF  <- obj$om$nF
  nA  <- obj$om$nA

  repObj <- obj$ctlList$opMod$histRpt
  nT  <- repObj$nT
  yrs <- seq(from = 1951, length.out = nT)


  # Now, pull conditioned ages
  condA_at <- repObj$A_apgt[,pIdx,fIdx,]
  condA_at[condA_at < 0] <- NA
  sel_a <- repObj$sel_apgt[,pIdx,fIdx,1]

  for( a in 1:nA )
    if(sel_a[a] == 0)
      condA_at[a,] <- 0

  # Simulated ages
  simA_at  <- obj$mp$data$A_iaxspft[iRep,,1,sIdx,pIdx,fIdx,1:nT]

  # Now, we need to go through and sweep out the
  # sample sizes and turn into proportions, and
  # then plot
  condSampSize_t  <- apply(X = condA_at, FUN = sum, MARGIN = 2)
  simSampSize_t   <- apply(X = simA_at, FUN = sum, MARGIN = 2)

  for( t in 1:nT )
  {
    condA_at[,t]  <- condA_at[,t] / condSampSize_t[t]
    simA_at[,t]   <- simA_at[,t] / simSampSize_t[t]
  }

  # Check which years have data
  simDatIdx   <- which(!is.na(simSampSize_t))
  condDatIdx  <- which(!is.na(condSampSize_t))
  dataTimeIdx <- union(simDatIdx,condDatIdx)

  nPlotYrs  <- length(dataTimeIdx) + 2
  nCols     <- ceiling(sqrt(nPlotYrs))
  nRows     <- ceiling(sqrt(nPlotYrs))

  par(  mfcol = c(nRows,nCols), 
        oma = c(4,4,3,3),
        mar = c(.1,.1,.1,.1)
      )

  for( tIdx in dataTimeIdx )
  {
    yr <- yrs[tIdx]

    plot(x = c(1,nA), y = c(0,1), type = "n",
          axes = FALSE )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3] | tIdx == max(dataTimeIdx))
        axis(side = 1)
      if(mfg[1] == 1)
        axis(side = 3)
      if( mfg[2] == 1 )
        axis(side = 2, las = 1)
      if(mfg[2] == mfg[4])
        axis(side = 4, las = 1)
      grid()
      box()
      rect( xleft = 1:nA - .3,
            xright = 1:nA,
            ybottom = rep(0,nA),
            ytop = condA_at[,tIdx],
            col = "grey40", border = NA )
      rect( xleft = 1:nA,
            xright = 1:nA + .3,
            ybottom = rep(0,nA),
            ytop = simA_at[,tIdx],
            col = "salmon", border = NA )
      text( x = nA/2, y = 0.9,
            label = yr, font = 2 )
  }
  plot(x = c(1,nA), y = c(0,1), type = "n",
        axes = FALSE)
  legend( x = 2, y = 0.8, bty = "n",
          legend = c( "Real Data",
                      "Simulated"),
          pch = 15, cex = 1.5, col = c("grey40","salmon"))
  mtext( side = 1, text = "Age", outer= TRUE, line = 2)
  mtext( side = 2, text = "Proportion caught-at-age", outer= TRUE, line = 2)
}

plotSSvsMSrefPts <- function( obj = blob )
{
  # Pull reference points and curves
  rp            <- obj$rp[[1]]
  refCurves     <- rp$refCurves
  EmsyRefPts    <- rp$EmsyRefPts
  EmsyMSRefPts  <- rp$EmsyMSRefPts
  FmsyRefPts    <- rp$FmsyRefPts

  nP  <- obj$om$nP
  nS  <- obj$om$nS
  tMP <- obj$om$tMP
  nT  <- obj$om$nT

  # Pull qF
  qF_sp       <- blob$om$qF_ispft[1,,,2,nT]
  
  speciesNames <- blob$om$speciesNames
  stockNames   <- blob$om$stockNames

  # Want to compare Umsy and Bmsy between
  # multi-stock and single-stock models
  BmsySS_sp <- FmsyRefPts$BeqFmsy_sp
  BmsyMS_sp <- EmsyMSRefPts$BeqEmsy_sp

  UmsySS_sp <- FmsyRefPts$Umsy_sp
  UmsyMS_sp <- EmsyMSRefPts$Umsy_sp  

  B0_sp     <- obj$ctlList$opMod$histRpt$B0_sp

  speciesCols <- RColorBrewer::brewer.pal(n = nS, "Set2")

  cols_sp   <- matrix(  speciesCols, nrow = nS, ncol = nP,
                        byrow = FALSE )
  pch_sp    <- matrix(  21:23, nrow = nS, ncol = nP, 
                        byrow = TRUE )


  par(  mfrow = c(1,2),
        mar = c(1.5,4,0,0),
        oma = c(2,1,1,1) )

  # First plot Umsy
  plot( x = c(0,max(UmsySS_sp, UmsyMS_sp) ),
        y = c(0,max(UmsySS_sp, UmsyMS_sp) ),
        type = "n", 
        xlab = "", 
        ylab = "", las = 1 )
    mtext(side =1, text = "Single-species Umsy", line = 2 )
    mtext(side =2, text = "Multi-species Umsy", line = 3 )
    points( x = UmsySS_sp, y = UmsyMS_sp,
            cex = 3 * B0_sp / max(B0_sp),
            bg = cols_sp, pch = pch_sp )
    abline( a = 0, b = 1, lty = 2 )

  # First plot Umsy
  plot( x = c(0,max(BmsySS_sp, BmsyMS_sp) ),
        y = c(0,max(BmsySS_sp, BmsyMS_sp) ),
        type = "n", las = 1,
        xlab = "", 
        ylab = "" )
    mtext(side =1, text = "Single-species Bmsy", line = 2 )
    mtext(side =2, text = "Multi-species Bmsy", line = 2 )
    points( x = BmsySS_sp, y = BmsyMS_sp,
            cex = 3 * B0_sp / max(B0_sp),
            bg = cols_sp, pch = pch_sp )
    abline( a = 0, b = 1, lty = 2 )


} # END plotSSvsMSrefPts


# plotEffYieldCurves()
# Function for plotting effort based
# yield curves for each species and
# the complex in a stock area.
plotEffYieldCurves <- function( obj = blob )
{
  # First, pull reference points and curves
  rp            <- obj$rp[[1]]
  refCurves     <- rp$refCurves
  EmsyRefPts    <- rp$EmsyRefPts
  EmsyMSRefPts  <- rp$EmsyMSRefPts
  FmsyRefPts    <- rp$FmsyRefPts

  nT  <- obj$om$nT
  nP  <- obj$om$nP
  nS  <- obj$om$nS
  tMP <- obj$om$tMP

  # Pull qF
  qF_sp       <- blob$om$qF_ispft[1,,,2,nT]
  
  speciesNames <- blob$om$speciesNames

  # Now compute MSY for SS and MS curves
  # MSY is still the same for SS, just need the effort
  # that corresponds to it, which is Fmsy/qF

  # get SS ref points
  EmsySS_sp <- rp$FmsyRefPts$Fmsy_sp / qF_sp
  MSYSS_sp  <- rp$FmsyRefPts$YeqFmsy_sp

  # Now get MS ref points
  Yeq_spe   <- rp$refCurves$EffCurves$Yeq_spe
  Yeq_spe[Yeq_spe < 0] <- NA
  Yeq_pe    <- apply(X = Yeq_spe, FUN = sum, MARGIN = c(2,3), na.rm = T )
  Yeq_pe[Yeq_pe == 0] <- NA
  Yeq_pe[,1] <- 0
  Eseq      <- rp$refCurves$EffCurves$E

  EmsyMS_p  <- EmsyMSRefPts$EmsyMS_p
  MSYMS_sp  <- EmsyMSRefPts$YeqEmsy_sp
  MSYMS_p   <- EmsyMSRefPts$YeqEmsy_p
  BmsyMS_sp <- EmsyMSRefPts$BeqEmsy_sp

  Yeq_e <- apply( X = Yeq_pe, FUN = sum, na.rm = T,
                  MARGIN = 2 )
  maxEval <- max(which(Yeq_e > 0))

  specCols <- RColorBrewer::brewer.pal(n = nS, "Dark2")

  par( mfrow = c(3,1), mar = c(1,1,1,1), oma = c(3,3,1,1) )
  for( p in 1:nP )
  {

    plot( x = c(0,Eseq[maxEval]), y = c(0, max(Yeq_pe[p,],na.rm = T) ),
          type = "n", xlab = "", ylab = "", axes = F )
      axis(side = 2, las = 1)
      box()
      grid()

      for( s in 1:nS )
      {
        lines( x = Eseq, y = Yeq_spe[s,p,],
               col = specCols[s], lty = 1, lwd = 2 )
        # lines( x = Eseq, y = Beq_spe[s,p,],
        #        col = specCols[s], lty = 2, lwd = 2 )

      }
      mfg <- par("mfg")
      if( mfg[1] == mfg[3])
        axis(side = 1)

      lines(  x = Eseq, y = Yeq_pe[p,],
              col = "black", lty = 1, lwd = 3 )

      segments( x0 = EmsyMS_p[p], col = "grey40",
                y0 = 0, y1 = MSYMS_p[p], lty = 2  )

      segments( x0 = 0, x1 = EmsyMS_p[p], col = "grey40",
                y0 = c(MSYMS_p[p],MSYMS_sp[,p]), lty = 2  )

      if(  p == 1 )
        legend( x = "topright", bty = "n",
                col = c(specCols,"black"),
                lwd = c(2,2,2,3),
                legend = c(speciesNames,"Complex") )

  }

  mtext( outer = TRUE, side = 1, text = "Total Effort Index", line = 2 )
  mtext( outer = TRUE, side = 2, text = "Equilibrium Yield (kt)", line = 2 )


  outList <- list(  EmsyMS_p = EmsyMS_p,
                    EmsySS_sp = EmsySS_sp,
                    MSYMS_p = MSYMS_p,
                    MSYSS_sp = MSYSS_sp,
                    MSYMS_sp = MSYMS_sp )
} # END plotEffYieldCurves()

# plotFYieldCurves()
# Function for plotting effort based
# yield curves for each species and
# the complex in a stock area.
plotFYieldCurves <- function( obj = blob,
                              maxF = 1,
                              xLim =c(0,1.2),
                              yLim = c(0,50) )
{
  # First, pull reference points and curves
  rp            <- obj$rp[[1]]
  refCurves     <- rp$refCurves
  FmsyRefPts    <- rp$FmsyRefPts

  nT  <- obj$om$nT
  nP  <- obj$om$nP
  nS  <- obj$om$nS
  tMP <- obj$om$tMP

  
  speciesNames <- blob$om$speciesNames

  # Now compute MSY for SS and MS curves
  # MSY is still the same for SS, just need the effort
  # that corresponds to it, which is Fmsy/qF

  # get SS ref points
  Fmsy_sp   <- rp$FmsyRefPts$Fmsy_sp
  MSY_sp    <- rp$FmsyRefPts$YeqFmsy_sp
  Bmsy_sp   <- rp$FmsyRefPts$BeqFmsy_sp


  # Now get MS ref points
  Yeq_spf   <- rp$refCurves$Yeq_spf
  Yeq_spf[Yeq_spf < 0] <- NA
  Beq_spf   <- rp$refCurves$Beq_spf
  Beq_spf[Yeq_spf < 0] <- NA
  Beq_spf[Beq_spf < 0] <- NA
  Fseq      <- rp$refCurves$F

  maxF <- Fseq[max(which(!is.na(Yeq_spf[1,1,])))]

  specCols <- RColorBrewer::brewer.pal(n = nS, "Dark2")

  par( mfrow = c(nP,nS), mar = c(1,1,1,1), oma = c(3,3,1,1) )
  for( s in 1:nS )
    for( p in 1:nP )
    {
      if(is.null(xLim))
        xLim <- c(0,maxF)
      if(is.null(yLim))
        yLim <- c(0, max(Yeq_spf[s,p,],Beq_spf[s,p,],na.rm = T) )

      plot( x = xLim, y = yLim,
            type = "n", xlab = "", ylab = "", axes = F,
            xaxs = "i" )
        axis(side = 2, las = 1)
        box()
        grid()

        
        mfg <- par("mfg")
        if( mfg[1] == mfg[3])
          axis(side = 1)

        lines(  x = Fseq, y = Yeq_spf[s,p,],
                col = "black", lty = 1, lwd = 3 )
        lines(  x = Fseq, y = Beq_spf[s,p,],
                col = "steelblue", lty = 1, lwd = 3 )


        segments( x0 = Fmsy_sp[s,p], col = "red",
                  y0 = 0, y1 = Bmsy_sp[s,p], lty = 2  )

        segments( x0 = 0, x1 = Fmsy_sp[p], col = "red",
                  y0 = c(MSY_sp[s,p],Bmsy_sp[s,p]), lty = 2  )


        legend( x = "topright",
                legend = c( "Yield",
                            paste0("MSY = ",round(MSY_sp[s,p],2)),
                            "Spawning Biomass",
                            paste0("Bmsy = ",round(Bmsy_sp[s,p],2)),
                            paste0("Fmsy = ",round(Fmsy_sp[s,p],2)) ),
                bty = "n", lwd = c(2,NA,2,NA,NA),
                col = c("black",NA,"steelblue",NA,NA) )

    }

  mtext( outer = TRUE, side = 1, text = "Fishing mortality (/yr)", line = 2 )
  mtext( outer = TRUE, side = 2, text = "Equilibrium Yield/Biomass (kt)", line = 2 )


  outList <- list(  Fmsy_sp = Fmsy_sp,
                    MSY_sp  = MSY_sp  )
} # END plotEffYieldCurves()


# Plot predator catchability vs prey biomass
plotPredqvB <- function(  obj = blob, showProj = FALSE, iRep = 1 )
{
  # Get pred gears
  predGears <- obj$ctlList$opMod$predGears
  predType <- obj$ctlList$opMod$predType

  # just the functional response ones
  predGears <- predGears[predType != 0]

  # Now pull pred effort
  tMP <- obj$om$tMP
  nP  <- obj$om$nP
  nF  <- obj$om$nF
  nT  <- obj$om$nT

  predE_ft <- obj$om$E_ipft[1,1,,1:(tMP-1)]
  histF_ft <- obj$ctlList$opMod$histRpt$apF_pgt[1,,]

  # age 2+ biomass
  B_t   <- obj$om$B_ispt[1,1,1,1:(tMP-1)]
  vB_ft <- obj$ctlList$opMod$histRpt$vulnB_pgt[1,,1:(tMP-1)]

  q_ft <- histF_ft/predE_ft

  Brange <- c(0,max(B_t))

  gearNames <- dimnames(vB_ft)[[1]]

  qF_ift      <- obj$om$qF_ispft[,1,1,,tMP:nT]
  vBproj_ift  <- obj$om$vB_ispft[,1,1,,tMP:nT]


  eta_f     <- obj$om$eta_if[1,]
  h_f       <- obj$om$h_if[1,]
  lambda_f  <- obj$om$lambda_if[1,]

  

  nPred <- length(predGears)
  par(mfcol = c(nPred,1), oma = c(3,4.5,1,1),
        mar = c(.1,1,.1,1) )
  for( f in predGears)
  {
    qrange <- c(0,max(q_ft[f,],na.rm = T))

    qBtable <- data.frame(q = q_ft[f,], B = vB_ft[f,])

    # Smoother
    smoothline <- loess.smooth(x = vB_ft[f,], y = q_ft[f,])

    # Fit Disc
    Bseq <- seq( from = 0.001, to = max(vB_ft[f,]),
                  length.out = 200 )
    qseq <- eta_f[f] * Bseq^(lambda_f[f]-1) / (1 + h_f[f] * Bseq^lambda_f[f])

    qrange <- range(qseq,qrange)

    if(showProj)
      qrange <- range(qrange,qF_ift[,f,])

    plot(x = Brange, y = qrange, las = 1,
          type = "n", axes = FALSE )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis(side = 1)
      if(mfg[2] == 1)
        axis(side = 2)
      box()

      points( x = vB_ft[f,], y = q_ft[f,],
              col = "grey75", pch = 16, cex = 1.2 )

      lines(x = Bseq, y = qseq, lwd = 2, col = "darkgreen")
      lines(smoothline, lwd = 2, lty = 2, col = "black" )

      if(showProj)
        points( x = vBproj_ift[,f,],
                y = qF_ift[,f,] )

      
      legend( x = "topright",  bty ="n",
              legend = c( gearNames[f],
                          paste0("eta = ", signif(eta_f[f],2)),
                          paste0("h = ", signif(h_f[f],2)),
                          paste0("lambda = ", signif(lambda_f[f],2)),
                          "Loess Smoother",
                          "Disc equation" ),
              lty = c(NA,NA,NA,NA,2,1),
              lwd = c(NA,NA,NA,NA,2,2),
              col = c(NA,NA,NA,NA,"black","darkgreen") )

  }

  mtext(side = 1, outer = TRUE, text = "Vulnerable Herring Biomass (kt)", line = 2 )
  mtext( side = 2, outer = TRUE, text = "Predation Mortality per unit of predator", line= 3 )


}


solveSimEqbria <- function(empRefCurves,
                            maxXspline,
                            sIdx = 1, pIdx = 1)
{



  # Pull empirical ref curves
  C_k <- empRefCurves$C_spk[sIdx,pIdx,]
  B_k <- empRefCurves$B_spk[sIdx,pIdx,]
  F_k <- empRefCurves$F_spk[sIdx,pIdx,]
  U_k <- empRefCurves$U_spk[sIdx,pIdx,]
  M_k <- empRefCurves$M_spk[sIdx,pIdx,]
  R_k <- empRefCurves$R_spk[sIdx,pIdx,]
  
  # Order
  Forder <- order(F_k)
  C_k <- C_k[Forder]
  B_k <- B_k[Forder]
  F_k <- F_k[Forder]
  U_k <- U_k[Forder]
  M_k <- M_k[Forder]
  R_k <- R_k[Forder]

  # Make spline
  CUspline <- splinefun(x = U_k, y = C_k)
  BUspline <- splinefun(x = U_k, y = B_k)

  maxU <- max(U_k)
  maxC <- max(C_k)
  maxB <- max(B_k)
  maxM <- max(M_k)
  maxR <- max(R_k)

  Umsy  <- try(uniroot(  interval = c(0.01,maxXspline),
                        f = CUspline,
                        deriv = 1 )$root)
  if( class(Umsy) == "try-error")
    Umsy <- 0
  
  MSY  <- CUspline(Umsy)
  Bmsy <- BUspline(Umsy)

  

  outList <- list(  C_k = C_k,
                    U_k = U_k,
                    F_k = F_k,
                    B_k = B_k,
                    M_k = M_k,
                    R_k = R_k,
                    Umsy  = Umsy,
                    Bmsy  = Bmsy,
                    MSY   = MSY )

  outList
}

plotSimEqbria <- function(  folder = "may24_CC_refPtsProcError",
                            saveEmpRefCurves = NULL,
                            maxXspline = 0.1,
                            sIdx = 1, pIdx = 1 )
{
  # Load ref curve list
  if(is.null(saveEmpRefCurves))
    load(here::here("Outputs",folder,"empRefCurves.RData"))

  simEqList <- solveSimEqbria(  empRefCurves = saveEmpRefCurves,
                                maxXspline = maxXspline)


  # ref curves
  C_k <- simEqList$C_k
  B_k <- simEqList$B_k
  F_k <- simEqList$F_k
  U_k <- simEqList$U_k
  M_k <- simEqList$M_k
  R_k <- simEqList$R_k

  maxU <- max(U_k)
  maxC <- max(C_k)
  maxB <- max(B_k)
  maxR <- max(R_k)
  maxM <- max(M_k)

  # Ref pts
  Umsy  <- simEqList$Umsy
  Bmsy  <- simEqList$Bmsy
  MSY   <- simEqList$MSY

  empUmsy <- round(Umsy,2)
  empBmsy <- round(Bmsy,2)
  empMSY  <- round(MSY,2)

  maxUidx <- max(which(B_k > 1e-2))
  Ucrash  <- (U_k[maxUidx])
  Bcrash  <- (B_k[maxUidx])
  Ccrash  <- (C_k[maxUidx])
  Mcrash  <- (M_k[maxUidx])
  Rcrash  <- (R_k[maxUidx])


  maxU <- 1.2*Ucrash

  par(mfrow = c(4,1),
      mar = c(.1,2,.1,2),
      oma = c(4,4,1,1) )


  # First, plot yield
  plot( x = c(0,maxU), y = c(0,1.2*maxC),
        axes = FALSE, type = "n" )
    axis(side = 2, las =1)
    grid()
    box()
    lines( x = U_k[1:maxUidx], y = C_k[1:maxUidx], col = "black", lwd = 3)
    points( x = empUmsy, y = empMSY, pch = 16, col = "grey40", cex = 2)
    points( x = Ucrash, y = Ccrash, pch = 21, bg = "red", cex = 2)
    # points( x = U_k[1:maxUidx], y = C_k[1:maxUidx], col = "black")
    mtext(side = 2, text = "Yield (kt)", line = 3)
    legend( "topright", bty = "n", cex = 2,
            legend = c( paste("Umsy   = ", empUmsy, sep = "" ),
                        paste(" MSY   = ", empMSY, sep = "") ,
                        paste("Ucrash = ", round(Ucrash,2),sep = "") ) )

  # Biomass
  plot( x = c(0,maxU), y = c(0,maxB),
        axes = FALSE, type = "n" )
    axis(side = 2, las =1)
    grid()
    box()
    lines( x = U_k[1:maxUidx], y = B_k[1:maxUidx], col = "black", lwd = 3)  
    points( x = empUmsy, y = empBmsy, pch = 16, col = "grey40", cex = 2)
    points( x = Ucrash, y = Bcrash, pch = 21, bg = "red", cex = 2)
    legend( "topright", bty = "n", cex = 2,
            legend = c( paste("Bmsy = ", empBmsy, sep = "") ) )
    mtext(side = 2, text = "Spawning\nBiomass (kt)", line = 2.5)

  # Mortality
  plot( x = c(0,maxU), y = c(0,maxM),
        axes = FALSE, type = "n" )
    axis(side = 2, las =1)
    grid()
    box()
    lines( x = U_k[1:maxUidx], y = M_k[1:maxUidx], col = "black", lwd = 3)
    points( x = Ucrash, y = Mcrash, pch = 21, bg = "red", cex = 2)
    mtext(side = 2, text = "Natural\nMortality (/yr)", line = 2.5)

  # Recruitment
  plot( x = c(0,maxU), y = c(0,maxR),
        axes = FALSE, type = "n" )
    axis(side = 2, las =1)
    axis(side = 1)
    grid()
    box()
    lines( x = U_k[1:maxUidx], y = R_k[1:maxUidx], col = "black", lwd = 3)
    points( x = Ucrash, y = Rcrash, pch = 21, bg = "red", cex = 2)
    mtext(side = 2, text = "Recruitment (1e6)", line = 3)

    mtext(side = 1, text = "Harvest rate C / (C + B) (/yr)", line = 2.5)

}

# # plotEmpYieldCurves
# # Function to plot median simuated yield 
# # resulting from a grid of constant fishing 
# # mortality rates - used for checking the
# # reference point calculations
# plotEmpYieldCurves <- function( sims = 1:20, 
#                                 folder = "SOG_detRefPts",
#                                 indepVar = "U",
#                                 maxXspline = 0.5,
#                                 xLim = c(0,.5),
#                                 plotYPRcurves = FALSE,
#                                 redoEmpRefCurves = FALSE)
# {
#   nSims <- length(sims)
#   blobList <- vector(mode = "list", length = nSims)

#   .loadSim(sims[1], folder = folder)

#   nS <- blob$om$nS
#   nP <- blob$om$nP
#   nT <- blob$om$nT

#   goodReps <- blob$goodReps


#   # Arrays to hold empirical eqbm yields
#   C_spk <- array( 0, dim = c(nS,nP,nSims) )
#   B_spk <- array( 0, dim = c(nS,nP,nSims) )
#   F_spk <- array( 0, dim = c(nS,nP,nSims) )
#   U_spk <- array( 0, dim = c(nS,nP,nSims) )
#   E_spk <- array( 0, dim = c(nS,nP,nSims) )
#   M_spk <- array( 0, dim = c(nS,nP,nSims) )
#   R_spk <- array( 0, dim = c(nS,nP,nSims) )
#   E_pk  <- array( 0, dim = c(nP,nSims) )

#   qF_sp <- blob$om$qF_ispft[1,,,2,nT]

#   if(!file.exists(here::here("Outputs",folder,"empRefCurves.RData")) |
#       redoEmpRefCurves )
#   {
#     for( x in sims )
#     {
#       .loadSim(x, folder = folder)
#       C_spk[,,x]  <- apply(X = blob$om$C_ispt[goodReps,,,nT-0:10,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
#       B_spk[,,x]  <- apply(X = blob$om$SB_ispt[goodReps,,,nT-0:10,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
#       F_spk[,,x]  <- apply(X = blob$om$F_ispft[1,,,,nT,drop = FALSE], FUN = sum, MARGIN = c(2,3))
#       U_spk[,,x]  <- C_spk[,,x]/(C_spk[,,x] + B_spk[,,x])
#       M_spk[,,x]  <- apply(X = blob$om$M_iaxspt[goodReps,2,1,,,nT-0:10,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
#       R_spk[,,x]  <- apply(X = blob$om$R_ispt[goodReps,,,nT-0:10,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
#       E_pk[,x]    <- apply(X = blob$om$E_ipft[goodReps,,2,nT-0:10,drop = FALSE], FUN = median, MARGIN = c(2), na.rm = T )
#       E_spk[,,x]  <- F_spk[,,x] / qF_sp

#       # Clean up
#       gc()
#     }

#     saveEmpRefCurves <- list( C_spk = C_spk,
#                               B_spk = B_spk,
#                               F_spk = F_spk,
#                               U_spk = U_spk,
#                               E_spk = E_spk,
#                               M_spk = M_spk,
#                               R_spk = R_spk,
#                               E_pk  = E_pk )


#     save( saveEmpRefCurves, file = here::here("Outputs",folder,"empRefCurves.RData"))
#   } else {
#     # Load ref curve list
#     load(here::here("Outputs",folder,"empRefCurves.RData"))

#     C_spk <- saveEmpRefCurves$C_spk
#     B_spk <- saveEmpRefCurves$B_spk
#     F_spk <- saveEmpRefCurves$F_spk
#     U_spk <- saveEmpRefCurves$U_spk
#     E_spk <- saveEmpRefCurves$E_spk
#     E_pk  <- saveEmpRefCurves$E_pk 
#   }


#   # Pull F based ref curves from RP object
#   refCurves       <- blob$rp[[1]]$refCurves
#   Fvec            <- refCurves$F
#   Uvec            <- refCurves$Ueq_spf
#   BeqRefCurve_spf <- refCurves$Beq_spf
#   YeqRefCurve_spf <- refCurves$Yeq_spf

#   # Pull effort based ref points
#   EmsyRefPts      <- blob$rp[[1]]$EmsyRefPts
#   Evec            <- refCurves$EffCurves$E

#   YeqRefCurve_spe <- refCurves$EffCurves$Yeq_spe
#   BeqRefCurve_spe <- refCurves$EffCurves$Beq_spe

#   Xmsy_sp <- array(0, dim = c(nS,nP))
#   Bmsy_sp <- array(0, dim = c(nS,nP))
#   MSY_sp <- array(0, dim = c(nS,nP))

#   par(  mfcol = c(nP,nS), 
#         mar = c(1,1.5,1,1.5),
#         oma = c(5,5,3,3) )

#   for( s in 1:nS )
#     for( p in 1:nP )
#     {
#       if( indepVar == "F" )
#         maxX <- max(F_spk[s,p,])
#       if( indepVar == "E" )
#         maxX <- max(E_pk[p,],E_spk[s,p,])

#       if( indepVar == "U" )
#         maxX <- max(U_spk[s,p,])

#       if(xLim[2] < maxX )
#         maxX <- xLim[2]

#       plot( x = c(0,maxX), y = c(0,max(B_spk[s,p,])),
#             type = "n", axes = FALSE,
#             xlab = "", ylab = "" )
#         axis( side = 1 )
#         axis( side = 2, las = 1 )
#         grid()
#         box()

#         F <- F_spk[s,p,]
#         U <- U_spk[s,p,]
#         C <- C_spk[s,p,]
#         B <- B_spk[s,p,]
#         E <- E_spk[s,p,]

#         actualOrder <- order(F)

#         C <- C[actualOrder]
#         B <- B[actualOrder]
#         U <- U[actualOrder]
#         F <- F[actualOrder]
#         E <- E[actualOrder]

#         if( indepVar == "F" )
#         {
#           X     <- F
#           Xvec  <- Fvec
#           Yeq   <- YeqRefCurve_spf[s,p,]
#           Beq   <- BeqRefCurve_spf[s,p,]

#           xLab  <- "Fishing Mortality"
#         }

#         if( indepVar == "U" )
#         {
#           X     <- U
#           Xvec  <- Uvec
#           Yeq   <- YeqRefCurve_spf[s,p,]
#           Beq   <- BeqRefCurve_spf[s,p,]

#           xLab  <- "Harvest rate C /(B+C) (/yr)"
#         }

#         if( indepVar == "E" )
#         {
          
#           X     <- E
#           Xvec  <- Evec
#           Yeq   <- YeqRefCurve_spe[s,p,]
#           Beq   <- BeqRefCurve_spe[s,p,]

#           xLab  <- "Fishing Effort"
#         }


#         ubX <- max(which(Yeq >= 0) )

#         CXspline <- splinefun(x = X, y = C)
#         BXspline <- splinefun(x = X, y = B)


#         empXmsy  <- try(uniroot(  interval = c(0.01,maxXspline),
#                               f = CXspline,
#                               deriv = 1 )$root)
#         if( class(empXmsy) == "try-error")
#           empXmsy <- 0
        
#         empMSY  <- CXspline(empXmsy)
#         empBmsy <- BXspline(empXmsy)

#         empXmsy <- round(empXmsy,2)
#         empBmsy <- round(empBmsy,2)
#         empMSY  <- round(empMSY,2)

        
        
#           lines( x = X, y = C,
#                   col = "steelblue", lwd = 2, lty = 1 )
#           lines( x = X, y = B,
#                   col = "black", lwd = 2, lty = 1 )        
#         if(plotYPRcurves)
#         {
#           lines( x = Xvec, y = Yeq,
#                   col = "salmon", lty = 2 )
#           lines( x = Xvec, y = Beq,
#                   col = "black", lty = 2 )
#         }


#         legend( "topright", bty = "n",
#                 legend = c( paste("Umsy = ", empXmsy, sep = "" ),
#                             paste(" MSY = ", empMSY, sep = ""),
#                             paste("Bmsy = ", empBmsy, sep = "") ) )


#         # Plot some guidelines
#         segments(  x0 = empXmsy, x1 = empXmsy,
#                     y0 = 0, y1 = empBmsy,
#                     lty = 2, col = "red" )
#         segments(  x0 = 0, x1 = empXmsy,
#                     y0 = empMSY, y1 = empMSY,
#                     lty = 2, col = "red" )
#         segments(  x0 = 0, x1 = empXmsy,
#                     y0 = empBmsy, y1 = empBmsy,
#                     lty = 2, col = "red" )
        
#         Xmsy_sp[s,p] <- empXmsy
#         Bmsy_sp[s,p] <- empBmsy
#         MSY_sp[s,p]  <- empMSY
#     }

#   mtext( side = 1, text = xLab, outer = T, line = 2)
#   mtext( side = 2, text = "Eqbm biomass and catch (kt)", outer = T, line = 2 )



#   out <- list(  Xmsy_sp = Xmsy_sp,
#                 Bmsy_sp = Bmsy_sp,
#                 MSY_sp  = MSY_sp )

#   out
# } # END plotEmpYieldCurves()


# plotEmpYieldCurves
# Function to plot median simuated yield 
# resulting from a grid of constant fishing 
# mortality rates - used for checking the
# reference point calculations
plotEmpYieldCurves <- function( sims = 1:50, 
                                folder = "SOG_OM1_refPts",
                                indepVar = "U",
                                maxXspline = 0.6,
                                plotYPRcurves = FALSE,
                                redoEmpRefCurves = FALSE)
{
  nSims <- length(sims)
  blobList <- vector(mode = "list", length = nSims)

  .loadSim(sims[1], folder = folder)

  nS  <- blob$om$nS
  nP  <- blob$om$nP
  nT  <- blob$om$nT
  nF  <- blob$om$nF
  tMP <- blob$om$tMP

  pT <- nT - tMP + 1

  goodReps <- which(blob$goodReps)
  nReps <- length(goodReps)


  # Arrays to hold empirical eqbm yields
  C_isptk      <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )   # Removals (incl dP)
  Y_isptk      <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )   # Yield ( C + K )
  uC_isptk     <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )   # Catch (incl. tptal ponded fish)
  B_isptk      <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )
  eB_isptk     <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )
  tB_isptk     <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )
  F_isptk      <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )
  U_isptk      <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )
  M_isptk      <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )
  R_isptk      <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )
  
  
  
  qF_sp <- blob$om$qF_ispft[1,,,2,nT]

  if(!file.exists(here::here("Outputs",folder,"empRefCurves.RData")) |
      redoEmpRefCurves )
  {
    for( x in sims )
    {
      .loadSim(x, folder = folder)
      for( s in 1:nS)
        for( p in 1:nP)
        {
          # Yield is catch - does this include SOK or ponded fish?
          Y_isptk[1:nReps,s,p,,x]      <- blob$om$C_ispt[goodReps,s,p,tMP:nT]


          if( nF == 5 )
          {
            # Removals, for SP calcs
            C_isptk[1:nReps,s,p,,x]      <- blob$om$C_ispt[goodReps,s,p,tMP:nT]
            uC_isptk[1:nReps,s,p,,x]     <- blob$om$C_ispt[goodReps,s,p,tMP:nT]
          }
            
          if( nF > 5 )
          {
            # Add dead ponded fish to removals for production
            # calculations, and total ponded fish for HR calcs
            uC_ispft  <- C_ispft <- blob$om$C_ispft

            dP_ispft  <- blob$om$dP_ispft
            P_ispft   <- blob$om$P_ispft
            
            C_ispft[,,,6,]  <- dP_ispft[,,,6,]
            uC_ispft[,,,6,] <- P_ispft[,,,6,]

            C_isptk[1:nReps,s,p,,x] <- apply(X = C_ispft[goodReps,s,p,,tMP:nT,drop = FALSE], FUN = sum, MARGIN = c(1,2,3,5))
            uC_isptk[1:nReps,s,p,,x] <- apply(X = uC_ispft[goodReps,s,p,,tMP:nT,drop = FALSE], FUN = sum, MARGIN = c(1,2,3,5))
          }

          # use end of year SSB so ponded fish are added back in
          B_isptk[1:nReps,s,p,,x]      <- blob$om$SB_ispt[goodReps,s,p,tMP:nT]
          eB_isptk[1:nReps,s,p,,x]     <- blob$om$endSB_ispt[goodReps,s,p,tMP:nT]
          tB_isptk[1:nReps,s,p,,x]     <- blob$om$B_ispt[goodReps,s,p,tMP:nT]
          F_isptk[1:nReps,s,p,,x]      <- apply(X = blob$om$F_ispft[goodReps,s,p,,tMP:nT,drop = FALSE], FUN = sum, MARGIN = c(1,5))
          U_isptk[1:nReps,s,p,,x]      <- uC_isptk[,,,,x]/(uC_isptk[,,,,x] + B_isptk[,,,,x])
          M_isptk[1:nReps,s,p,,x]      <- blob$om$M_iaxspt[goodReps,2,1,s,p,tMP:nT]
          R_isptk[1:nReps,s,p,,x]      <- blob$om$R_ispt[goodReps,s,p,tMP:nT]
          
        }
      
      # Clean up
      gc()
    }

    B_sptk  <- apply(X = B_isptk, FUN = median, MARGIN = 2:5 )
    C_sptk  <- apply(X = C_isptk, FUN = median, MARGIN = 2:5 )
    uC_sptk <- apply(X = uC_isptk, FUN = median, MARGIN = 2:5 )
    Y_sptk  <- apply(X = Y_isptk, FUN = median, MARGIN = 2:5 )
    tB_sptk <- apply(X = tB_isptk, FUN = median, MARGIN = 2:5 )
    eB_sptk <- apply(X = eB_isptk, FUN = median, MARGIN = 2:5 )
    U_sptk  <- apply(X = U_isptk, FUN = median, MARGIN = 2:5 )
    F_sptk  <- apply(X = F_isptk, FUN = median, MARGIN = 2:5 )
    M_sptk  <- apply(X = M_isptk, FUN = median, MARGIN = 2:5 )
    R_sptk  <- apply(X = R_isptk, FUN = median, MARGIN = 2:5 )

    

    C_spk   <- apply(X = C_sptk[,,pT-(150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    uC_spk  <- apply(X = uC_sptk[,,pT-(150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    Y_spk   <- apply(X = Y_sptk[,,pT-(150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    B_spk   <- apply(X = B_sptk[,,pT-(150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    tB_spk  <- apply(X = tB_sptk[,,pT-(150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    eB_spk  <- apply(X = eB_sptk[,,pT-(150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    U_spk   <- apply(X = U_sptk[,,pT-(150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    F_spk   <- apply(X = F_sptk[,,pT - (150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    M_spk   <- apply(X = M_sptk[,,pT - (150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    R_spk   <- apply(X = R_sptk[,,pT-(150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)

    saveEmpRefCurves <- list( B_isptk     = B_isptk,
                              C_isptk     = C_isptk,
                              uC_isptk    = uC_isptk,
                              tB_isptk    = tB_isptk,
                              U_isptk     = U_isptk,
                              F_isptk     = F_isptk,
                              Y_isptk     = Y_isptk,
                              M_isptk     = M_isptk,
                              R_isptk     = R_isptk,
                              C_sptk      = C_sptk,
                              uC_sptk     = uC_sptk,
                              B_sptk      = B_sptk,
                              tB_sptk     = tB_sptk,
                              eB_sptk     = eB_sptk,
                              F_sptk      = F_sptk,
                              U_sptk      = U_sptk,
                              M_sptk      = M_sptk,
                              R_sptk      = R_sptk,
                              Y_sptk      = Y_sptk,
                              C_spk       = C_spk,
                              uC_spk      = uC_spk,
                              B_spk       = B_spk,
                              tB_spk      = tB_spk,
                              eB_spk      = eB_spk,
                              U_spk       = U_spk,
                              F_spk       = F_spk,
                              R_spk       = R_spk,
                              M_spk       = M_spk )


    save( saveEmpRefCurves, file = here::here("Outputs",folder,"empRefCurves.RData"))
  } else {
    # Load ref curve list
    load(here::here("Outputs",folder,"empRefCurves.RData"))

    C_sptk  <- saveEmpRefCurves$C_sptk
    uC_sptk <- saveEmpRefCurves$uC_sptk
    Y_sptk  <- saveEmpRefCurves$Y_sptk
    B_sptk  <- saveEmpRefCurves$B_sptk
    tB_sptk <- saveEmpRefCurves$tB_sptk
    F_sptk  <- saveEmpRefCurves$F_sptk
    U_sptk  <- saveEmpRefCurves$U_sptk
    M_sptk  <- saveEmpRefCurves$M_sptk
    R_sptk  <- saveEmpRefCurves$R_sptk
  
  }
  C_spk   <- apply(X = C_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  Y_spk   <- apply(X = Y_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  B_spk   <- apply(X = B_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  tB_spk  <- apply(X = tB_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  U_spk   <- apply(X = U_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  F_spk   <- apply(X = F_sptk[,,pT - (50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)

  # Pull F based ref curves from RP object
  refCurves       <- blob$rp[[1]]$refCurves
  Fvec            <- refCurves$F
  Uvec            <- refCurves$Ueq_spf
  BeqRefCurve_spf <- refCurves$Beq_spf
  YeqRefCurve_spf <- refCurves$Yeq_spf

  # browser()
  Dstar <- blob$ctlList$opMod$histCtl$hypo$Dstar - 1950
  Dscalar <- blob$ctlList$opMod$histCtl$hypo$Dscalar
  SB_ispt <- blob$om$SB_ispt

  B_usr <- Dscalar*mean(SB_ispt[,,,Dstar])

  Xmsy_sp <- array(0, dim = c(nS,nP))
  Bmsy_sp <- array(0, dim = c(nS,nP))
  MSY_sp <- array(0, dim = c(nS,nP))

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(5,5,3,3) )

  for( s in 1:nS )
    for( p in 1:nP )
    {
      if( indepVar == "F" )
        maxX <- max(F_spk[s,p,])
      if( indepVar == "E" )
        maxX <- max(E_pk[p,],E_spk[s,p,])

      if( indepVar == "U" )
        maxX <- max(U_spk[s,p,])

      plot( x = c(0,maxX), y = c(0,max(B_spk[s,p,])),
            type = "n", axes = FALSE,
            xlab = "", ylab = "" )
        axis( side = 1 )
        axis( side = 2, las = 1 )
        grid()
        box()

        F   <- F_spk[s,p,]
        U   <- U_spk[s,p,]
        Y   <- Y_spk[s,p,]
        C   <- C_spk[s,p,]
        uC  <- uC_spk[s,p,]
        B   <- B_spk[s,p,]
        # E <- E_spk[s,p,]

        actualOrder <- order(F)

        C   <- C[actualOrder]
        Y   <- Y[actualOrder]
        uC  <- uC[actualOrder]
        B   <- B[actualOrder]
        U   <- U[actualOrder]
        F   <- F[actualOrder]
      

        if( indepVar == "F" )
        {
          X     <- F
          Xvec  <- Fvec
          Yeq   <- YeqRefCurve_spf[s,p,]
          Beq   <- BeqRefCurve_spf[s,p,]

          xLab  <- "Fishing Mortality"
        }

        if( indepVar == "U" )
        {
          X     <- U
          Xvec  <- Uvec
          Yeq   <- YeqRefCurve_spf[s,p,]
          Beq   <- BeqRefCurve_spf[s,p,]

          xLab  <- "Harvest rate C /(B+C) (/yr)"
        }


        ubX <- max(which(Yeq >= 0) )

        CXspline    <- splinefun(x = X, y = C)
        YXspline    <- splinefun(x = X, y = Y)
        uCXspline   <- splinefun(x = X, y = uC)
        BXspline    <- splinefun(x = X, y = B)
        UsrXspline  <- splinefun(x = X, y = B-B_usr)

        empXmsy  <- try(uniroot(  interval = c(0.01,maxXspline),
                                  f = YXspline,
                                  deriv = 1 )$root)


        empUsr  <- try(uniroot(   interval = c(0,maxXspline),
                                  f = UsrXspline,
                                  deriv = 0)$root)
        if( class(empXmsy) == "try-error")
          empXmsy <- 0
        if( class(empUsr) == "try-error")
          empUsr <- 0
        
        empMSY  <- YXspline(empXmsy)
        empBmsy <- BXspline(empXmsy)

        empXmsy <- round(empXmsy,2)
        empBmsy <- round(empBmsy,2)
        empMSY  <- round(empMSY,2)

        
        
        lines( x = X, y = Y,
                col = "steelblue", lwd = 2, lty = 1 )
        lines( x = X, y = C,
                col = "steelblue", lwd = 2, lty = 1 )
        lines( x = X, y = uC,
                col = "steelblue", lwd = 2, lty = 1 )
        lines( x = X, y = B,
                col = "black", lwd = 2, lty = 1 )        
      

        legend( "topright", bty = "n",
                legend = c( paste("Umsy = ", empXmsy, sep = "" ),
                            paste(" MSY = ", empMSY, sep = ""),
                            paste("Bmsy = ", empBmsy, sep = ""),
                            paste("Busr = ", round(B_usr,4), sep = ""),
                            paste("Uusr = ", round(empUsr,4), sep = "") ) )


        # Plot some guidelines
        segments(  x0 = empXmsy, x1 = empXmsy,
                    y0 = 0, y1 = empBmsy,
                    lty = 2, col = "red" )
        segments(  x0 = 0, x1 = empXmsy,
                    y0 = empMSY, y1 = empMSY,
                    lty = 2, col = "red" )
        segments(  x0 = 0, x1 = empXmsy,
                    y0 = empBmsy, y1 = empBmsy,
                    lty = 2, col = "red" )

        segments(  x0 = 0, x1 = empUsr,
                    y0 = B_usr, y1 = B_usr,
                    lty = 2, col = "green3" )
        segments(  x0 = empUsr, x1 = empUsr,
                    y0 = 0, y1 = B_usr,
                    lty = 2, col = "green3" )
        
        Xmsy_sp[s,p] <- empXmsy
        Bmsy_sp[s,p] <- empBmsy
        MSY_sp[s,p]  <- empMSY
    }

  mtext( side = 1, text = xLab, outer = T, line = 2)
  mtext( side = 2, text = "Eqbm biomass and catch (kt)", outer = T, line = 2 )



  out <- list(  Xmsy_sp = Xmsy_sp,
                Bmsy_sp = Bmsy_sp,
                MSY_sp  = MSY_sp )

  out
} # END plotEmpYieldCurves()

# plotEmpYieldCurves
# Function to plot median simuated yield 
# resulting from a grid of constant fishing 
# mortality rates - used for checking the
# reference point calculations
plot2dimEmpYieldCurves <- function( sims = 1:150, 
                                    folder = "CC_msyAlloc",
                                    indepVar = "U",
                                    plotYPRcurves = TRUE,
                                    redoEmpRefCurves = FALSE )
{
  nSims <- length(sims)
  blobList <- vector(mode = "list", length = nSims)

  .loadSim(sims[1], folder = folder)

  nS  <- blob$om$nS
  nP  <- blob$om$nP
  nT  <- blob$om$nT
  nF  <- blob$om$nF
  tMP <- blob$om$tMP

  pT <- nT - tMP + 1

  goodReps <- which(blob$goodReps)
  nReps <- length(goodReps)


  # Arrays to hold empirical eqbm yields
  C_isptk       <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )   # Removals (incl dP)
  Y_isptk       <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )   # Yield ( C + K )
  uC_isptk      <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )   # Catch (incl. tptal ponded fish)
  B_isptk       <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )
  eB_isptk      <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )
  tB_isptk      <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )
  F_isptk       <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )
  U_isptk       <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )
  M_isptk       <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )
  R_isptk       <- array( 0, dim = c(nReps,nS,nP,pT,nSims) )
  alloc_k       <- c()

  # Did we really get to allocation we were after?
  
  
  qF_sp <- blob$om$qF_ispft[1,,,2,nT]


  if(!file.exists(here::here("Outputs",folder,"empRefCurves.RData")) |
      redoEmpRefCurves )
  {
    for( x in sims )
    {

      .loadSim(x, folder = folder)

      alloc_k[x] <- blob$ctlList$opMod$projAlloc_f[2]
      for( s in 1:nS)
        for( p in 1:nP)
        {
          # Yield is catch - does this include SOK or ponded fish?
          Y_isptk[1:nReps,s,p,,x]      <- blob$om$C_ispt[goodReps,s,p,tMP:nT]


          if( nF == 5 )
          {
            # Removals, for SP calcs
            C_isptk[1:nReps,s,p,,x]      <- blob$om$C_ispt[goodReps,s,p,tMP:nT]
            uC_isptk[1:nReps,s,p,,x]     <- blob$om$C_ispt[goodReps,s,p,tMP:nT]
          }
            
          if( nF > 5 )
          {
            # Add dead ponded fish to removals for production
            # calculations, and total ponded fish for HR calcs
            uC_ispft  <- C_ispft <- blob$om$C_ispft

            dP_ispft  <- blob$om$dP_ispft
            P_ispft   <- blob$om$P_ispft
            
            C_ispft[,,,6,]  <- P_ispft[,,,6,] * (1 - exp(-0.315))
            uC_ispft[,,,6,] <- P_ispft[,,,6,]

            # if(alloc_k[x] == 0.75)
            # {
            #   F <- apply(X = blob$om$F_ispft[goodReps,s,p,,tMP:nT,drop = FALSE], FUN = sum, MARGIN = c(1,5))
            #   message("alloc = 0.75, F = ", round(F[1],2))
            #   browser()
            # }
            

            C_isptk[1:nReps,s,p,,x] <- apply(X = C_ispft[goodReps,s,p,,tMP:nT,drop = FALSE], FUN = sum, MARGIN = c(1,2,3,5))
            uC_isptk[1:nReps,s,p,,x] <- apply(X = uC_ispft[goodReps,s,p,,tMP:nT,drop = FALSE], FUN = sum, MARGIN = c(1,2,3,5))
          }

          # use end of year SSB so ponded fish are added back in
          B_isptk[1:nReps,s,p,,x]      <- blob$om$SB_ispt[goodReps,s,p,tMP:nT]
          eB_isptk[1:nReps,s,p,,x]     <- blob$om$endSB_ispt[goodReps,s,p,tMP:nT]
          tB_isptk[1:nReps,s,p,,x]     <- blob$om$B_ispt[goodReps,s,p,tMP:nT]
          F_isptk[1:nReps,s,p,,x]      <- apply(X = blob$om$F_ispft[goodReps,s,p,,tMP:nT,drop = FALSE], FUN = sum, MARGIN = c(1,5))
          U_isptk[1:nReps,s,p,,x]      <- uC_isptk[,,,,x]/(uC_isptk[,,,,x] + B_isptk[,,,,x])
          M_isptk[1:nReps,s,p,,x]      <- blob$om$M_iaxspt[goodReps,2,1,s,p,tMP:nT]
          R_isptk[1:nReps,s,p,,x]      <- blob$om$R_ispt[goodReps,s,p,tMP:nT]
          
        }
      
      # Clean up
      message("Seine F = ", round(F_isptk[1,1,1,100,x],2),"; alloc = ",alloc_k[x],"\n")
      gc()
    }

    B_sptk  <- apply(X = B_isptk, FUN = median, MARGIN = 2:5 )
    C_sptk  <- apply(X = C_isptk, FUN = median, MARGIN = 2:5 )
    uC_sptk <- apply(X = uC_isptk, FUN = median, MARGIN = 2:5 )
    Y_sptk  <- apply(X = Y_isptk, FUN = median, MARGIN = 2:5 )
    tB_sptk <- apply(X = tB_isptk, FUN = median, MARGIN = 2:5 )
    eB_sptk <- apply(X = eB_isptk, FUN = median, MARGIN = 2:5 )
    U_sptk  <- apply(X = U_isptk, FUN = median, MARGIN = 2:5 )
    F_sptk  <- apply(X = F_isptk, FUN = median, MARGIN = 2:5 )
    M_sptk  <- apply(X = M_isptk, FUN = median, MARGIN = 2:5 )
    R_sptk  <- apply(X = R_isptk, FUN = median, MARGIN = 2:5 )

    

    C_spk   <- apply(X = C_sptk[,,pT-(150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    uC_spk  <- apply(X = uC_sptk[,,pT-(150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    Y_spk   <- apply(X = Y_sptk[,,pT-(150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    B_spk   <- apply(X = B_sptk[,,pT-(150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    tB_spk  <- apply(X = tB_sptk[,,pT-(150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    eB_spk  <- apply(X = eB_sptk[,,pT-(150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    U_spk   <- apply(X = U_sptk[,,pT-(150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    F_spk   <- apply(X = F_sptk[,,pT - (150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    M_spk   <- apply(X = M_sptk[,,pT - (150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
    R_spk   <- apply(X = R_sptk[,,pT-(150:50),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)

    saveEmpRefCurves <- list( B_isptk     = B_isptk,
                              C_isptk     = C_isptk,
                              uC_isptk    = uC_isptk,
                              tB_isptk    = tB_isptk,
                              U_isptk     = U_isptk,
                              F_isptk     = F_isptk,
                              Y_isptk     = Y_isptk,
                              M_isptk     = M_isptk,
                              R_isptk     = R_isptk,
                              C_sptk      = C_sptk,
                              uC_sptk     = uC_sptk,
                              B_sptk      = B_sptk,
                              tB_sptk     = tB_sptk,
                              eB_sptk     = eB_sptk,
                              F_sptk      = F_sptk,
                              U_sptk      = U_sptk,
                              M_sptk      = M_sptk,
                              R_sptk      = R_sptk,
                              Y_sptk      = Y_sptk,
                              C_spk       = C_spk,
                              uC_spk      = uC_spk,
                              B_spk       = B_spk,
                              tB_spk      = tB_spk,
                              eB_spk      = eB_spk,
                              U_spk       = U_spk,
                              F_spk       = F_spk,
                              R_spk       = R_spk,
                              M_spk       = M_spk,
                              alloc_k     = alloc_k )


    save( saveEmpRefCurves, file = here::here("Outputs",folder,"empRefCurves.RData"))
  } else {
    # Load ref curve list
    load(here::here("Outputs",folder,"empRefCurves.RData"))

    C_sptk  <- saveEmpRefCurves$C_sptk
    uC_sptk <- saveEmpRefCurves$uC_sptk
    Y_sptk  <- saveEmpRefCurves$Y_sptk
    B_sptk  <- saveEmpRefCurves$B_sptk
    tB_sptk <- saveEmpRefCurves$tB_sptk
    F_sptk  <- saveEmpRefCurves$F_sptk
    U_sptk  <- saveEmpRefCurves$U_sptk
    M_sptk  <- saveEmpRefCurves$M_sptk
    R_sptk  <- saveEmpRefCurves$R_sptk
    alloc_k <- saveEmpRefCurves$alloc_k
  
  }
  C_spk   <- apply(X = C_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  uC_spk  <- apply(X = uC_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  Y_spk   <- apply(X = Y_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  B_spk   <- apply(X = B_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  tB_spk  <- apply(X = tB_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  U_spk   <- apply(X = U_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  F_spk   <- apply(X = F_sptk[,,pT - (50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)



  nCurves <- length(unique(alloc_k))

  allAllocs <- unique(alloc_k)
  allAllocs <- allAllocs[order(allAllocs)]

  Xmsy_spc <- array(0, dim = c(nS,nP,nCurves))
  Bmsy_spc <- array(0, dim = c(nS,nP,nCurves))
  MSY_spc <- array(0, dim = c(nS,nP,nCurves))
  

  par(  mfcol = c(nCurves,1), 
        mar = c(.1,1.5,.1,1.5),
        oma = c(5,5,3,3) )

  if( indepVar == "U" )
    maxX <- max(U_spk[1,1,])

  # Loop over histFmult factors
  for( j in 1:nCurves)
  {
    thisAlloc <- allAllocs[j]

    simIdx <- which(alloc_k == thisAlloc)

    plot( x = c(0,maxX), y = c(0,max(Y_spk)),
          type = "n", axes = FALSE,
          xlab = "", ylab = "" )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis( side = 1, cex.axis = 1.5 )
      axis( side = 2, las = 1, cex.axis = 1.5 )
      grid()
      box()

      F   <- F_spk[1,1,simIdx]
      U   <- U_spk[1,1,simIdx]
      Y   <- Y_spk[1,1,simIdx]
      C   <- C_spk[1,1,simIdx]
      uC  <- uC_spk[1,1,simIdx]
      B   <- B_spk[1,1,simIdx]
      # E <- E_spk[1,1,]

      actualOrder <- order(F)

      C   <- C[actualOrder]
      Y   <- Y[actualOrder]
      uC  <- uC[actualOrder]
      B   <- B[actualOrder]
      U   <- U[actualOrder]
      F   <- F[actualOrder]
    

      if(indepVar == "F")
        X <- F
      if(indepVar == "U")
      {
        X <- U
        xLab  <- "Whole fish harvest rate"
      }

      CXspline <- splinefun(x = X, y = Y)
      BXspline <- splinefun(x = X, y = B)

      maxXidx <- max(which(Y > 1e-2))

      empXmsy  <- try(uniroot(  interval = c(0,X[maxXidx]),
                            f = CXspline,
                            deriv = 1 )$root)
      if( class(empXmsy) == "try-error")
        empXmsy <- 0
      
      empMSY  <- CXspline(empXmsy)
      empBmsy <- BXspline(empXmsy)

      empXmsy <- round(empXmsy,2)
      empBmsy <- round(empBmsy,2)
      empMSY  <- round(empMSY,2)
    
      lines( x = X, y = Y,
              col = "steelblue", lwd = 2, lty = 1 )
      

      legend( "topright", bty = "n",cex = 2,
              legend = c( paste( "     ",100*(1 - thisAlloc), "% SOK"),
                          bquote("      "*U[MSY]*" =   " ~ .(format(empXmsy,nsmall=2))),
                          bquote("      "*MSY*" =   " ~ .(format(empMSY,nsmall=2))),
                          bquote("      "*B[MSY]*" = " ~ .(format(empBmsy,nsmall=2))) ) )


      # Plot some guidelines
      segments(  x0 = empXmsy, x1 = empXmsy,
                  y0 = 0, y1 = empMSY,
                  lty = 2, col = "red" )
      segments(  x0 = 0, x1 = empXmsy,
                  y0 = empMSY, y1 = empMSY,
                  lty = 2, col = "red" )
      
      Xmsy_spc[1,1,j] <- empXmsy
      Bmsy_spc[1,1,j] <- empBmsy
      MSY_spc[1,1,j]  <- empMSY
  }

  mtext( side = 1, text = xLab, outer = T, line = 3, cex = 2)
  mtext( side = 2, text = "Equilibrium Yield (kt)", outer = T, line = 2, cex = 2 )



  out <- list(  Xmsy_spc  = Xmsy_spc,
                Bmsy_spc  = Bmsy_spc,
                MSY_spc   = MSY_spc,
                alloc_c   = allAllocs )

  out
} # END plotEmpYieldCurves()

# plotCapVsHarv
# Function to plot median simuated yield 
# resulting from a grid of constant fishing 
# mortality rates - used for checking the
# reference point calculations
plotCapVsHarv <- function(  sims = 1:100, 
                            folder = "PRD_CheckAlloc",
                            indepVar = "U",
                            xLim = c(0,.5),
                            plotYPRcurves = TRUE,
                            redoEmpRefCurves = FALSE )
{
  nSims <- length(sims)
  blobList <- vector(mode = "list", length = nSims)

  .loadSim(sims[1], folder = folder)

  nS  <- blob$om$nS
  nP  <- blob$om$nP
  nT  <- blob$om$nT
  nF  <- blob$om$nF
  tMP <- blob$om$tMP

  pT <- nT - tMP + 1

  goodReps <- which(blob$goodReps)
  nReps <- length(goodReps)

  # Did we really get to allocation we were after?
  # Load ref curve list
  load(here::here("Outputs",folder,"empRefCurves.RData"))

  C_sptk  <- saveEmpRefCurves$C_sptk
  uC_sptk <- saveEmpRefCurves$uC_sptk
  Y_sptk  <- saveEmpRefCurves$Y_sptk
  B_sptk  <- saveEmpRefCurves$B_sptk
  eB_sptk <- saveEmpRefCurves$eB_sptk
  tB_sptk <- saveEmpRefCurves$tB_sptk
  F_sptk  <- saveEmpRefCurves$F_sptk
  U_sptk  <- saveEmpRefCurves$U_sptk
  M_sptk  <- saveEmpRefCurves$M_sptk
  R_sptk  <- saveEmpRefCurves$R_sptk
  alloc_k <- saveEmpRefCurves$alloc_k

  C_spk   <- apply(X = C_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  uC_spk  <- apply(X = uC_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  Y_spk   <- apply(X = Y_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  B_spk   <- apply(X = B_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  eB_spk  <- apply(X = eB_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  tB_spk  <- apply(X = tB_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  U_spk   <- apply(X = U_sptk[,,pT-(50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)
  F_spk   <- apply(X = F_sptk[,,pT - (50:150),,drop = FALSE], FUN = mean, MARGIN = c(1,2,4),na.rm = T)

  H_spk   <- C_spk/(uC_spk + B_spk)


  nCurves <- length(unique(alloc_k))

  allAllocs <- unique(alloc_k)
  allAllocs <- allAllocs[order(allAllocs)]

  Xmsy_spc <- array(0, dim = c(nS,nP,nCurves))
  Bmsy_spc <- array(0, dim = c(nS,nP,nCurves))
  MSY_spc <- array(0, dim = c(nS,nP,nCurves))


  par(  mfcol = c(1,1), 
        mar = c(.1,1.5,.1,1.5),
        oma = c(5,5,3,3) )

  
  plot( x = xLim, y = xLim,
        type = "n", axes = FALSE,
        xlab = "", ylab = "" )
    axis( side = 1 )
    axis( side = 2, las = 1 )
    grid()
    box()

    


  # Loop over alloc factors
  for( j in 1:nCurves)
  {
    thisAlloc <- allAllocs[j]

    simIdx <- which(alloc_k == thisAlloc)

    U <- U_spk[1,1,simIdx]   
    H <- H_spk[1,1,simIdx]   
    F   <- F_spk[1,1,simIdx]
    actualOrder <- order(F)

    U   <- U[actualOrder]
    H   <- H[actualOrder]
    F   <- F[actualOrder]
  

    lines(  x = U, y = H,
            col = j, lwd = 2, lty = 1 )
      
  }

  legTxt <- paste("Seine Allocation = ", allAllocs, sep = "")

  abline(a = 0, b = 1, lty = 2, col = "grey75" )

  legend( "topleft", bty = "n",cex = 2,
          legend = legTxt, lwd = 2,
          col = 1:nCurves )

  mtext( side = 1, text = "Capture Rate", outer = T, line = 3, cex = 2)
  mtext( side = 2, text = "Harvest Rate", outer = T, line = 2, cex = 2 )

} # END plotCapVsHarv()



# compareRefPts()
# Compares OM and conditioning fit reference points
# for a sanity check
compareYieldCurves <- function( obj )
{
  refPtsOM <- obj$rp[[1]]
  refPtsAM <- obj$ctlList$opMod$histRpt$refPts

  # Pull model dimensions
  nS <- obj$om$nS
  nP <- obj$om$nP

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames

  # Now, loop and plot reference points
  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(5,5,3,3) )

  for( s in 1:nS )
    for( p in 1:nP )
    {
      refCurvesOM <- refPtsOM$refCurves
      refCurvesAM <- refPtsAM$refCurves
      plot( x = range( refCurvesOM$F ),
            y = c(0, max(refCurvesOM$Yeq_spf[s,p,], refCurvesOM$Beq_spf[s,p,]  )),
            type = "n", axes = FALSE )
        axis( side = 1 )
        axis( side = 2, las = 1 )

        lines(  x = refCurvesOM$F, y = refCurvesOM$Yeq_spf[s,p,],
                col = "black", lwd = 3 )
        lines(  x = refCurvesOM$F, y = refCurvesOM$Beq_spf[s,p,],
                col = "steelblue", lwd = 3 )

        lines(  x = refCurvesAM$F, y = refCurvesAM$Yeq_spf[s,p,],
                col = "black", lty = 2, lwd = 3 )
        lines(  x = refCurvesAM$F, y = refCurvesAM$Beq_spf[s,p,],
                col = "black", lty = 2, lwd = 3 )


        points( x = refPtsOM$FmsyRefPts$Fmsy_sp[s,p],
                y = refPtsOM$FmsyRefPts$YeqFmsy_sp[s,p], 
                col = "red", pch = 16, cex = 1.5 )
        points( x = refPtsAM$FmsyRefPts$Fmsy_sp[s,p],
                y = refPtsAM$FmsyRefPts$YeqFmsy_sp[s,p], 
                col = "steelblue", pch = 21, cex = 1.5 )

    }
} 


# plotCvsB()
# Fishing mortality as a function of biomass,
# used for showing how well an MP meets the 
# proposed HCR, and for determining the HCR
# implied by the omniscient manager sim
plotCvsB <- function( obj )
{
  goodReps <- obj$goodReps

  SB_ispt   <- obj$om$SB_ispt[goodReps,,,,drop = FALSE]
  C_ispt    <- obj$om$C_ispt[goodReps,,,,drop = FALSE]
  F_ispt    <- obj$om$F_ispft[goodReps,,,2,]
  vB_ispt   <- obj$om$vB_ispft[goodReps,,,2,]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nReps   <- sum(goodReps)

  # Now, let's start plotting. We can fix it later
  speciesNames  <- (obj$om$speciesNames)
  stockNames    <- (obj$om$stockNames)


  projYrs <- tMP:nT

  dep_ispt <- SB_ispt
  Bmsy_sp  <- obj$rp[[1]]$FmsyRefPts$BeqFmsy_sp
  Fmsy_sp  <- obj$rp[[1]]$FmsyRefPts$Fmsy_sp
  MSY_sp   <- obj$rp[[1]]$FmsyRefPts$YeqFmsy_sp

  ctlList <- obj$ctlList

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(5,5,3,3) )

  U_ispt <- C_ispt / vB_ispt

  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      # Convert biomass to Bmsy depletion
      dep_ispt[,s,p,] <- SB_ispt[,s,p,] / Bmsy_sp[s,p]

      # Depletion and F
      maxDep <- max(dep_ispt[,s,p,projYrs], na.rm = T)
      maxC   <- max(C_ispt[,s,p,projYrs]/MSY_sp[s,p], na.rm = T)

      # Plot window
      plot( x = c(0,3), y = c(0,2),
            type = "n", axes = F)
        mfg <- par("mfg")
        # axes
        axis( side = 1 )

        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2, line = 0.2, cex = 1.5 )

        if( mfg[2] == mfg[4] )
        {
          corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          par(xpd = TRUE) #Draw outside plot area
          text( x = corners[2]+0.5, 
                y = mean(corners[3:4]), 
                labels = stockNames[p], srt = 270,
                font = 2, cex = 1.5 )
          par(xpd = FALSE)
        }
        
        axis( side = 2, las = 1 )
        box()

        ptCol <- scales::alpha("grey70", alpha = .3)

        points( x = dep_ispt[,s,p,projYrs], y = C_ispt[,s,p,projYrs]/MSY_sp[s,p],
                col = ptCol, pch = 1 ) 

        for( i in 1:nReps )
        { 
          # Plot a smoother for each replicate
          lineCol <- scales::alpha("red", alpha = .3)
          smoother <- loess.smooth( x = dep_ispt[i,s,p,projYrs], 
                                    y = C_ispt[i,s,p,projYrs]/MSY_sp[s,p] )
          lines( smoother, col = "grey30", lwd = .8 )
          flB <- dep_ispt[i,s,p,projYrs[c(1,length(projYrs))]]
          flC <- C_ispt[i,s,p,projYrs[c(1,length(projYrs))]]/MSY_sp[s,p]
          points( x = flB,
                  y = flC,
                  col = c("blue","red"), cex = .5,
                  pch = 16 )
        }


        abline( v = 1, lty = 2, lwd = .8)
        abline( h = 1, lty = 2, lwd = .8)

    }
  }

  mtext( side = 1, text = expression(B/B[MSY]), outer = TRUE, line = 2, font = 2)
  mtext( side = 2, text = expression(C/MSY), outer = TRUE, line = 2, font = 2)

}


# plotFvsB()
# Fishing mortality as a function of biomass,
# used for showing how well an MP meets the 
# proposed HCR, and for determining the HCR
# implied by the omniscient manager sim
plotFvsB <- function( obj )
{
  goodReps <- obj$goodReps

  SB_ispt   <- obj$om$SB_ispt[goodReps,,,,drop = FALSE]
  C_ispt    <- obj$om$C_ispt[goodReps,,,,drop = FALSE]
  F_ispt    <- obj$om$F_ispft[goodReps,,,2,]
  vB_ispt   <- obj$om$vB_ispft[goodReps,,,2,]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nReps   <- sum(goodReps)

  # Now, let's start plotting. We can fix it later
  speciesNames  <- (obj$om$speciesNames)
  stockNames    <- (obj$om$stockNames)


  projYrs <- tMP:nT

  dep_ispt <- SB_ispt
  Bmsy_sp  <- obj$rp[[1]]$FmsyRefPts$BeqFmsy_sp
  Fmsy_sp  <- obj$rp[[1]]$FmsyRefPts$Fmsy_sp

  ctlList <- obj$ctlList

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,2.5),
        oma = c(5,5,3,3) )

  U_ispt <- C_ispt / vB_ispt

  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      # Convert biomass to Bmsy depletion
      dep_ispt[,s,p,] <- SB_ispt[,s,p,] / Bmsy_sp[s,p]

      # Depletion and F
      maxDep <- max(dep_ispt[,s,p,projYrs], na.rm = T)
      maxF   <- max(F_ispt[,s,p,projYrs]/Fmsy_sp[s,p], na.rm = T)

      # Plot window
      plot( x = c(0,3), y = c(0,2),
            type = "n", axes = F)
        mfg <- par("mfg")
        # axes
        axis(side =1)

        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2, line = 0.2, cex = 1.5 )

        if( mfg[2] == mfg[4] )
        {
          corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          par(xpd = TRUE) #Draw outside plot area
          text( x = corners[2]+0.5, 
                y = mean(corners[3:4]), 
                labels = stockNames[p], srt = 270,
                font = 2, cex = 1.5 )
          par(xpd = FALSE)
        }
        
        axis( side = 2, las = 1 )
        box()

        ptCol <- scales::alpha("grey70", alpha = .3)

        points( x = dep_ispt[,s,p,projYrs], y = F_ispt[,s,p,projYrs]/Fmsy_sp[s,p],
                col = ptCol, pch = 1 ) 

        for( i in 1:nReps )
        { 
          # Plot a smoother for each replicate
          lineCol <- scales::alpha("red", alpha = .3)
          smoother <- loess.smooth( x = dep_ispt[i,s,p,projYrs], 
                                    y = F_ispt[i,s,p,projYrs]/Fmsy_sp[s,p] )
          lines( smoother, col = "grey30", lwd = .8 )
          flB <- dep_ispt[i,s,p,projYrs[c(1,length(projYrs))]]
          flF <- F_ispt[i,s,p,projYrs[c(1,length(projYrs))]]/Fmsy_sp[s,p]
          points( x = flB,
                  y = flF,
                  col = c("blue","red"), cex = .5,
                  pch = 16 )
        }


        abline( v = 1, lty = 2, lwd = .8)
        abline( h = 1, lty = 2, lwd = .8)

    }
  }

  mtext( side = 1, text = expression(B/B[MSY]), outer = TRUE, line = 2, font = 2)
  mtext( side = 2, text = expression(F/F[MSY]), outer = TRUE, line = 2, font = 2)

}

# plotBatchPerf_sp()
# Plots multipanels (faceted by species/stock)
# of performance statistics on the y axis, with respect
# to the OM grid on the x axis
# CAUTION: currently only works for numeric xAxis and yAxis
plotBatchPerf_sp <- function( batchFolder = "fourthBatch",
                              xAxis = "projObsErrMult",
                              yAxis = "PBtGt.8Bmsy",
                              yRangeIn = NULL )
{
  # First load full stats table from the batch folder
  statsFile <- here::here(  "Outputs",batchFolder,"statistics",
                            "fullStatTable.csv")
  statTable <- read.csv(statsFile, header = TRUE, stringsAsFactors = FALSE )

  # Now, let's start plotting. We can fix it later
  speciesNames  <- unique(statTable$species)
  stockNames    <- unique(statTable$stock)

  scenarios <- unique(statTable$scenario)
  mps       <- unique(statTable$mp)

  nS <- length(speciesNames)
  nP <- length(stockNames)

  nScen <- length(scenarios)
  nMPs  <- length(mps)



  xLabs <- unique(statTable[,xAxis])
  xLabs <- xLabs[order(xLabs)]
  xMax  <- max(xLabs)
  xMin  <- min(xLabs)
  
  cols <- RColorBrewer::brewer.pal(n = nMPs, "Dark2")


  par(  mfcol = c( nP, nS ), 
        mar = c( 0,2,0,2 ), 
        oma = c(5,7,3,4) )

  yRange <- yRangeIn

  mpJitter <- seq( from = -.3, to = .3, length = nMPs )

  xDiff <- mean(diff(xLabs))
  mpJitter <- mpJitter * xDiff


  for( s in 1:nS )
    for( p in 1:nP )
    {
      subTable <- statTable %>%
                  filter( species == speciesNames[s],
                          stock == stockNames[p] )
      
      if(is.null(yRangeIn))            
        yRange <- range(subTable[,yAxis])      
      
      
      plot( x = c(-.5 + xMin,xMax + .5), y = yRange,
            type = "n", axes = FALSE )
        mfg <- par("mfg")
        # axes
        if( mfg[1] == mfg[3])
          axis( side = 1, at = xLabs )

        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
        
        axis( side = 2, las = 1 )
        
        box()
        # Add vertical grid lines to split OMs
        grid()

        # Now, plot the stuff
        for( m in 1:nMPs )
        {
          mpTable <- subTable %>% filter( mp == mps[m] )
          mpTable <- mpTable[order(mpTable[,xAxis]),]

          points( x = mpTable[,xAxis] + mpJitter[m],
                  y = mpTable[,yAxis],
                  col = cols[m], pch = 14 + m )
          lines(  x = mpTable[,xAxis] + mpJitter[m],
                  y = mpTable[,yAxis],
                  col = cols[m], lwd = .8 )
        }

        if( mfg[2] == mfg[4] )
        {
          corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          par(xpd = TRUE) #Draw outside plot area
          text( x = corners[2]+0.2, 
                y = mean(corners[3:4]), 
                labels = stockNames[p], srt = 270,
                font = 2, cex = 1.5 )
          par(xpd = FALSE)
        }



    }
  legend( "topleft",
          legend = mps,
          col = cols,
          pch = 14 + 1:nMPs,
          lwd = 1 )

  

  mtext( side = 1, text = xAxis,
         outer = TRUE, cex = 1.5, 
         font = 2, line = 3 )

  mtext( side = 2, text = yAxis,
         outer = TRUE, cex = 1.5, font = 2, line = 2 )



}

# plotHCRex()
# Plots a generic hockey-stick
# harvest control rule, with given
# lower, upper control points, and
# low/high Fs. Stock status is calculated
# as a proportion of Bmsy, and fishing
# mortality rate as a proportion of
# Fmsy.
plotHCRex <- function(  LCP = .3,
                        UCP = .6,
                        lowF = .0,
                        highF = 1,
                        maxSeq = seq(from = 0.02, to = 0.2, by = 0.02),
                        language = "English")
{
  x <- seq(0,1.3, length.out = 100 )
  y <- rep(lowF, length(x))


  par( mar = c(4,5,1,1), oma = c(2,2,1,1) )

  plot( x = range(x), y = c(0,1.2*max(maxSeq)), type = "n",
        las = 1,
        xlab = "",
        ylab = "",
        cex.lab = 1.5,
        cex.axis = 1.5,
        axes = FALSE)
  if(language=="English"){
    mtext(side = 2, text = "Target Harvest Rate", line = 4, cex = 1.5 )
    mtext(side = 1, text = expression(B/B[0]), line = 4, cex = 1.5 )
    box()
    axis(side=1, at=seq(0.0, max(x), by=0.2), labels = sprintf("%.1f", seq(0.0, max(x), by=0.2)))
    axis(side=2, at=seq(0.0, 0.2, by=0.02), labels = sprintf("%.2f", seq(0.0, 0.2, by=0.02)), las=1)
    # axis(side=4, at=seq(0.0, 0.2, by=0.02), labels = sprintf("%.2f", seq(0.0, 0.2, by=0.02)), las=1)
  }else{ #French
    mtext(side = 2, text = "Taux de récolte cible", line = 4, cex = 1.5 )
    mtext(side = 1, text = expression(B/B[0]), line = 4, cex = 1.5 )
    box()
    axis(side=1, at=seq(0.0, max(x), by=0.2), labels = chartr(".", ",", sprintf("%.1f", seq(0.0, max(x), by=0.2))))
    axis(side=2, at=seq(0.0, 0.2, by=0.02), labels = chartr(".", ",", sprintf("%.2f", seq(0.0, 0.2, by=0.02))), las=1)
  }
    segments( x0 = 0, x1 = LCP,
              y0 = lowF, y1 = lowF,
              col = "grey50", lwd = 1.5 )
    for( k in 1:length(maxSeq))
    {
      lineWd <- 1
      lineCol <- "grey50"
      if( maxSeq[k] %in% c(0.1))
      {
        lineWd <- 3
        lineCol <- "black"

        if(maxSeq[k] == 0.064)
          lineCol <- "royalblue"
      }
      segments( x0 = LCP, x1 = UCP,
                y0 = lowF, y1 = highF*maxSeq[k],
                col = lineCol, lwd = lineWd )
      segments( x0 = UCP, x1 = max(x),
                y0 = highF*maxSeq[k], y1 = highF*maxSeq[k],
                col = lineCol, lwd = lineWd )
    }
    abline( v = c(LCP, UCP),
            col = c("red","darkgreen"),
            lty = 2, lwd = 2 )
    # abline( h = highF*maxSeq[1], lty = 2, col = "grey70" )


} # END plotHCRex()




# plotHCR()
# Plots a generic hockey-stick
# harvest control rule, with given
# lower, upper control points, and
# low/high Fs. Stock status is calculated
# as a proportion of Bmsy, and fishing
# mortality rate as a proportion of
# Fmsy.
plotHCR <- function(  LCP = .4,
                      UCP = .6,
                      lowF = .1,
                      highF = 1 )
{
  x <- seq(0,1.3, length.out = 100 )
  y <- rep(lowF, length(x))


  par( mar = c(4,5,1,1), oma = c(2,2,1,1) )

  plot( x = range(x), y = c(0,1.2), type = "n",
        las = 1,
        xlab = expression(B/B[MSY]),
        ylab = expression(F/F[MSY]),
        cex.lab = 1.5,
        cex.axis = 1.5 )
    segments( x0 = 0, x1 = LCP,
              y0 = lowF, y1 = lowF,
              col = "grey50", lwd = 3 )
    segments( x0 = LCP, x1 = UCP,
              y0 = lowF, y1 = highF,
              col = "grey50", lwd = 3 )
    segments( x0 = UCP, x1 = max(x),
              y0 = highF, y1 = highF,
              col = "grey50", lwd = 3 )
    abline( v = c(LCP, UCP),
            col = c("red","darkgreen"),
            lty = 2, lwd = 2 )
    abline( h = highF, lty = 2, col = "grey70" )


} # END plotHCR()

plotUtilityFunction <- function( )
{
  x <- seq(1e-1,1, length = 1000)

  y <- (5*x - 1)/4/x

  plot( x = range(x), y = range(y),
        xlab = expression(C[spt]/TAC[spt]),
        ylab = "Utility", type = "n", las = 1 )
    mtext( side = 3, text = "TAC utilisation utility", font = 2)
    lines( x = x, y = y, lwd = 3, col = "grey40" )
    abline(h = c(0,1), lty = 2, col = "black" )
    abline( v = .2, lty = 2, col = "red", lwd = .8 )

}


# plotConvStats()
plotConvStats <- function( obj = blob )
{
  goodReps      <- obj$goodReps

  # Pull max gradient value and hessian indicator
  maxGrad_itsp  <- obj$mp$assess$maxGrad_itsp[goodReps,,,,drop = FALSE]
  pdHess_itsp   <- obj$mp$assess$pdHess_itsp[goodReps,,,,drop = FALSE]

  nReps <- dim(maxGrad_itsp)[1]

  # Now we want to get the mean and SD
  # of these values over the replicates
  quantsMaxGrad_qtsp <- apply(  X = maxGrad_itsp, FUN = quantile,
                                MARGIN = 2:4, probs = c(0.05, 0.5, 0.95),
                                na.rm = T )

  propPDHess_tsp  <- apply( X = pdHess_itsp, FUN = mean,
                            MARGIN = 2:4, na.rm = T )

  pT <- obj$ctlList$opMod$pT
  nS <- obj$om$nS
  nP <- obj$om$nP

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames

  xJitter     <- seq( from = -.3, to = .3, length.out = nP )
  stockCols   <- RColorBrewer::brewer.pal( nP, "Dark2" )
  stockPts    <- seq( from = 21, by = 1, length.out = nP )
  rectWidth   <- .6 / nP

  par( mfcol = c(nS,2), mar = c(0,2,0,1), oma = c(3,3,3,2) )

  # First, plot the maxGrads, using 
  # different colours for each species, and 
  # pch for each stock... Jitter!
  for( s in 1:nS )
  {
    plot( x = c(1,pT), 
          y = range(quantsMaxGrad_qtsp, na.rm = T),
          type = "n",
          axes = FALSE )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3])
        axis( side = 1 )
      if(mfg[1] == 1 )
        mtext( side = 3, text = "Max Gradient Component")
      
      axis( side = 2, las = 1 )
      box()
      grid()

      for( p in 1:nP )
      {
        points( x = 1:pT + xJitter[p], y = quantsMaxGrad_qtsp[2,,s,p],
                col = stockCols[p], pch = stockPts[p], bg = stockCols[p] )
        segments( x0 = 1:pT + xJitter[p], x1 = 1:pT + xJitter[p],
                  y0 = quantsMaxGrad_qtsp[1,,s,p],
                  y1 = quantsMaxGrad_qtsp[3,,s,p],
                  col = stockCols[p], lty = 1 )
      }
      if( s == 1 & !is.null(stockNames) )
        legend( "topright", bty = "n",
                legend = stockNames,
                col = stockCols,
                pch = stockPts,
                pt.bg = stockCols )
  }

  # now do proportion of PD hessians
  for( s in 1:nS )
  {
    plot( x = c(0,pT+1), y = c(0,1.3),
          axes = FALSE, type = "n" )
      # Axes
      if( mfg[1] == mfg[3])
        axis( side = 1 )
      axis( side = 2, las = 1 )
      mtext( side = 4, text = speciesNames[s], font = 2, line = 2)
      if(mfg[1] == 1 )
        mtext( side = 3, text = "Proportion of PD Hessians")

      box()
      abline( h = 1.0, lty = 2, lwd = 2, col = "grey40" )

      # Plot rectangles of PD Hessians
      for( p in 1:nP )
      {
        rect( xleft = 1:pT + xJitter[p] - rectWidth/2,
              xright = 1:pT + xJitter[p] + rectWidth/2,
              ybottom = 0,
              ytop = propPDHess_tsp[,s,p],
              col = stockCols[p] )
      }
      if( s == 1 & (!is.null(stockNames)) )
        legend( "topright", bty = "n",
                legend = stockNames,
                col = stockCols,
                pch = 22, pt.bg = stockCols )

  }

}

# plotRefPtSeries()
# A copy-forward of the mseR plot of the same
# name. Shows a time-series of control points
# as biomass estimates, true SSB, and target
# vs effective harvest rates.
plotRefPtSeries <- function(  obj = blob,
                              iRep = 1, 
                              sIdx = 1, pIdx = 1 )
{
  # Get dims
  tMP       <- obj$om$tMP
  pT        <- obj$ctlList$opMod$pT
  nT        <- obj$om$nT
  fYear     <- obj$ctlList$opMod$fYear

  fleetType_f <- obj$om$fleetType_f
  sokFleets   <- which(fleetType_f %in% 2:3)

  # Pull SSB and U
  SB_t      <- obj$om$SB_ispt[iRep,sIdx,pIdx,]
  TAC_t     <- obj$mp$hcr$TAC_ispt[iRep,sIdx,pIdx,]

  # Gotta pull ponded fish and landings together
  C_t       <- apply(X = obj$om$C_ispft[iRep,sIdx,pIdx,fleetType_f == 1,,drop =FALSE], FUN = sum, MARGIN = 5)
  P_t       <- apply(X = obj$om$P_ispft[iRep,sIdx,pIdx,fleetType_f %in% c(2,3),,drop =FALSE], FUN = sum, MARGIN = 5)

  effU_t    <- TAC_t/(SB_t + TAC_t)
  effU2_t   <- C_t/(SB_t + C_t)
  targU_t   <- obj$mp$hcr$targetF_ispt[iRep,sIdx,pIdx,]



  # Now pull retrospective values
  if(pT > 1)
    retroSB_t     <- diag(obj$mp$assess$retroSB_itspt[iRep,1:pT,sIdx,pIdx,tMP:nT])
  else retroSB_t <- obj$mp$assess$retroSB_itspt[iRep,1:pT,sIdx,pIdx,tMP:nT]
  
  retroB0_t     <- obj$mp$assess$retroB0_itsp[iRep,,sIdx,pIdx]

  if(all(is.na(retroB0_t)))
    retroB0_t <- rep(obj$rp[[iRep]]$B0_sp[1,sIdx,pIdx],pT)


  # Need CPs and Fref/targetF
  LCP_t     <- obj$mp$hcr$LCP_ispt[iRep,sIdx,pIdx,]
  UCP_t     <- obj$mp$hcr$UCP_ispt[iRep,sIdx,pIdx,]
  Fref_t    <- obj$mp$hcr$Fref_ispt[iRep,sIdx,pIdx,]
  Bref_t    <- obj$mp$hcr$Bref_ispt[iRep,sIdx,pIdx,]

  
  # Full years
  yrs     <- seq( from = fYear, by = 1, length.out = nT)
  par(  mfrow = c(2,1),
        oma = c(3,3,2,2),
        mar = c(.1,2,.1,2) )
  # First plot B series, over projection period
  maxB <- max(SB_t[tMP:nT],UCP_t[tMP:nT],retroSB_t,na.rm = T)
  plot( x = range(yrs[tMP:nT]), y = c(0,maxB),
        axes = FALSE, type = "n")
    axis( side = 2, las = 1)
    mtext( side = 2, text = "OM and AM Biomass (kt)", line = 3)
    grid()
    box()

    lines(  x = yrs, y = SB_t, col = "salmon",
            lwd = 2 )
    points( x = yrs, y = LCP_t, pch = 16, col = "red" )
    points( x = yrs, y = UCP_t, pch = 16, col = "darkgreen" )
    points( x = yrs[tMP:nT], y = retroB0_t, pch = 16, col = "royalblue" )
    lines(  x = yrs[tMP:nT], y = retroSB_t, col = "grey40", lwd = 2 )
    legend( x = "topright",
            bty = "n",
            legend = c("OM SSB", "AM Forecast","UCP","LCP","AM B0"),
            lwd = c(2,2,NA,NA,NA),
            pch = c(NA,NA,16,16,16),
            col = c("salmon","grey40","darkgreen","red","royalblue"))

  maxU <- max(effU_t[tMP:nT],effU2_t[tMP:nT],Fref_t,na.rm = T)

  plot( x = range(yrs[tMP:nT]), y = c(0,maxU), 
        axes = FALSE, type = "n")
    axis( side = 2, las = 1)
    axis( side = 1 )
    mtext( side = 1, text = "Year", line = 3)
    mtext( side = 2, text = "Harvest Rate", line = 3)
    grid()
    box()
    segments( x0 = yrs,
              y1 = Fref_t,
              y0 = targU_t,
              col = "grey40",
              lwd = 2 )
    lines( x = yrs, y = effU2_t, lwd = 2, col = "royalblue")
    # lines( x = yrs, y = legU_t, lwd = 1, col = "royalblue")
    points( x = yrs, y = Fref_t,
            col = "salmon", pch = 16 )
    points( x = yrs, y = targU_t,
            col = "salmon", pch = 21, bg = NA )

    legend( x = "topright",
            bty = "n",
            legend = c("Effective HR","Target HR","Reference HR"),
            lwd = c(2,NA,NA),
            pch = c(NA,21,16),
            col = c("royalblue","salmon","salmon"),
            pt.bg = c(NA,NA,NA) )

}


plotTulipHR <- function( obj = blob, nTrace = 3 )
{
  # Model dimensions
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF

  # Good replicates
  goodReps <- obj$goodReps

  # Pull target HR
  Uref_p <- obj$ctlList$mp$hcr$Uref_p

  fleetType_f <- obj$om$fleetType_f
  sokFleets   <- which(fleetType_f %in% 2:3)
  predGears   <- obj$ctlList$opMod$predGears
  
  if(is.null(predGears))
    predGears <- c()

  commGears   <- (1:nF)[!1:nF %in% predGears]
    
  # Catch
  C_ispft  <- obj$om$C_ispft[goodReps,,,,,drop = FALSE]
  C_ispt       <- apply(X = C_ispft[,,,-sokFleets,,drop = FALSE], FUN = sum, MARGIN = c(1:3,5))

  predC_ispt   <- apply(X = C_ispft[,,,predGears[!predGears %in% sokFleets],,drop = FALSE], FUN = sum, MARGIN = c(1:3,5))
  commC_ispt   <- apply(X = C_ispft[,,,commGears[!commGears %in% sokFleets],,drop = FALSE], FUN = sum, MARGIN = c(1:3,5))

  # Ponded fish
  P_ispft     <- obj$om$P_ispft[goodReps,,,,,drop = FALSE]
  P_ispt      <- apply(X = obj$om$P_ispft[goodReps,,,sokFleets,,drop = FALSE], FUN=sum, MARGIN=c(1:3,5))
  predP_ispt   <- apply(P_ispft[,,,intersect(sokFleets,predGears),,drop=FALSE], FUN=sum, MARGIN=c(1,2,3,5))
  commP_ispt   <- apply(P_ispft[,,,intersect(sokFleets,commGears),,drop=FALSE], FUN=sum, MARGIN=c(1,2,3,5))



  SB_ispt  <- obj$om$SB_ispt[goodReps,,,,drop = FALSE]
  TAC_ispt <- obj$mp$hcr$TAC_ispt[goodReps,,,,drop = FALSE]
  effHR_ispt <- TAC_ispt/(SB_ispt + TAC_ispt)

  targetF_ispt <- obj$mp$hcr$targetF_ispt

  # Calculate aggregate harvest rate
  commU_ispt   <- commC_ispt/(SB_ispt + commC_ispt + commP_ispt + predP_ispt + predC_ispt)
  predU_ispt   <- predC_ispt/(SB_ispt + commC_ispt + commP_ispt + predP_ispt + predC_ispt)

  # Calculate fleet specific harvest rate
  U_ispft <- array(NA, dim=c(sum(goodReps),nS,nP,nF,nT))
  for (s in 1:nS)
    for(f in 1:nF)
    {
      if(f %in% sokFleets)
        U_ispft[,s,,f,] <- P_ispft[,s,,f,]/(SB_ispt[,s,,] + C_ispt[,s,,] + P_ispt[,s,,])

      if(!f %in% sokFleets )
        U_ispft[,s,,f,] <- C_ispft[,s,,f,]/(SB_ispt[,s,,] + C_ispt[,s,,] + P_ispt[,s,,])        
    }

  # Harvest rate envelopes
  U_qspft <- apply(  X = U_ispft, FUN = quantile,
                    MARGIN = c(2,3,4,5), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  # Quantile for aggregate harvest rate
  predU_qspt <- apply(  X = predU_ispt, FUN = quantile,
                        MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                        na.rm = T )

  commU_qspt <- apply(  X = commU_ispt, FUN = quantile,
                        MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                        na.rm = T )

  # effective U
  effU_qspt <- apply(  X = effHR_ispt, FUN = quantile,
                        MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                        na.rm = T )

  targU_qspt <- apply(  X = targetF_ispt, FUN = quantile,
                        MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                        na.rm = T )


  # fishing fleets
  fishG  <- obj$ctlList$opMod$commGears
  fleets <- obj$ctlList$opMod$fleets
  fColrs <- brewer.pal(nF, 'Dark2')

  nReps   <- dim(C_ispft)[1]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,2),
        oma = c(5,5,3,3) )

  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(predU_qspt[,s,p,], commU_qspt[,s,p,], na.rm = T) ),
            type = "n", axes = FALSE)

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
        mtext( side = 4, text = stockNames[p], line = 2)
      
      box()
      grid()
      
      # plot one polygon for aggregate comm and pred U
      polygon(  x = c(yrs, rev(yrs)),
                y = c(commU_qspt[1,s,p,], rev(commU_qspt[3,s,p,])),
                col = "grey65", border = NA )
      
      lines( x = yrs, y = commU_qspt[2,s,p,], lwd = 1, lty=2 )

      polygon(  x = c(yrs, rev(yrs)),
                y = c(predU_qspt[1,s,p,], rev(predU_qspt[3,s,p,])),
                col = scales::alpha("darkgreen",.5), border = NA )
      # lines( x = yrs, y = predU_qspt[2,s,p,], lwd = 1, lty=2,
      #        col = "darkgreen" )

      # plot individual lines for each fleet
      for( f in fishG )
      {
        
        if(sum(U_qspft[2,s,p,f,], na.rm=T)==0)
          next()

        lines( x = yrs, y = U_qspft[2,s,p,f,], col=fColrs[f], lwd = 1.5)
        # for( tIdx in traces )
        #   lines( x = yrs, y = C_ispft[tIdx,s,p,f,], lwd = .8 )
      }
      # lines(x = yrs, y = targU_qspt[2,s,p,], lwd = 3)

      # abline(h=0, lwd=1.5)
      abline( v = yrs[tMP], col = "grey30", lty = 3 )
      segments( x0 = yrs[tMP], x1 = yrs[nT],
                y0 = Uref_p[p], lty = 2, lwd = 1, col = "red")

      legend('topright', bty='n', cex=0.8,
              legend=c('total HR', fleets[fishG],'Ref. HR'),
              lwd=c(1, rep(1,length(fishG)),1),
              lty=c(2, rep(1,length(fishG)),3),
              col=c('black',fColrs[fishG],'red')
              )

    }
  }

  mtext( side = 2, outer = TRUE, text = "Harvest Rate",
          line = 2, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)
}




plotTulipF <- function( obj = blob, nTrace = 3 )
{
  # Model dimensions
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF

  # Good replicates
  goodReps <- obj$goodReps

  # Pull reference points
  Fmsy_sp <- obj$rp[[1]]$FmsyRefPts$Fmsy_sp
    
  # Fishing mortality series
  F_ispft <- obj$om$F_ispft[goodReps,,,,,drop = FALSE]

  # Fishing mortality envelopes
  F_qspft <- apply(  X = F_ispft, FUN = quantile,
                    MARGIN = c(2,3,4,5), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  # Aggregate fishing mortality
  F_ispt <- apply(  X = F_ispft, FUN = sum,
                    MARGIN = c(1,2,3,5),
                    na.rm = T )
  F_qspt <- apply(  X = F_ispt, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  # fishing fleets
  fishG  <- obj$ctlList$opMod$commGears
  fleets <- obj$ctlList$opMod$fleets
  fColrs <- brewer.pal(nF, 'Dark2')
  

  nReps   <- dim(F_ispft)[1]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,2),
        oma = c(5,5,3,3) )


  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(F_qspt[,s,p,], na.rm = T) ),
            type = "n", axes = F, ylim=c(0,max(F_qspt[,s,p,], na.rm = T) ))

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }

      box()
      grid()
      
      # plot one polygon for aggregate F
      polygon(  x = c(yrs, rev(yrs)),
                y = c(F_qspt[1,s,p,], rev(F_qspt[3,s,p,])),
                col = "grey65", border = NA )
      # lines( x = yrs, y = F_qspt[2,s,p,], lwd = 2 )

      # plot individual lines for each fleet
      for( f in fishG )
      {
        
        if(sum(F_qspft[2,s,p,f,], na.rm=T)==0)
          next()

        lines( x = yrs, y = F_qspft[2,s,p,f,], col=fColrs[f], lwd = 1.5)
        # for( tIdx in traces )
        #   lines( x = yrs, y = F_ispft[tIdx,s,p,f,], lwd = .8 )
      }

      # abline(h=0, lwd=1.5)
      abline( v = yrs[tMP], col = "grey30", lty = 3 )
      abline( h = Fmsy_sp[s,p], lty = 2, lwd = 1, col = "red")

      legend('topright', bty='n', cex=0.8,
              legend=c('total F', fleets[fishG],'Fmsy'),
              lwd=c(2, rep(1,length(fishG)),1),
              lty=c(1, rep(1,length(fishG)),3),
              col=c('black',fColrs[fishG],'red')
              )

    }
  }

  mtext( side = 2, outer = TRUE, text = "Fishing mortality (/yr)",
          line = 2, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)
}

# TAC utilisation envelopes
plotTulipTACu <- function( obj = blob, nTrace = 3, fIdx=1 )
{

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT

  goodReps <- obj$goodReps
  

  # Get catch for trawl fleet, in projections only
  C_ispft     <- obj$om$C_ispft[goodReps,,,2,tMP:nT,drop = FALSE]
  TAC_ispft   <- obj$mp$hcr$TAC_ispft[goodReps,,,2,tMP:nT,drop = FALSE]

  nReps   <- dim(TAC_ispft)[1]

  TACu_ispft <- C_ispft / TAC_ispft


  TACu_qspt <- apply( X = TACu_ispft, FUN = quantile,
                      MARGIN = c(2,3,5), probs = c(0.025, 0.5, 0.975),
                      na.rm = T )

  TACu_ispt <- array(0, dim = c(nReps,nS,nP,(nT - tMP + 1)))
  TACu_ispt[1:nReps,,,] <- TACu_ispft[1:nReps,,,fIdx,]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)
  yrs <- yrs[tMP:nT]


  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,2),
        oma = c(5,5,3,3) )


  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(TACu_qspt, na.rm = T) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      polygon(  x = c(yrs, rev(yrs)),
                y = c(TACu_qspt[1,s,p,], rev(TACu_qspt[3,s,p,])),
                col = "grey65", border = NA )
      lines( x = yrs, y = TACu_qspt[2,s,p,], lwd = 3 )
      for( tIdx in traces )
        lines( x = yrs, y = TACu_ispt[tIdx,s,p,], lwd = .8 )

      abline( v = yrs[tMP], col = "grey30", lty = 3 )
    }
  }
  mtext( side = 2, outer = TRUE, text = expression(C[spt]/TAC[spt]),
          line = 2, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)

} # END plotTulipTACu


# plotDistN_SOK()
# Distribution of N_SOK values across
# all simulation replicates.
plotDistN_SOK <- function( obj = blob )
{
  goodReps <- obj$goodReps

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nReps   <- length(goodReps)

  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]

  stockNames <- c(stockNames,"Agg")

  fleetType_f <- obj$om$fleetType_f
  sokFleets <- which(fleetType_f == 2)

  P_ispft   <- obj$om$P_ispft[goodReps,,,sokFleets,tMP:nT,drop = FALSE]
  P_ispt    <- apply(X = P_ispft[,1,,,,drop = FALSE], FUN = sum, MARGIN = c(1,2,3,5) )

  P_it      <- apply(X = P_ispt, FUN = sum, MARGIN = c(1,4))

  P_ispt[P_ispt > 0] <- 1
  P_it[P_it > 0] <- 1

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")


  # Get granular dist
  Nsok_isp <- apply( X = P_ispt, FUN = sum, MARGIN = c(1,2,3), na.rm = T)

  Nsok_qsp <- apply(X = Nsok_isp, FUN = quantile, probs = c(0.025,.25, .5, .75, .975),
                      MARGIN = c(2,3))

  # Get aggregate distribution as well
  Nsok_i <- apply(X = P_it, FUN = sum, MARGIN = 1)
  Nsok_q <- quantile(Nsok_i, probs = c(.025, .25, .5, .75, .975) )

  # Make a bar plot

  plot( x = c(.5,4.5), y = c(0,max(Nsok_i, Nsok_isp)),
        type = "n", axes = FALSE, xlab = "Population", ylab = "Number of years with SOK fishing" )
    axis( side = 1, at = 1:4, labels = stockNames )
    axis( side = 2, las = 1 )
    box()
    grid()
    for( p in 1:nP )
    {
      # rect( xleft = p - .3, xright= p + .3,
      #       ybottom = 0, ytop = Nsok_qsp[2,1,p], border = NA,
      #       col = "grey60" )
      
      segments( x0 = p, x1 = p, 
                y0 = Nsok_qsp[1,1,p],
                y1 = Nsok_qsp[5,1,p], col = "black", lwd = 2 )
      rect( xleft = p-.1, xright = p+.1, 
            ybottom = Nsok_qsp[2,1,p], border = "black",
            ytop = Nsok_qsp[4,1,p], col = "white", lwd = 2 )
      segments( x0 = p-.1, x1 = p + .1, y0 = Nsok_qsp[3,1,p], col = "black", lwd = 2 )
    }
    segments( x0 = nP + 1, x1 = nP + 1, 
              y0 = Nsok_q[1],
              y1 = Nsok_q[5], col = "black", lwd = 2 )
    rect( xleft = nP + 1 - .1, xright= nP + 1 + .1,
          ybottom = Nsok_q[2], ytop = Nsok_q[4], border = "black",
          col = "white", lwd = 2 )
    segments( x0 = nP + 1-.1, x1 = nP + 1 + .1, y0 = Nsok_q[3], col = "black", lwd = 2 )


    mtext(  side = 1, adj = .9, line = 2, cex = .6,
            text = stamp, col = "grey60" )

} # END plotDistN_SOK

plotTulipPondedFish <- function(  obj = blob,
                                  nTrace = 3,
                                  traces = NULL,
                                  leg = TRUE,
                                  proj = FALSE )
{
  goodReps <- obj$goodReps

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nReps   <- length(goodReps)

  fleetType_f <- obj$om$fleetType_f
  sokFleets <- which(fleetType_f == 2)

  P_ispft   <- obj$om$P_ispft[goodReps,,,sokFleets,,drop = FALSE]
  P_ispt    <- apply(X = P_ispft[,1,,,,drop = FALSE], FUN = sum, MARGIN = c(1,2,3,5) )

  P_qspt    <- apply( X = P_ispt, FUN = quantile, probs = c(0.025, 0.5, 0.975),
                      MARGIN = c(2,3,4) )

  if( is.null(traces))
    traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)
  tMin <-1

  if( proj )
    tMin <- tMP - 1

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(4,3,3,3) )


  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs[tMin:nT]),
            y = c(0,max(P_qspt[,s,p,tMin:nT], na.rm = T) ),
            type = "n", axes = F)

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
        rmtext( outer = TRUE, cex = 1.5, txt = stockNames[p],
                font = 2, line = 1)
      box()
      grid()
      polygon(  x = c(yrs, rev(yrs)),
                y = c(P_qspt[1,s,p,], rev(P_qspt[3,s,p,])),
                col = "grey65", border = NA )
      lines( x = yrs, y = P_qspt[2,s,p,], lwd = 3 )

      for( tIdx in traces )
        lines( x = yrs, y = P_ispt[tIdx,s,p,], lwd = .8 )

      abline( v = yrs[tMP], col = "grey30", lty = 3 )
    

      if( mfg[1] == 1 & mfg[2] == 1 & leg )
        legend( x = "topleft", bty = "n",
                legend = c( "Median Ponded Fish", 
                            "Central 95%",
                            "Replicate Traces"),
                col = c(  "black", "grey65", "black"),
                pch = c(NA,15, NA ),
                lty = c(1, NA, 1),
                lwd = c(3, NA, .8 ) )
    }
  }
  mtext( side = 2, outer = TRUE, text = "Ponded Fish (kt)",
          line = 1.5, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)

  mtext(  outer = TRUE, side = 1, adj = .8, line = 3, cex = .6,
            text = stamp, col = "grey60" )

} # END plotTulipPondedFish

# plotTulipMultiPanel()
# Mulit-panel plot showing simulation envelopes for Spawning 
# Stock Biomass (end of year), total start-of-year biomass,
# catch, and natural mortality.
plotTulipMultiPanel <- function(  obj = blob,
                                  nTrace = 3,
                                  traces = NULL,
                                  traceSeed = 123,
                                  Bmsy = 41.26,
                                  B0   = 95.13 )
{
  goodReps <- obj$goodReps

  D_ispt <- SB_ispt   <- obj$om$SB_ispt[goodReps,,,,drop = FALSE]
  C_ispt              <- obj$om$C_ispt[goodReps,,,,drop = FALSE]

  SB_ispt[SB_ispt == 0] <- NA

  tMP       <- obj$om$tMP
  nS        <- obj$om$nS
  nP        <- obj$om$nP 
  nT        <- obj$om$nT
  nReps     <- sum(goodReps)
  pT        <- obj$ctlList$opMod$pT
  tInit_sp  <- obj$om$tInit_sp

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear


  yrs <- seq( from = fYear, by = 1, length.out = nT)

  # Get reference points, calc depletion
  B0_isp     <- array(NA, dim = c(nReps,nS,nP))

  for(i in 1:nReps)
  {
    if(is.null(B0))
      B0_isp[i,,]     <- obj$rp[[i]]$B0_sp[,,1]
    else B0_isp[i,,]  <- B0

    D_ispt[i,,,]    <- SB_ispt[i,,,]/B0_isp[i,,]
    C_ispt[i,,,]    <- C_ispt[i,,,]/B0_isp[i,,]
  }

  Dmsy <- Bmsy/B0

  # Axis label for SSB depletion
  sbDepAxisLab <- expression(B[t]/B[0])

  # Now take quantiles
  SB_qspt <- apply( X = SB_ispt, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  D_qspt  <- apply( X = D_ispt, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )



  C_qspt <- apply( X = C_ispt, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  SB_qspt[SB_qspt == 0] <- NA 

  if( is.null(traces))
  {
    set.seed(traceSeed)
    traces <- sample( 1:sum(goodReps), size = min(nTrace,sum(goodReps))  )
  }

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")
  tPlotIdx <- tInit_sp[1,1]:nT

  par(  mfcol = c(2,1), 
        mar = c(.5,1.5,.5,1.5),
        oma = c(4,3,3,3) )

  # SSB depletion
  s <- 1
  p <- 1

  plot( x = range(yrs),
            y = c(0,max(D_qspt[,s,p,], na.rm = T) ),
            type = "n", axes = F)
    axis(side = 2, las = 1)
    mtext(side = 2, text = sbDepAxisLab, line = 3)
    box()
    grid()

    polygon(  x = c(yrs[tPlotIdx], rev(yrs[tPlotIdx])),
              y = c(D_qspt[1,s,p,tPlotIdx], rev(D_qspt[3,s,p,tPlotIdx])),
              col = "grey75", border = NA )

    rect( xleft = yrs[tPlotIdx] - .3, xright = yrs[tPlotIdx] + .3,
          ybottom = 0, ytop = C_qspt[2,s,p,tPlotIdx], col = "grey40",
          border = NA )


    abline( h = Dmsy, lty = 3, col = "darkgreen", lwd = 2)
    abline( h = 0.3, lty = 3, col = "red", lwd = 2)
    

    segments( x0 = yrs[tMP:nT], x1 = yrs[tMP:nT],
              y0 = C_qspt[1,s,p,tMP:nT], y1 = C_qspt[3,s,p,tMP:nT],
              col = "black" )

    for( tIdx in traces )
      lines( x = yrs[tPlotIdx], y = D_ispt[tIdx,s,p,tPlotIdx], lwd = .8 )

    lines( x = yrs[tPlotIdx], y = D_qspt[2,s,p,tPlotIdx], lwd = 3 )
    
    abline( v = yrs[tMP], col = "black", lty = 3 )
  
  plot( x = range(yrs),
            y = c(0,max(SB_qspt[,s,p,], na.rm = T) ),
            type = "n", axes = F)
    axis(side = 2, las = 1 )
    mtext(side = 2, line = 3, text = "Spawning Biomass (kt)")
    box()
    grid()

    polygon(  x = c(yrs[tPlotIdx], rev(yrs[tPlotIdx])),
              y = c(SB_qspt[1,s,p,tPlotIdx], rev(SB_qspt[3,s,p,tPlotIdx])),
              col = "grey75", border = NA )

    for( tIdx in traces )
      lines( x = yrs[tPlotIdx], y = SB_ispt[tIdx,s,p,tPlotIdx], lwd = .8 )

    abline( h = Bmsy, lty = 3, col = "darkgreen", lwd = 2)
    abline( h = 0.3*B0, lty = 3, col = "red", lwd = 2)

    lines( x = yrs[tPlotIdx], y = SB_qspt[2,s,p,tPlotIdx], lwd = 3 )
    
    abline( v = yrs[tMP], col = "black", lty = 3 )
  
  axis(side = 1)
  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)

}

# Biomass envelopes - with catch
# if chosen
plotTulipBt <- function(  obj = blob, nTrace = 3,
                          traces = NULL,
                          dep = FALSE,
                          ref = "B0",
                          var = "SB_ispt",
                          Ct  = FALSE,
                          leg = TRUE,
                          proj = FALSE,
                          tMin = NULL,
                          traceSeed = 123 )
{
  goodReps <- obj$goodReps

  SB_ispt   <- obj$om[[var]][goodReps,,,,drop = FALSE]
  C_ispt    <- obj$om$C_ispt[goodReps,,,,drop = FALSE]

  SB_ispt[SB_ispt == 0] <- NA

  
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nReps   <- sum(goodReps)
  pT      <- obj$ctlList$opMod$pT

  # Get reference points
  B0_isp     <- array(NA, dim = c(nReps,nS,nP))
  

  tInit_sp  <- obj$om$tInit_sp

  for(i in 1:nReps)
    B0_isp[i,,]     <- obj$rp[[which(goodReps)[i]]]$B0_sp[,,1]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(is.null(tMin))
    tMin <-1

  if( proj )
    tMin <- tMP - 1

  if( dep )
  {

    if( ref == "B0" )
    {

      for( s in 1:nS )
        for( p in 1:nP )
        {
          for(t in 1:nT)
          {
            SB_ispt[,s,p,t] <- SB_ispt[,s,p,t] / B0_isp[,s,p]
            C_ispt[,s,p,t]  <- C_ispt[,s,p,t] / B0_isp[,s,p]
          }
        }
      # LCP_p     <- LCP_p / B0_sp
      B0_isp     <- B0_isp / B0_isp

    }

    if( ref == "Bmsy" )
    {
      for( s in 1:nS )
        for( p in 1:nP )
        {
          SB_ispt[,s,p,] <- SB_ispt[,s,p,] / Bmsy_sp[s,p]
          C_ispt[,s,p,]  <- C_ispt[,s,p,] / Bmsy_sp[s,p]
        }

      # LCP_p     <- LCP_p / Bmsy_sp[1,]
      B0_sp     <- B0_sp / Bmsy_sp
      BmsySS_sp <- BmsySS_sp / Bmsy_sp
      
    }
  }

  if( !dep )
    yAxisLab <- "Biomass (kt)"

  if( dep )
  {
    if( ref == "Bmsy" )
      yAxisLab <- expression(B[t]/B[MSY])

    if( ref == "B0" )
      yAxisLab <- expression(B[t]/B[0])
  }

  # Now take quantiles
  SB_qspt <- apply( X = SB_ispt, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  C_qspt <- apply( X = C_ispt, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  # Estimate SB0 estimate for years 100-500
  # if(pT>100)
  # {
  #     t1 <- tMP + 50
  #     t2 <- tMP+pT - 1
  #     simB0_sp <- apply( X = SB_ispt[,,,t1:t2,drop=F], FUN = mean,
  #                       MARGIN = c(2,3), na.rm = T )
  # }  


  SB_qspt[SB_qspt == 0] <- NA 

  if( is.null(traces))
  {
    set.seed(traceSeed)
    traces <- sample( 1:sum(goodReps), size = min(nTrace,sum(goodReps))  )
  }

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(4,3,3,3) )


  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs[tMin:nT]),
            y = c(0,max(SB_qspt[,s,p,tMin:nT], na.rm = T) ),
            type = "n", axes = F)

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      # if( mfg[2] == mfg[4] )
      #   rmtext( outer = TRUE, cex = 1.5, txt = stockNames[p],
      #           font = 2, line = 10)
      box()
      grid()
      tPlotIdx <- tInit_sp[s,p]:nT
      polygon(  x = c(yrs[tPlotIdx], rev(yrs[tPlotIdx])),
                y = c(SB_qspt[1,s,p,tPlotIdx], rev(SB_qspt[3,s,p,tPlotIdx])),
                col = "grey65", border = NA )
      lines( x = yrs[tPlotIdx], y = SB_qspt[2,s,p,tPlotIdx], lwd = 3 )

      for( tIdx in traces )
        lines( x = yrs[tPlotIdx], y = SB_ispt[tIdx,s,p,tPlotIdx], lwd = .8 )

      if( Ct )
      {
        rect( xleft = yrs[tPlotIdx] - .3, xright = yrs[tPlotIdx] + .3,
              ybottom = 0, ytop = C_qspt[2,s,p,tPlotIdx], col = "grey65",
              border = NA )
        segments( x0 = yrs[tMP:nT], x1 = yrs[tMP:nT],
                  y0 = C_qspt[1,s,p,tMP:nT], y1 = C_qspt[3,s,p,tMP:nT],
                  col = "black" )
      }

      abline( v = yrs[tMP], col = "grey30", lty = 3 )
      # abline( h = B0_sp[s,p], lty = 3, col = "grey50", lwd = 2  )
      # abline( h = 0.3*B0_sp[s,p], lty = 2, col = "red", lwd = 3)
      # abline( h = LCP_p[p], lty = 2, col = "salmon", lwd = 3)
      # abline( h = BmsySS_sp[s,p], lty = 3, col = "darkgreen", lwd = 2)

      if( mfg[1] == 1 & mfg[2] == 1 & leg )
        legend( x = "topright", bty = "n",
                legend = c( "Median Spawning Biomass", 
                            "Central 95%",
                            "Replicate Traces"),
                            # "Unfished Biomass",
                            # expression(paste("LRP = 0.3",B[0]))),
                            # "Median Projection SB" ),
                            # expression(B[MSY,MS])),
                col = c(  "black", "grey65", "black",
                          "grey50", "red","salmon" ),
                pch = c(NA,15, NA ),# NA, NA,NA,NA),
                lty = c(1, NA, 1), # 3, 2, 2),
                lwd = c(3, NA, .8) )#, 2, 3, 3 ) )


    }
  }
  mtext( side = 2, outer = TRUE, text = yAxisLab,
          line = 1.5, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)

  mtext(  outer = TRUE, side = 1, adj = .8, line = 3, cex = .6,
            text = stamp, col = "grey60" )
}




# Catch envelopes
plotTulipCt <- function(  obj = blob, nTrace = 3,
                          ref = "B0",
                          leg = TRUE,
                          tMin = NULL,
                          plotFleets=FALSE )
{
  goodReps <- obj$goodReps

  C_ispt    <- obj$om$C_ispt[goodReps,,,,drop = FALSE]
  C_ispft   <- obj$om$C_ispft[goodReps,,,,,drop = FALSE]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(C_ispt)[1]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  # fishing fleets
  fishG  <- obj$ctlList$opMod$commGears
  fleets <- obj$ctlList$opMod$fleets
  fColrs <- brewer.pal(nF, 'Dark2')

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(is.null(tMin))
    tMin <-1
 
  C_qspt <- apply( X = C_ispt, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  C_qspft <- apply( X = C_ispft, FUN = quantile,
                    MARGIN = c(2,3,4,5), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )


  traces <- sample( 1:dim(C_ispt)[1], size = min(nTrace,nReps)  )

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(4,3,3,3) )

  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs[tMin:nT]),
            y = c(0,max(C_qspt[,s,p,tMin:nT], na.rm = T) ),
            type = "n", axes = F)

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+0.1, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      polygon(  x = c(yrs, rev(yrs)),
                y = c(C_qspt[1,s,p,], rev(C_qspt[3,s,p,])),
                col = "grey65", border = NA )
      lines( x = yrs, y = C_qspt[2,s,p,], lwd = 3 )


      for( tIdx in traces )
        lines( x = yrs, y = C_ispt[tIdx,s,p,], lwd = .8 )

      abline( v = yrs[tMP], col = "grey30", lty = 3 )

      if( mfg[1] == 1 & mfg[2] == 1 & leg )
        legend( x = "topleft", bty = "n",
                legend = c( "Median Catch (kt)", 
                            "Central 95%",
                            "Replicate Traces"),
                col = c(  "black", "grey65", "black"),
                pch = c(NA,15, NA),
                lty = c(1, NA, 1),
                lwd = c(3, NA, .8) )


      # plot individual lines for each fleet
      if(plotFleets)
      {
        for( f in fishG )
        {
          
          if(sum(C_qspft[2,s,p,f,],na.rm=T)==0)
            next()

          lines( x = yrs, y = C_qspft[2,s,p,f,], col=fColrs[f], lwd = 1.5)
          # for( tIdx in traces )
          #   lines( x = yrs, y = F_ispft[tIdx,s,p,f,], lwd = .8 )
        }
        
        if( mfg[1] == 1 & mfg[2] == 1 & leg & !is.null(fleets) )
        legend('bottomleft', bty='n', cex=0.8,
                legend=c(fleets[fishG]),
                lwd=c(rep(1,length(fishG))),
                lty=c(rep(1,length(fishG))),
                col=c(fColrs[fishG])
                )
      }
    }
  }
  mtext( side = 2, outer = TRUE, text = 'Catch (1000 t)',
          line = 1.5, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)

  mtext(  outer = TRUE, side = 1, adj = .8, line = 3, cex = .6,
            text = stamp, col = "grey60" )
} # END plotTulipCt()

# Catch envelopes
plotTulipPredC_ft <- function( obj = blob, nTrace = 3,
                            leg = TRUE,
                            tMin = NULL )
{
  goodReps <- obj$goodReps

  C_ispft   <- obj$om$C_ispft[goodReps,,,,,drop = FALSE]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(C_ispft)[1]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  # fishing fleets
  fishG  <- obj$ctlList$opMod$commGears
  predG  <- obj$ctlList$opMod$predGears

  nPreds <- length(predG)
  
  fleets <- obj$ctlList$opMod$fleets
  fColrs <- brewer.pal(nPreds, 'Dark2')

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(is.null(tMin))
    tMin <-1
 
  C_qspft <- apply( X = C_ispft, FUN = quantile,
                    MARGIN = c(2,3,4,5), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")



  par(  mfcol = c(nPreds,nP), 
        mar = c(.1,1.5,.1,1.5),
        oma = c(4,3,3,3) )

  for( p in 1:nP )
  {
    for( fIdx in 1:nPreds )
    {
      f <- predG[fIdx]
      plot( x = range(yrs[tMin:nT]),
            y = c(0,max(C_qspft[,1,p,f,tMin:nT], na.rm = T) ),
            type = "n", axes = F)

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      
      if(mfg[2] == 1)
        axis( side = 2, las = 1 )
      box()
      grid()
      polygon(  x = c(yrs, rev(yrs)),
                y = c(C_qspft[1,1,p,f,], rev(C_qspft[3,1,p,f,])),
                col = fColrs[fIdx], border = NA )
      lines( x = yrs, y = C_qspft[2,1,p,f,], lwd = 3 )


      for( tIdx in traces )
        lines( x = yrs, y = C_ispft[tIdx,1,p,f,], lwd = .8 )

      abline( v = yrs[tMP], col = "grey30", lty = 3 )

      legend(x = "topleft", legend = fleets[f], bty = "n")


    }
    
  }

  mtext( side = 2, outer = TRUE, text = 'Consumption (kt)',
          line = 1.5, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)

  mtext(  outer = TRUE, side = 1, adj = .8, line = 3, cex = .6,
            text = stamp, col = "grey60" )
} # END plotTulipCt()

# Catch envelopes
plotTulipPredVulnB_ft <- function(  obj = blob, nTrace = 3,
                                    leg = TRUE,
                                    tMin = NULL )
{
  goodReps <- obj$goodReps

  vB_ispft   <- obj$om$vB_ispft[goodReps,,,,,drop = FALSE]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(vB_ispft)[1]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  # fishing fleets
  fishG  <- obj$ctlList$opMod$commGears
  predG  <- obj$ctlList$opMod$predGears

  nPreds <- length(predG)
  
  fleets <- obj$ctlList$opMod$fleets
  fColrs <- brewer.pal(nPreds, 'Dark2')

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(is.null(tMin))
    tMin <-1
 
  vB_qspft <- apply( X = vB_ispft, FUN = quantile,
                    MARGIN = c(2,3,4,5), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")



  par(  mfcol = c(nPreds,nP), 
        mar = c(.1,1.5,.1,1.5),
        oma = c(4,3,3,3) )

  for( p in 1:nP )
  {
    for( fIdx in 1:nPreds )
    {
      f <- predG[fIdx]
      plot( x = range(yrs[tMin:nT]),
            y = c(0,max(vB_qspft[,1,p,f,tMin:nT], na.rm = T) ),
            type = "n", axes = F)

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      
      if(mfg[2] == 1)
        axis( side = 2, las = 1 )
      box()
      grid()
      polygon(  x = c(yrs, rev(yrs)),
                y = c(vB_qspft[1,1,p,f,], rev(vB_qspft[3,1,p,f,])),
                col = fColrs[fIdx], border = NA )
      lines( x = yrs, y = vB_qspft[2,1,p,f,], lwd = 3 )


      for( tIdx in traces )
        lines( x = yrs, y = vB_ispft[tIdx,1,p,f,], lwd = .8 )

      abline( v = yrs[tMP], col = "grey30", lty = 3 )

      legend(x = "topleft", legend = fleets[f], bty = "n")


    }
    
  }

  mtext( side = 2, outer = TRUE, text = 'Vulnerable Herring Biomass (kt)',
          line = 1.5, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)

  mtext(  outer = TRUE, side = 1, adj = .8, line = 3, cex = .6,
            text = stamp, col = "grey60" )
} # END plotTulipCt()

# Catch envelopes
plotTulipPredF_ft <- function(  obj = blob, nTrace = 3,
                                leg = TRUE,
                                tMin = NULL )
{
  goodReps <- obj$goodReps

  F_ispft   <- obj$om$F_ispft[goodReps,,,,,drop = FALSE]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(F_ispft)[1]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  # fishing fleets
  fishG  <- obj$ctlList$opMod$commGears
  predG  <- obj$ctlList$opMod$predGears
  nPreds <- length(predG)
  
  fleets <- obj$ctlList$opMod$fleets
  fColrs <- brewer.pal(nPreds, 'Dark2')

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(is.null(tMin))
    tMin <-1
 
  F_qspft <- apply( X = F_ispft, FUN = quantile,
                    MARGIN = c(2,3,4,5), probs = c(0.25, 0.5, 0.75),
                    na.rm = T )

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")


  par(  mfcol = c(nPreds,nP), 
        mar = c(.1,1.5,.1,1.5),
        oma = c(4,3,3,3) )

  for( p in 1:nP )
  {
    for( fIdx in 1:nPreds )
    {
      f <- predG[fIdx]
      plot( x = range(yrs[tMin:nT]),
            y = c(0,max(F_qspft[,1,p,f,tMin:nT], na.rm = T) ),
            type = "n", axes = F)

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      
      if(mfg[2] == 1)
        axis( side = 2, las = 1 )
      box()
      grid()
      polygon(  x = c(yrs, rev(yrs)),
                y = c(F_qspft[1,1,p,f,], rev(F_qspft[3,1,p,f,])),
                col = fColrs[fIdx], border = NA )
      lines( x = yrs, y = F_qspft[2,1,p,f,], lwd = 3 )


      for( tIdx in traces )
        lines( x = yrs, y = F_ispft[tIdx,1,p,f,], lwd = .8 )

      abline( v = yrs[tMP], col = "grey30", lty = 3 )

      legend(x = "topleft", legend = fleets[f], bty = "n")


    }
    
  }

  mtext( side = 2, outer = TRUE, text = 'Predation Mortality (/yr)',
          line = 1.5, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)

  mtext(  outer = TRUE, side = 1, adj = .8, line = 3, cex = .6,
            text = stamp, col = "grey60" )
} # END plotTulipCt()


# Catch envelopes
plotTulipMt <- function(  obj = blob, nTrace = 3,
                          leg = TRUE, tMin = NULL,
                          traceSeed = 123)
{
  goodReps <- obj$goodReps

  M_iaxspt <- obj$om$M_iaxspt[goodReps,,,,,,drop = FALSE]


  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(M_iaxspt)[1]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(is.null(tMin))
    tMin <-1
 
  M_qspt <- apply( X = M_iaxspt[,5,1,,,,drop=FALSE], FUN = quantile,
                    MARGIN = c(4,5,6), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  set.seed(traceSeed)
  traces <- sample( which(goodReps), size = min(nTrace,nReps)  )

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(4,3,3,3) )

  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs[tMin:nT]),
            y = c(0,max(M_qspt[,s,p,tMin:nT], na.rm = T) ),
            type = "n", axes = F)

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+1, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      polygon(  x = c(yrs, rev(yrs)),
                y = c(M_qspt[1,s,p,], rev(M_qspt[3,s,p,])),
                col = "grey65", border = NA )
      lines( x = yrs, y = M_qspt[2,s,p,], lwd = 3 )


      for( tIdx in traces )
        lines( x = yrs, y = M_iaxspt[tIdx,5,1,s,p,], lwd = .8 )

      abline( v = yrs[tMP], col = "grey30", lty = 3 )

      
        legend( x = "topleft", bty = "n",
                legend = c( "Median M", 
                            "Central 95%",
                            "Replicate Traces"),
                col = c(  "black", "grey65", "black"),
                pch = c(NA,15, NA),
                lty = c(1, NA, 1),
                lwd = c(3, NA, .8) )
      
    }
  }
  mtext( side = 2, outer = TRUE, text = 'Natural mortality',
          line = 1.5, font = 2)

  mtext( side = 1, outer = TRUE, text = "",
          line = 2, font = 2)

  mtext(  outer = TRUE, side = 1, adj = .8, line = 3, cex = .6,
            text = stamp, col = "grey60" )
} # END plotTulipMt()

# plotTulipEffort_p()
# Effort over time gridded
# by stock area and fleet (right now predators) - envelopes
plotTulipEffort_p <- function(  obj = blob, 
                                nTrace = 3,
                                pIdx = 1,
                                fIdx = 7:13 )
{
  goodReps <- obj$goodReps

  E_ipft <- obj$om$E_ipft[goodReps,pIdx,,,drop = FALSE ]

  E_ipft[E_ipft == 0] <- NA

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(E_ipft)[1]

  E_qpft <- apply(  X = E_ipft, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  E_qpft[E_qpft == 0] <- NA

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  fleetCols <- RColorBrewer::brewer.pal( length(fIdx), "Dark2" )

  par(  mfcol = c(length(fIdx),nP), 
        mar = c(.1,1.5,.1,1.5),
        oma = c(3,4,3,3) )

  for(p in pIdx)
    for( i in 1:length(fIdx) )
    {
      f <- fIdx[i]
      plot( x = range(yrs),
            y = c(0,max(E_qpft[,,f,],na.rm = T) ),
            type = "n", axes = F )
        mfg <- par("mfg")
        # Plot axes and facet labels
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        axis( side = 2, las = 1 )
        if( mfg[2] == mfg[4] )
        {
          rmtext( txt =  fleetNames[f], font = 2, line = .05,
                  cex = 1.5, outer = TRUE )
        }
        box()
        grid()

        polygon(  x = c(yrs,rev(yrs)), y = c(E_qpft[1,p,f,],rev(E_qpft[3,p,f,])), 
                  col = scales::alpha(fleetCols[i], alpha = .3), border = NA )
        
        for( tIdx in traces )
          lines( x = yrs, y = E_ipft[tIdx,p,f,], col = fleetCols[i], lwd = .8 )
        
        lines( x = yrs, y = E_qpft[2,p,f,], col = fleetCols[i], lwd = 3)
        abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }
  mtext( side = 2, outer = TRUE, text = "Predator Effort (biomass/abundance)",
          line = 2, font = 2)
} # END plotTulipEffort_p()


# plotTulipCatch_pft()
# Catch over time gridded
# by stock area and fleet (right now predators) - envelopes
plotTulipCatch_pft <- function( obj = blob, 
                                nTrace = 3,
                                pIdx = 1,
                                fIdx = 7:13 )
{
  goodReps <- obj$goodReps

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(obj$om$C_ispft)[1]


  C_ipft <- array(0,dim = c(length(goodReps),nP,nF,nT) )

  C_ipft[,pIdx,fIdx,] <- obj$om$C_ispft[goodReps,1,pIdx,fIdx,]
  

  C_ipft[C_ipft == 0] <- NA

  

  C_qpft <- apply(  X = C_ipft, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  C_qpft[C_qpft == 0] <- NA

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  fleetCols <- RColorBrewer::brewer.pal( length(fIdx), "Dark2" )

  par(  mfcol = c(length(fIdx),nP), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,4,3,3) )

  for(p in pIdx)
    for( i in 1:length(fIdx) )
    {
      f <- fIdx[i]
      plot( x = range(yrs),
            y = c(0,max(C_qpft[,,f,],na.rm = T) ),
            type = "n", axes = F )
        mfg <- par("mfg")
        # Plot axes and facet labels
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        axis( side = 2, las = 1 )
        if( mfg[2] == mfg[4] )
        {
          corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          par(xpd = TRUE) #Draw outside plot area
          text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
                font = 2, cex = 1.5 )
          par(xpd = FALSE)
        }
        box()
        grid()

        polygon(  x = c(yrs,rev(yrs)), y = c(C_qpft[1,p,f,],rev(C_qpft[3,p,f,])), 
                  col = scales::alpha(fleetCols[i], alpha = .3), border = NA )
        
        for( tIdx in traces )
          lines( x = yrs, y = C_ipft[tIdx,p,f,], col = fleetCols[i], lwd = .8 )
        
        lines( x = yrs, y = C_qpft[2,p,f,], col = fleetCols[i], lwd = 3)
        abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }
  mtext( side = 2, outer = TRUE, text = "Trawl Effort (fishing hours?)",
          line = 2, font = 2)
} # END plotTulipCatch_pft()

# plotTulipPondedBio_pft()
# Catch over time gridded
# by stock area and fleet (right now predators) - envelopes
plotTulipPondedBio_pft <- function( obj = blob, 
                                    nTrace = 3,
                                    pIdx = 1 )
{
  goodReps <- obj$goodReps

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(obj$om$C_ispft)[1]

  # Get SOK fleets
  browser()
  sokFleets <- c(6,13)

  P_ipft <- array(0,dim = c(length(goodReps),nP,nF,nT) )

  P_ipft[,pIdx,sokFleets,] <- obj$om$P_ispft[goodReps,1,pIdx,sokFleets,]
  

  P_ipft[P_ipft == 0] <- NA

  

  P_qpft <- apply(  X = P_ipft, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  P_qpft[P_qpft == 0] <- NA

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  fleetCols <- RColorBrewer::brewer.pal( length(sokFleets), "Dark2" )

  par(  mfcol = c(length(sokFleets),nP), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,4,3,3) )

  for(p in pIdx)
    for( i in 1:length(sokFleets) )
    {
      f <- sokFleets[i]
      plot( x = range(yrs),
            y = c(0,max(P_qpft[,,f,],na.rm = T) ),
            type = "n", axes = F )
        mfg <- par("mfg")
        # Plot axes and facet labels
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        axis( side = 2, las = 1 )
        if( mfg[2] == mfg[4] )
        {
          corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          par(xpd = TRUE) #Draw outside plot area
          text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
                font = 2, cex = 1.5 )
          par(xpd = FALSE)
        }
        box()
        grid()

        polygon(  x = c(yrs,rev(yrs)), y = c(P_qpft[1,p,f,],rev(P_qpft[3,p,f,])), 
                  col = scales::alpha(fleetCols[i], alpha = .3), border = NA )
        
        for( tIdx in traces )
          lines( x = yrs, y = P_ipft[tIdx,p,f,], col = fleetCols[i], lwd = .8 )
        
        lines( x = yrs, y = P_qpft[2,p,f,], col = fleetCols[i], lwd = 3)
        abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }
  mtext( side = 2, outer = TRUE, text = "Trawl Effort (fishing hours?)",
          line = 2, font = 2)
} # END plotTulipPondedBio_pft()



# plotTulipSOKEffort_p()
# Effort over time gridded
# by stock area - envelopes
plotTulipSOKEffort_p <- function( obj = blob, 
                                  nTrace = 3,
                                  proj = TRUE,
                                  traces = NULL )
{
  goodReps <- obj$goodReps

  E_ipft <- obj$mp$hcr$sokEff_ispft[goodReps,1,,,]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(E_ipft)[1]

  Eagg_ift <- apply(X = E_ipft, FUN = sum, MARGIN = c(1,3,4), na.rm = T)

  Etot_ipft <- array(0, dim = dim(E_ipft) + c(0,1,0,0) )
  Etot_ipft[,1:nP,,] <- E_ipft
  Etot_ipft[,nP+1,,] <- Eagg_ift

  
  E_ipft <- Etot_ipft


  # Pull scenario and MP labels
  scenarioName <- obj$ctlList$ctl$scenarioName
  mpName <- obj$ctlList$ctl$mpName

  stamp <- paste0(scenarioName,":",mpName)

  E_qpft <- apply(  X = E_ipft, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = TRUE )

  if( is.null(traces))
    traces <- sample( 1:nReps, size = min(nTrace,nReps)  )


  speciesNames  <- obj$om$speciesNames
  stockNames    <- c("C/S","JP/S","Lou","Aggregate")
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(proj)
    tdx <- (tMP - 1):nT
  else tdx <- 1:nT

  fleetCols <- RColorBrewer::brewer.pal( nF, "Dark2" )

  par(  mfcol = c(nP+1,1), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,4,3,3) )

  for( p in 1:(nP+1) )
  {
    plot( x = range(yrs[tdx]),
          y = c(0,max(E_qpft[,p,,tdx],na.rm = T) ),
          type = "n", axes = F )
      mfg <- par("mfg")
      # Plot axes and facet labels
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      axis( side = 2, las = 1 )
      box()
      grid()
      for( f in 6:7 )
      {
        # Skip if the fleet doesn't fish
        if( all(E_qpft[,,f,] == 0) )
          next

        polygon( x = c(yrs,rev(yrs)), y = c(E_qpft[1,p,f,],rev(E_qpft[3,p,f,])), 
                  col = scales::alpha(fleetCols[f], alpha = .3), border = NA )
          for( tIdx in traces )
            lines( x = yrs, y = E_ipft[tIdx,p,f,], col = fleetCols[f], lwd = .8 )

        lines( x = yrs, y = E_qpft[2,p,f,], col = fleetCols[f], lwd = 3)
        
      }
      if( mfg[2] == mfg[4] )
      {
        rmtext( line = 0.25, txt = stockNames[p], font = 2, outer = FALSE, cex = 1.5 )
      }
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
  }
  mtext( side = 2, outer = TRUE, text = "Active SOK licenses",
          line = 2, font = 2)
  mtext( side = 1, text = stamp, col = "grey60", outer = TRUE, 
          line = 1, adj = .75, cex = .5 )
} # END plotTulipSOKEffort_p()

# rmtext()
# Refactored procedure to plot right hand inner
# margin mtext with the bottom towards the middle
# of the plot
rmtext <- function( line = 1, 
                    txt = "Sample", 
                    font = 1,
                    cex = 1,
                    outer = FALSE,
                    yadj = .5)
{
  corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
  if( outer )
    par(xpd = NA) #Draw outside the figure region
  if( !outer )
    par(xpd = TRUE)

  xRange <- corners[2] - corners[1]

  text( x = corners[2]+line, 
        y = yadj * sum(corners[3:4]), 
        labels = txt, srt = 270,
        font = font, cex = cex )
  par(xpd = FALSE)
} # END rmtext()


# plotEffort_p()
# Effort over time gridded
# by stock area
plotEffort_p <- function( obj = blob,
                          iRep = 1 )
{
  E_pft <- obj$om$E_ipft[iRep,,,]

  E_pft[E_pft == 0] <- NA

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  fleetCols <- RColorBrewer::brewer.pal( nF, "Dark2" )

  par(  mfcol = c(nP,1), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,4,3,3) )

  for( p in 1:nP )
  {
    plot( x = range(yrs),
          y = c(0,max(E_pft[p,,],na.rm = T) ),
          type = "n", axes = F )
      mfg <- par("mfg")
      # Plot axes and facet labels
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      for( f in 1:nF )
        lines( x = yrs, y = E_pft[p,f,], col = fleetCols[f], lwd = 3 )
        
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
  }
  mtext( side = 2, outer = TRUE, text = "Trawl Effort (fishing hours?)",
          line = 2, font = 2)
} # END plotEffort_p()


# plotCtTACt_sp()
# Comparison plot of TAC and realized catch for each 
# species/stock
plotCtTACt_sp <- function(  obj = blob,
                            iRep = 1,
                            fleets = 1:2 )
{
  C_spft   <- obj$om$C_ispft[iRep,,,,]
  TAC_spft <- obj$mp$hcr$TAC_ispft[iRep,,,,]

  C_spft[C_spft == 0] <- NA

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  fleetCols <- RColorBrewer::brewer.pal( nF, "Dark2" )


  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )

  for(s in 1:nS)
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(C_spft[s,p,,], TAC_spft[s,p,,], na.rm = T) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      for( f in fleets )
      {
        rect( xleft = yrs - .3,
              xright = yrs + .3,
              ybottom = 0,
              ytop = TAC_spft[s,p,f,],
              border = NA, col = fleetCols[f] )
        
        lines(  x = yrs[], y = C_spft[s,p,f,],
                lty = 2, lwd = 2, col = "grey50" )
      }
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }
    mtext( outer = TRUE, side = 2, text = "Catch and TAC (kt)" )
} # END plotCtTACt_sp()




# plotBtCtRt_p()
# Biomass and catch plots
# by stock for the historical
# and projection period
plotBtCtRt_p <- function( obj = blob, iRep = 1, sIdx=1)
{
  SB_ispt  <- obj$om$SB_ispt
  B_ispt   <- obj$om$B_ispt
  C_ispt   <- obj$om$C_ispt
  R_ispt   <- obj$om$R_ispt
  TAC_ispt <- obj$mp$hcr$TAC_ispt

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT


  # stockNames    <- obj$om$stockNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  par(  mfcol = c(3,nP), 
        mar = c(1,1.5,1,1.5), mgp=c(1,0.4,0),
        oma = c(2,3,0,0), tck=-.01 )

    for( p in 1:nP )
    {
      
      # biomass plots
      plot( x = range(yrs),
            y = c(0,max(SB_ispt[iRep,sIdx,p,]) ),
            type = "n", axes = T, las=1, xlab='', ylab='' )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = stockNames[p], font = 2, line = 0 )
      # axis( side = 2, las = 1 )

      box()
      grid()
      lines( x = yrs, y = SB_ispt[iRep,sIdx,p,], col = "red", lwd = 1 )
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )

      # catch plots
      plot( x = range(yrs),
            y = c(0,max(C_ispt[iRep,sIdx,p,],na.rm=T) ),
            type = "n", axes = T, las=1, xlab='', ylab='' )  
      # axis( side = 2, las = 1 )
      box()
      grid()
      lines( x = yrs, y = C_ispt[iRep,sIdx,p,], col = "black", lwd = 1 )
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )

      # recruitment plots
      plot( x = range(yrs),
            y = c(0,max(R_ispt[iRep,sIdx,p,]) ),
            type = "n", axes = T, las=1, xlab='', ylab='' )
      # axis( side = 2, las = 1 )
      box()
      grid()      
      lines( x = yrs, y = R_ispt[iRep,sIdx,p,], col = "orange", lwd = 1 )
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )

    }

  mtext( side =2, outer = TRUE, text = "Stock SB (kt), catch (kt), and age-1 recruitment (1000s)", line = 1.5)

}

# plotBtCtRt_p()
# Biomass and catch plots
# by stock for the historical
# and projection period
plotBtCtRtMt_p <- function( obj = blob, iRep = 1, sIdx=1)
{
  SB_ispt  <- obj$om$endSB_ispt
  B_ispt   <- obj$om$B_ispt
  C_ispt   <- obj$om$C_ispt
  C_ispft  <- obj$om$C_ispft
  R_ispt   <- obj$om$R_ispt
  M_iaxspt <- obj$om$M_iaxspt

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT

  # replace zeros with NAs for Mortality object
  M_iaxspt[M_iaxspt==0] <- NA

  # Get predator fleets
  predGears <- obj$ctlList$opMod$predGears
  predC_ispt <- array(NA,dim = dim(C_ispt))

  if(length(predGears) > 0 )
  {
    C_ispt <- apply( X = C_ispft[,,,-predGears,,drop = FALSE], 
                     FUN = sum, MARGIN = c(1:3,5), na.rm = T )
    predC_ispt <- apply(  X = C_ispft[,,,predGears,,drop = FALSE], 
                          FUN = sum, MARGIN = c(1:3,5), na.rm = T )

    # Need to add "dead" ponded fish here
  }

  # stockNames    <- obj$om$stockNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  par(  mfcol = c(4,nP), 
        mar = c(1,1.5,1,1.5), mgp=c(1,0.4,0),
        oma = c(2,3,0.5,0), tck=-.01 )

    for( p in 1:nP )
    {
      
      # biomass plots
      plot( x = range(yrs),
            y = c(0,max(SB_ispt[iRep,sIdx,p,]) ),
            type = "n", axes = T, las=1, xlab='', ylab='' )
        mfg <- par("mfg")
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        if( mfg[1] == 1 )
          mtext( side = 3, text = stockNames[p], font = 2, line = 0 )
        # axis( side = 2, las = 1 )
        if(p==1)
          mtext( side = 2, text = 'Spawning Biomass (kt)', line=2.5, cex=0.8 )

        box()
        grid()
        lines( x = yrs, y = SB_ispt[iRep,sIdx,p,], col = "red", lwd = 2 )
        abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 1 )
      # legend('topright',bty='n', cex=0.8,
      #        legend='SBt',lty=1, col='red')

      # catch plots
      plot( x = range(yrs),
            y = c(0,max(C_ispt[iRep,sIdx,p,],
                        predC_ispt[iRep,sIdx,p,],na.rm=T) ),
            type = "n", axes = T, las=1, xlab='', ylab='' )  
      # axis( side = 2, las = 1 )
        if(p==1)
          mtext( side = 2, text = 'Catch (kt)', line=2.5, cex=0.8 )
        box()
        grid()
        lines( x = yrs, y = C_ispt[iRep,sIdx,p,], col = "black", lwd = 1 )
        abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 1 )
        if(length(predGears) > 0)
          lines(  x = yrs, 
                  y = predC_ispt[iRep,sIdx,p,], col = "salmon", lwd = 2)
      # legend('topright',bty='n', cex=0.8,
      #        legend='Ct',lty=1, col='black')

      # recruitment plots
      plot( x = range(yrs),
            y = c(0,max(R_ispt[iRep,sIdx,p,]) ),
            type = "n", axes = T, las=1, xlab='', ylab='' )
      # axis( side = 2, las = 1 )
      if(p==1)
        mtext( side = 2, text = 'Recruits (1000s)', line=2.5, cex=0.8 )

      box()
      grid()      
      lines( x = yrs, y = R_ispt[iRep,sIdx,p,], col = "orange", lwd = 1 )
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 1 )
      legend('topleft',bty='n', cex=0.8,
             legend=c('Age 1 Rt'),lty=1, col='orange')


      # natural mortality plots
      plot( x = range(yrs),
            y = c(0,max(M_iaxspt[iRep,,1,sIdx,p,]) ),
            type = "n", axes = T, las=1, xlab='', ylab='' )
      # axis( side = 2, las = 1 )
      if(p==1)
        mtext( side = 2, text = 'Natural Mortality', line=2.5, cex=0.8 )

      box()
      grid()

      # age1 and age2 mortality
      M1 <- M_iaxspt[iRep,1,1,sIdx,p,]
      M2 <- M_iaxspt[iRep,2,1,sIdx,p,]
      M1[M1==0] <- NA
      M2[M2==0] <- NA

      lines( x = yrs, y = M1, col = "blue", lwd = 1, lty=3 )
      lines( x = yrs, y = M2, col = "blue", lwd = 1, lty=1 )
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 1 )

      legend('topleft',bty='n', cex=0.8,
              legend=c('age 1 M',"age 2 M"),
              lty=c(3,1), col='blue')

    }

  # mtext( side =2, outer = TRUE, text = "Stock SB (kt), catch (kt), age-1 recruitment (1000s), and natural mortality", line = 1.5)

}


# plotBtCt_sp()
# Biomass and catch plots
# by species/stock for the historical
# and projection period
plotBtCt_sp <- function(  obj = blob,
                          iRep = 1 )
{
  SB_ispt  <- obj$om$SB_ispt
  B_ispt   <- obj$om$B_ispt
  C_ispt   <- obj$om$C_ispt
  TAC_ispt <- obj$mp$hcr$TAC_ispt

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )

  for(s in 1:nS)
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(B_ispt[iRep,s,p,]) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      rect( xleft = yrs - .3, xright = yrs + .3,
            ybottom = 0, ytop = C_ispt[iRep,s,p,], col = "grey70", border = NA )
      segments( x0 = yrs, x1 = yrs,
                y0 = 0, y1 = TAC_ispt[iRep,s,p,], col = "black" )
      lines( x = yrs, y = SB_ispt[iRep,s,p,], col = "red", lwd = 1 )
      lines( x = yrs, y = B_ispt[iRep,s,p,], col = "black", lwd = 1, lty = 2 )
      

      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )

    }


  mtext( side =1, outer = TRUE, text = "Year", line = 2)
  mtext( side =2, outer = TRUE, text = "Stock biomass and catch (kt)", line = 1.5)

}


# plotRetroSB()
# Retrospective plots of AM fits for a given replicate
plotRetroSB <- function(  obj  = blob, 
                          iRep = 1, 
                          vB_f = 5,
                          plotTB = FALSE,
                          plotVB = FALSE,
                          plotData = TRUE  )
{
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nT  <- obj$om$nT
  tMP <- obj$om$tMP
  pT  <- nT - tMP + 1
  nF  <- obj$om$nF

  # Need to adjust this for blended index
  blendIdx      <- obj$ctlList$opMod$blendIdx
  # sum( matAge_asp[,s,p] * vB_axspft[,nX,s,p,f,t])
  # browser()
  # Get biomass arrays
  SB_spt       <- ClimProjDiags::Subset(  x= obj$om$SB_ispt[iRep,,,,drop=FALSE],
                                          along = 2:4, indices = list( 1:nS, 1:nP, 1:nT),
                                          drop = "non-selected" )
  
  VB_spft      <- ClimProjDiags::Subset(  x = obj$om$vB_ispft[iRep,,,,,drop=FALSE],
                                          along = 2:5, indices = list(1:nS,1:nP,vB_f,1:nT),
                                          drop = "non-selected")
  
  totB_spt     <- ClimProjDiags::Subset(  x = obj$om$B_ispt[iRep,,,,drop=F],
                                          along = 2:4,
                                          indices = list(1:nS,1:nP,1:nT),
                                          drop = "non-selected")

  q_spft      <- ClimProjDiags::Subset( x = obj$om$q_ispft[iRep,,,,,drop = FALSE],
                                        along = 2:5,
                                        indices = list(1:nS,1:nP,1:nF,1:(nT)),
                                        drop = "non-selected" )


  I_spft       <- ClimProjDiags::Subset(  x = obj$mp$data$I_ispft[iRep,,,,,drop = FALSE],
                                          along = 2:5,
                                          indices = list(1:nS,1:nP,1:nF,1:(nT)),
                                          drop = "non-selected" )

  retroSB_tspt  <- ClimProjDiags::Subset( x = obj$mp$assess$retroSB_itspt[iRep,,,,,drop = FALSE],
                                          along = c(2:5),
                                          indices = list(1:pT,1:nS,1:nP,1:nT),
                                          drop = "non-selected" )

  retroSB_tspt[retroSB_tspt < 0] <- NA

  # Pull blended index stuff
  rI_spft       <- ClimProjDiags::Subset( x = obj$om$rI_ispft[iRep,,,,,drop = FALSE],
                                          along = 2:5,
                                          indices = list(1:nS,1:nP,1:nF,1:(nT)),
                                          drop = "non-selected" )

  qComb_tspt    <- ClimProjDiags::Subset( x = obj$mp$assess$retroqComb_itspt[iRep,,,,,drop = FALSE],
                                          along = 2:5,
                                          indices = list(1:pT,1:nS,1:nP,1:(nT)),
                                          drop = "non-selected" )

  # browser()
  # Get proportion of TACs for splitting aggregated biomass
  propTAC_spt   <- obj$mp$hcr$propTAC_ispt[iRep,,,]

  I_spft[I_spft < 0] <- NA


  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- obj$ctlList$opMod$stock
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )
  for(s in 1:nS)
    for( p in 1:nP )
    {
      maxSB <- 1.2*max(SB_spt[s,p,], plotTB * totB_spt[s,p,], plotVB*VB_spft[s,p,,], na.rm = T)
      plot( x = range(yrs),
            y = c(0,maxSB ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+1.5, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      lines( x = yrs, y = SB_spt[s,p,], col = "red", lwd = 3 )
      if( plotData )
      {
          for( f in 1:nF )
            points(x = yrs, y = I_spft[s,p,f, ]/q_spft[s,p,f,] )

        # if( blendIdx )
        #   points( x = yrs, y = I_spft[s,p,5,] / omqComb_spt[s,p,])
      }

      if( plotVB )
        for( f in vB_f )
          lines( x = yrs, y = VB_spft[s,p,f,], col = "grey40", lwd = 2, lty = 3 )
      if( plotTB )
        lines( x = yrs, y = totB_spt[s,p,], col = "black", lwd = 2 )
      
      # Plot retro fits
      for( tt in 1:pT )
      {
        # propTAC <- propTAC_spt[s,p,tMP + tt - 1]
        lines( x = yrs, y = retroSB_tspt[tt,s,p,], col = "grey60", lwd = 1 )
      }
  
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }

} # plotRetroSB()

# plotBtIt_p()
# Plot Biomass and index
plotBtIt_p <- function( obj = blob, iRep = 1, f=5,
                      plotVB=FALSE, addCatch=TRUE,
                      parArg=TRUE, YlabOn= TRUE,
                      legdOn=TRUE)
{
  # Get biomass arrays
  SB_ispt     <- obj$om$SB_ispt
  VB_ispft    <- obj$om$vB_ispft
  I_ispt      <- obj$mp$data$I_ispft[,,,f,] 
  # nS+1 and nP+1 for I_ispt represent aggregates

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT

  # Hack to plot empirical assessment method for Herring projections
  if(addCatch)
  {
    C_ispt <- obj$om$C_ispt
    
    for(s in 1:nS)
      for(p in 1:nP)
        I_ispt[,s,p,]<- I_ispt[,s,p,]+ C_ispt[,s,p,]

  }  

  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear
  pT            <- obj$ctlList$opMod$pT

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(parArg)
  par(  mfcol = c(nP,nS), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )
  for(s in 1:nS)
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(VB_ispft[iRep,s,p,f,],
                        SB_ispt[iRep,s,p,],
                        I_ispt[iRep,s,p,], na.rm = T) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )

      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+1.5, y = mean(corners[3:4]), stockNames[p], srt = 270,font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      lines( x = yrs, y = SB_ispt[iRep,s,p,], col = "red", lwd = 2 ) 

      # plot vulnerable biomass
      if(plotVB)
      lines( x = yrs, y = VB_ispft[,s,p,,], col = "grey40", lwd = 2, lty = 3 )

      # plot index
      points( x = yrs, y = I_ispt[iRep,s,p,], 
        bg = "green", pch=21, cex=1 )
      
      abline( v = yrs[tMP], lty = 2, lwd = 0.5 )

      if(s==1 & p==1 & !addCatch & legdOn)
      legend('topright',bty='n',
             legend=c('Spawning Biomass',
                      'Dive Survey Index'),
             lwd=c(2,NA), pch=c(NA,21),
             col=c('red','black'),
             pt.bg=c(NA, 'green'))

      if(s==1 & p==1 & addCatch & legdOn)
      legend('topright',bty='n',
             legend=c('Spawning Biomass',
                      'It + Ct'),
             lwd=c(2,NA), pch=c(NA,21),
             col=c('red','black'),
             pt.bg=c(NA, 'green'))
    }

  if(YlabOn)
  mtext( side =2, outer = TRUE, text = "Spawning biomass kt)", line = 1.5)

} # plotBtIt_p()



# plotRetroSBagg()
# Retrospective plots of AM fits for a given replicate
# with biomasses aggregated to match the scale of the AM
plotRetroSBagg <- function( obj = blob, iRep = 1, Ct = TRUE )
{
  # Get biomass arrays
  # Model dims
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  pT      <- obj$ctlList$opMod$pT

  SB_spt        <- array(NA, dim = c(nS,nP,nT))
  VB_spt        <- array(NA, dim = c(nS,nP,nT))
  retroSB_tspt  <- array(NA, dim = c(pT,nS,nP,nT))


  SB_spt[1:nS,,]        <- obj$om$SB_ispt[iRep,1:nS,,]
  VB_spt[1:nS,,]        <- obj$om$vB_ispft[iRep,1:nS,,2,]
  retroSB_tspt[1:pT,,,] <- obj$mp$assess$retroSB_itspt[iRep,1:pT,,,]

  ctlList <- obj$ctlList

  retroSB_tspt[retroSB_tspt < 0] <- NA


  C_spft   <- array(NA, dim = c(nS,nP,nF,nT))
  TAC_spft <- array(NA, dim = c(nS,nP,nF,nT))

  C_spt    <- array(NA, dim = c(nS,nP,nT))
  TAC_spt  <- array(NA, dim = c(nS,nP,nT))

  C_spft[1:nS,,,]   <- obj$om$C_ispft[iRep,,,,]
  TAC_spft[1:nS,,,] <- obj$mp$hcr$TAC_ispft[iRep,,,,]

  C_spt[1:nS,,]    <- apply(X = C_spft, FUN = sum, MARGIN = c(1,2,4))
  TAC_spt[1:nS,,]  <- apply(X = TAC_spft, FUN = sum, MARGIN = c(1,2,4))

  # Aggregate OM biomasses to match AM biomass
  if( ctlList$mp$data$spatialPooling )
  {
    newSB_spt        <- apply( X = SB_spt, FUN = sum, MARGIN = c(1,3), na.rm = T )
    newVB_spt        <- apply( X = VB_spt, FUN = sum, MARGIN = c(1,3), na.rm = T )

    SB_spt    <- SB_spt[,1,,drop = FALSE]
    SB_spt[,1,] <- newSB_spt
    VB_spt    <- VB_spt[,1,,drop = FALSE]
    VB_spt[,1,] <- newVB_spt

    sumC_spt          <- C_spft[,1,1,,drop = FALSE]
    sumC_spt[,1,1,]   <- apply(X = C_spft, FUN = sum, MARGIN = c(1,4))
    newC_spt          <- array(NA, dim = c(nS,1,nT))
    newC_spt[,1,]     <- sumC_spt[,1,1,]
    C_spt             <- newC_spt

    sumTAC_spt        <- C_spft[,1,,,drop = FALSE]
    sumTAC_spt[,1,1,] <- apply(X = TAC_spft, FUN = sum, MARGIN = c(1,4))
    newTAC_spt        <- array(NA, dim = c(nS,1,nT))
    newTAC_spt[,1,]   <- sumTAC_spt[,1,1,]
    TAC_spt           <- newTAC_spt        

  }

  if( ctlList$mp$data$speciesPooling )
  {
    newSB_spt        <- apply( X = SB_spt, FUN = sum, MARGIN = c(2,3), na.rm = T )
    newVB_spt        <- apply( X = VB_spt, FUN = sum, MARGIN = c(2,3), na.rm = T )

    SB_spt    <- SB_spt[1,,,drop = FALSE]
    SB_spt[1,,] <- newSB_spt
    VB_spt    <- VB_spt[1,,,drop = FALSE]
    VB_spt[1,,] <- newVB_spt

    sumC_spt          <- C_spft[1,,1,,drop = FALSE]
    sumC_spt[1,,1,]   <- apply(X = C_spft, FUN = sum, MARGIN = c(2,4))
    newC_spt          <- array(NA, dim = c(1,nP,nT))
    newC_spt[1,,]     <- sumC_spt[1,,1,]
    C_spt             <- newC_spt
    

    sumTAC_spt        <- TAC_spft[1,,1,,drop = FALSE]
    sumTAC_spt[1,,,]  <- apply(X = TAC_spft, FUN = sum, MARGIN = c(2,4))
    newTAC_spt        <- array(NA, dim = c(1,nP,nT))
    newTAC_spt[1,,]   <- sumTAC_spt[1,,1,]
    TAC_spt           <- newTAC_spt
  }

  SB_spt[SB_spt == 0]     <- NA
  VB_spt[VB_spt == 0]     <- NA

  nSS     <- dim( SB_spt)[1]
  nPP     <- dim( SB_spt)[2]

  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- obj$ctlList$opMod$stock
  fYear         <- obj$ctlList$opMod$fYear
  pT            <- obj$ctlList$opMod$pT

  if( nPP == 1 )
    stockNames <- "Coastwide"

  if( nSS == 1 )
    speciesNames <- "Data Pooled"

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,":",iRep,sep = "")

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  par(  mfcol = c(nPP,nSS), 
        mar = c(1,1.5,1,1.5),
        oma = c(4,3,3,3) )
  for(s in 1:nSS)
    for( p in 1:nPP )
    {
      plot( x = range(yrs),
            y = c(0,1.2*max(VB_spt[s,p,],SB_spt[s,p,],na.rm = T) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+1.5, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      if( Ct )
      {
        # Plot actual catch
        rect( xleft = yrs - .3, xright = yrs + .3, 
              ytop = C_spt[s,p,], ybottom = 0, col = "grey40",
              border = NA )
        # Plot a rectangle of TACs
        rect( xleft = yrs[tMP:nT] - .3, xright = yrs[tMP:nT] + .3, 
              ytop = TAC_spt[s,p,tMP:nT], ybottom = 0, col = NA,
              border = "black" )

      }
      lines( x = yrs, y = SB_spt[s,p,], col = "red", lwd = 3 )
      lines( x = yrs, y = VB_spt[s,p,], col = "grey40", lwd = 2, lty = 3 )
      for( tt in 1:pT )
      {
        lines( x = yrs, y = retroSB_tspt[tt,s,p,], col = "grey60", lwd = 1 )
      }
      
  
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }
  mtext( side = 1, text = "Year", outer = TRUE, line = 2, font = 2)
  mtext( side = 2, text = "Spawning Biomass, TAC, and Catch (kt)", 
          outer = TRUE, line = 2, font = 2)
  mtext(  outer = TRUE, side = 1, adj = .8, line = 3, cex = .6,
            text = stamp, col = "grey60" )

}

# plotEBSBratio()
# Plots of ratio between vulnerable
# biomass and exploitable biomass
# for given depletion levels.
plotEBSBratio <- function( obj = blob )
{
  # Pull model dims
  nS <- obj$om$nS
  nP <- obj$om$nP
  nF <- obj$om$nF

  fleetNames <- obj$om$fleetNames

  # Pull reference points
  rp <- obj$rp[[1]]

  nu_spfk <- array(NA, dim = c(nS+1,nP+1,nF,3))
  P1_spf  <- array(NA, dim = c(nS+1,nP+1,nF))
  P2_spf  <- array(NA, dim = c(nS+1,nP+1,nF))

  # Fill!
  # nu pars
  nu_spfk[1:nS,1:nP,1:nF,] <- rp$EBSBpars$stockSpec$nu_spfk
  nu_spfk[nS+1,1:nP,1:nF,] <- rp$EBSBpars$coastWide$nu_spfk
  nu_spfk[1:nS,1+nP,1:nF,] <- rp$EBSBpars$dataPooled$nu_spfk
  nu_spfk[1+nS,1+nP,1:nF,] <- rp$EBSBpars$totalAgg$nu_spfk

  # P1
  P1_spf[1:nS,1:nP,1:nF] <- rp$EBSBpars$stockSpec$P1_spf
  P1_spf[nS+1,1:nP,1:nF] <- rp$EBSBpars$coastWide$P1_spf
  P1_spf[1:nS,1+nP,1:nF] <- rp$EBSBpars$dataPooled$P1_spf
  P1_spf[1+nS,1+nP,1:nF] <- rp$EBSBpars$totalAgg$P1_spf

  # P2
  P2_spf[1:nS,1:nP,1:nF] <- rp$EBSBpars$stockSpec$P2_spf
  P2_spf[nS+1,1:nP,1:nF] <- rp$EBSBpars$coastWide$P2_spf
  P2_spf[1:nS,1+nP,1:nF] <- rp$EBSBpars$dataPooled$P2_spf
  P2_spf[1+nS,1+nP,1:nF] <- rp$EBSBpars$totalAgg$P2_spf

  depSeq <- seq(0,1,length.out = 100 )

  par( mfcol = c(nS+1,nP+1),
        mar = c(.5,.5,.5,.5),
        oma = c(3,3,2,2) )


  fleetCols <- RColorBrewer::brewer.pal(nF,"Set1")

  for( s in 1:(nS+1) )
    for( p in 1:(nP+1) )
    {
      plot( x = c(0,1), y = c(0,2), axes = FALSE,
            type = "n" )
        mfg <- par( "mfg" )

        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        if( mfg[2] == 1 )
          axis( side = 2, las = 1 )

        box()

        abline( h = 1, lwd = .8, lty = 2 )

        for( f in 1:nF )
        {
          numerator   <- (1 - exp( -nu_spfk[s,p,f,3]*( depSeq - P1_spf[s,p,f]))) 
          denominator <- (1 - exp( -nu_spfk[s,p,f,3]*( P2_spf[s,p,f] - P1_spf[s,p,f]))) 

          ratioSeq    <- nu_spfk[s,p,f,1] + (nu_spfk[s,p,f,2] - nu_spfk[s,p,f,1]) * numerator/denominator


          lines( x = depSeq, y = ratioSeq, lwd = 2, col = fleetCols[f] )
        }
        
    }

  legend( x = "topleft",
          legend = fleetNames,
          col = fleetCols,
          lwd = 2, bty = "n" )

} # END plotEBSBratio()


# Basically, a stacked up plot showing 
# simulated index data:
# 1. plotScaledIndices
# 2. q_t series (r_t as rectangles)
# 3. delta values
plotSimIdx <- function( obj = blob,
                        iRep = 1,
                        sIdx = 1, pIdx = 1)
{
  tMP       <- obj$om$tMP
  nT        <- obj$om$nT

  I_t       <- obj$mp$data$I_ispft[iRep,sIdx,pIdx,5,1:nT]
  delta_t   <- obj$om$errors$delta_ispft[iRep,sIdx,pIdx,5,1:nT]
  SB_t      <- obj$om$SB_ispt[iRep,sIdx,pIdx,1:nT]
  q_t       <- obj$om$q_ispft[iRep,sIdx,pIdx,5,1:nT]
  rI_ft     <- obj$om$rI_ispft[iRep,sIdx,pIdx,4:5,1:nT]

  tauObs_f  <- obj$ctlList$opMod$histRpt$tauComb_pg[pIdx,4:5]

  tau_t     <- c()
  for(t in 1:nT)
    tau_t[t] <- sum(rI_ft[,t] * tauObs_f)

  # recalc delta
  delta_t[1:(tMP-1)] <- (log(I_t[1:(tMP-1)]) - log(q_t[1:(tMP-1)]) - log(SB_t[1:(tMP-1)]))/tau_t[1:(tMP-1)]

  yrs <- seq(from = 1951, by = 1, length.out = nT)

  par(mfrow = c(3,1), 
      mar = c(.5,1,.5,1),
      oma = c(5,5,1,5))

  plot( x = range(yrs), y = c(0,max(SB_t,I_t)),
        type = "n", axes = FALSE)
    axis(side = 2, las = 1, cex.axis = 1.5)
    mtext(side = 2, text = "Spawning Biomass and Index (kt)", line = 3)
    box()
    grid()

    
    lines(x = yrs, y = SB_t, lwd = 3, col = "red")
    points(x = yrs, y = I_t/q_t, pch = 16, cex = 1.5, col = "grey40")

  plot( x = range(yrs), y = c(0,max(q_t,rI_ft,na.rm =T)),
        type = "n", axes = FALSE)
    axis(side = 2, las = 1, cex.axis = 1.5)
    axis(side = 4, las = 1, cex.axis = 1.5 )
    mtext(side = 2, text = "Catchability", line = 3)
    mtext(side = 4, text = expression(eta), line = 3)
    box()
    grid()

    rect( xleft = yrs - .3, xright = yrs + .3,
          ybottom = 0, ytop = rI_ft[1,], border = NA,
          col = "grey60" )
    lines(x = yrs, y = q_t, lwd = 1, col = "black")
    points(x = yrs, y = q_t, lwd = 1, pch = 16)

  plot( x = range(yrs), y = range(delta_t),
        type = "n", axes = FALSE)
    axis(side = 2, las = 1, cex.axis = 1.5)
    mtext(side = 2, text = "Std. Observation Residual", line = 3)
    axis(side = 1, cex.axis = 1.5)
    mtext( side = 1, text = "Year", line = 2)
    box()
    grid()

    points(x = yrs[1:(tMP-1)], y = delta_t[1:(tMP-1)],
            pch = 16, cex = 1.5, col = "grey40" )
    histMean <- mean(delta_t[1:(tMP-1)])
    histSD   <- sd(delta_t[1:(tMP-1)])
    segments(y0 = histMean, x0 = yrs[1], x1 = yrs[tMP-1], lwd = 1, lty = 1, col = "grey30")
    segments(y0 = histMean + c(1,-1)*histSD, x0 = yrs[1], x1 = yrs[tMP-1], lwd = 1, lty = 2, col = "grey30")

    points(x = yrs[tMP:nT], y = delta_t[tMP:nT],
            pch = 16, cex = 1.5, col = "salmon" )
    projMean <- mean(delta_t[tMP:nT])
    projSD   <- sd(delta_t[tMP:nT])
    segments(y0 = projMean, x0 = yrs[tMP], x1 = yrs[nT], lwd = 1, lty = 1, col = "salmon")
    segments(y0 = projMean + c(1,-1)*projSD, x0 = yrs[tMP], x1 = yrs[nT], lwd = 1, lty = 2, col = "salmon")

  
}

# Basically, a stacked up plot showing 
# simulated index data:
# 1. plotScaledIndices
# 2. Recruitment
# 3. Mortality
plotSimProj <- function(  obj = blob,
                          iRep = 1,
                          sIdx = 1, pIdx = 1)
{
  tMP       <- obj$om$tMP
  nT        <- obj$om$nT

  I_t       <- obj$mp$data$I_ispft[iRep,sIdx,pIdx,5,1:nT]
  SB_t      <- obj$om$SB_ispt[iRep,sIdx,pIdx,1:nT]
  M_t       <- obj$om$M_iaxspt[iRep,2,1,sIdx,pIdx,1:nT]
  R_t       <- obj$om$R_ispt[iRep,sIdx,pIdx,1:nT]
  q_t       <- obj$om$q_ispft[iRep,sIdx,pIdx,5,1:nT]

  yrs <- seq(from = 1951, by = 1, length.out = nT)

  par(mfrow = c(3,1), 
      mar = c(.5,1,.5,1),
      oma = c(5,5,1,2))

  plot( x = range(yrs), y = c(0,max(SB_t,I_t)),
        type = "n", axes = FALSE)
    axis(side = 2, las = 1, cex.axis = 1.5)
    mtext(side = 2, text = "Spawning Biomass and Index (kt)", line = 4)
    box()
    grid()
    points(x = yrs, y = I_t/q_t, pch = 16, cex = 1, col = "grey60")
    lines(x = yrs, y = SB_t, lwd = 3, col = "red")
    abline(v = yrs[tMP]-.5, lty = 2)

  plot( x = range(yrs), y = c(0,max(M_t)),
        type = "n", axes = FALSE)
    axis(side = 2, las = 1, cex.axis = 1.5)
    mtext(side = 2, text = "Age-2+ Mortality (/yr)", line = 4)
    box()
    grid()

    lines(x = yrs, y = M_t, lwd = 3, col = "salmon")
    abline(v = yrs[tMP]-.5, lty = 2)

  plot( x = range(yrs), y = c(0,max(R_t)),
        type = "n", axes = FALSE)
    axis(side = 2, las = 1, cex.axis = 1.5)
    mtext(side = 2, text = "Age-1 Recruits (1e6)", line = 4)
    axis(side = 1, cex.axis = 1.5)
    mtext( side = 1, text = "Year", line = 2)
    box()
    grid()
    abline(v = yrs[tMP]-.5, lty = 2)
    points(x = yrs, y = R_t, pch = 16, col = "grey40")
  
}

# plotRetroCatchability()
# Plot of retrospective catchability estimates for each 
# fleet compared to the mean taken from the OM
plotRetroCatchability <- function(  obj = blob,
                                    iRep = 1 )
{
  # Get model dims
  tMP <- obj$om$tMP
  nT  <- obj$om$nT
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nF  <- obj$om$nF

  # Get model state arrays and AM fits
  SB_spt        <- obj$om$SB_ispt[iRep,,,1:t]
  VB_spt        <- obj$om$vB_ispft[iRep,,,2,1:t]
  totB_spt      <- obj$om$B_ispt[iRep,,,1:t]
  fitSB_spt     <- obj$mp$assess$retroSB_itspt[iRep,projt,,,1:t]
  fitVB_spft    <- obj$mp$assess$retroVB_itspft[iRep,projt,,,,1:t]
  fitq_spf      <- obj$mp$assess$retroq_itspf[iRep,,,,]
  fitq_spft     <- obj$mp$assess$retroq_itspft[iRep,,,,,]

  # Get ctlList
  ctlList <- obj$ctlList
  pT      <- ctlList$opMod$pT

  # Get mean catchability (OM fits)
  mq_spf      <- ctlList$opMod$histRpt$q_spf
  mq_sf       <- array(1, dim = c(nS,nF))
  mq_sf[,3:4] <- ctlList$opMod$histRpt$qSurv_sf
  sdlnq_f     <- ctlList$mp$assess$spsdlnq_f


}

# Plot comp fits
plotAgeDataYrs <- function( obj = blob,
                            iRep = 1,
                            sIdx = 1,
                            pIdx = 1,
                            fIdx = 1 )
{
  A_at  <- obj$mp$data$A_iaxspft[iRep,,1,sIdx,pIdx,fIdx,]
  A_at[is.na(A_at)] <- -1
  # A_at[L_lxt < 0] <- NA

  pA_at <- A_at
  # Pull fleet names
  fleetNames_f <- obj$om$fleetNames_f

  # Pull model dims
  nA      <- obj$om$nA  
  nF      <- obj$om$nF
  tMP     <- obj$om$tMP
  nLenX   <- obj$om$nLenX
  nT      <- obj$om$nT
  fYear   <- obj$ctlList$opMod$histRpt$fYear
  
  # Make colours vector
  cols_f    <- brewer.pal( n = nF, "Paired" )

  # Make years vector
  years   <- seq(fYear, length = nT, by = 1)

  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  gearTimes <- c()
  for(t in 1:nT)
  {
    if( any(A_at[,t] > 0))
    {
      gearTimes <- c(gearTimes,t)

      pA_at[,t] <- A_at[,t]/max(sum(A_at[,t],na.rm = T),1)
    }

    if( all(A_at[,t] <= 0))
      pA_at[,t] <- NA          

  }
  pA_at[pA_at <= 0] <- NA
  A_at[A_at <= 0] <- NA

  Nobs_t <- apply(X = A_at, FUN = sum, MARGIN =2, na.rm = T)

  nObs <- length(gearTimes)
  nCols <- round(sqrt(nObs))
  nRows <- ceiling(nObs/nCols)

  par(  mfcol = c(nRows,nCols), 
        mar = c(0,0.1,0,0.1),
        oma = c(3,3,3,3) )

  yMax <- max( pA_at, na.rm = T )


  for( tIdx in gearTimes )
  { 
    plot( x = c(1,nA), y = c(0,yMax),
          xlab = "", ylab = "", type = "n", las = 1,
          axes = FALSE )
    box()
    grid()
    mfg <- par("mfg")
    if( mfg[1] == mfg[3] )
      axis( side = 1 )

    if( mfg[2] == 1 )
      axis( side = 2, las = 1 )

    if(tIdx <= tMP - 1)
      barCol <- "grey40"

    if(tIdx > tMP - 1)
      barCol <- "salmon"

    rect( xleft = 1:nA - .3, xright = 1:nA +.3,
          ybottom = 0, ytop = pA_at[,tIdx],
          col = barCol, border = NA )
    
    legend( x = "topright",
            legend = c( years[tIdx],
                        paste("N = ", round(Nobs_t[tIdx]), sep = "" )),
            bty = "n" )

  }
  mtext(  side = 3, 
          outer = T, 
          text = fleetNames_f[fIdx], 
          line = 2, font = 2)
  mtext(  side = 1, outer = T, text = "Age (yrs)", line = 2 )
  mtext(  side = 2, outer = T, text = "Proportion-at-Age", 
          line = 2 )

} # END plotCompFitYrs()



# plotScaledIndices()
# Scaled indices with model fits for a given
# replicate and year. 
plotScaledIndices <- function(  obj = blob, 
                                iRep = 1, 
                                Ct = TRUE,
                                t = blob$om$tMP )
{
  # Get model dims
  tMP <- obj$om$tMP
  nT  <- obj$om$nT
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nF  <- obj$om$nF

  # Calculate the projection time step
  projt <- t - tMP + 1

  # Blended index?
  blendIdx      <- obj$ctlList$opMod$blendIdx

  SB_spt        <- array(NA, dim = c(nS,nP,t))
  VB_spt        <- array(NA, dim = c(nS,nP,t))
  fitSB_spt     <- array(NA, dim = c(nS,nP,t))
  fitVB_spft    <- array(NA, dim = c(nS,nP,nF,t))
  fitq_spf      <- array(NA, dim = c(nS,nP,nF))
  fitq_spft     <- array(NA, dim = c(nS,nP,nF,t))
  fitqComb_spt  <- array(NA, dim = c(nS,nP,t))
  C_spft        <- array(NA, dim = c(nS,nP,nF,t))
  C_spt         <- array(NA, dim = c(nS,nP,t))

  # Get biomass arrays and catchabilities
  SB_spt[1:nS,,]        <- obj$om$SB_ispt[iRep,1:nS,,1:t]
  VB_spt[1:nS,,]        <- obj$om$vB_ispft[iRep,,,2,1:t]
  fitSB_spt[1:nS,,]     <- obj$mp$assess$retroSB_itspt[iRep,projt,,,1:t]
  fitVB_spft[1:nS,,,]   <- obj$mp$assess$retroVB_itspft[iRep,projt,,,,1:t]
  fitq_spf[1:nS,,]      <- obj$mp$assess$retroq_itspf[iRep,projt,,,]
  fitq_spft[1:nS,,,]    <- obj$mp$assess$retroq_itspft[iRep,projt,,,,1:t]
  fitqComb_spt[1:nS,,]  <- obj$mp$assess$retroqComb_itspt[iRep,projt,,,1:t]

  ctlList <- obj$ctlList
  spFleets <- ctlList$mp$assess$spFleets


  fitSB_spt[fitSB_spt < 0] <- NA

  C_spft[1:nS,,,]   <- obj$om$C_ispft[iRep,,,,1:t]
  C_spt[1:nS,,]     <- apply(X = C_spft, FUN = sum, MARGIN = c(1,2,4))

  # Now pull indices
  I_spft <- obj$mp$data$I_ispft[iRep,,,,1:t]
  I_spft[I_spft < 0] <- NA

  # Get ratios for blending index
  rI_spft <- obj$om$rI_ispft[iRep,,,,1:t]

  fleetPCH  <- 20 + 1:nF
  fleetBG   <- RColorBrewer::brewer.pal(nF, "Set1")
  stockCol  <- RColorBrewer::brewer.pal(nP, "Spectral")


  spTVqFleets <- ctlList$mp$assess$spTVqFleets

  SB_spt[SB_spt == 0]     <- NA
  VB_spt[VB_spt == 0]     <- NA

  nSS     <- dim( SB_spt)[1]
  nPP     <- dim( SB_spt)[2]



  scaledIdx_spft <- array( NA, dim = c(nSS, nP, nF, t ) )
  scaledIdx_spt  <- array( NA, dim = c(nSS, nP, t ) )
  for( s in 1:nSS )
    for( p in 1:nPP )
    {
      for( f in 1:nF )
      {
        if( f %in% spTVqFleets)
          scaledIdx_spft[s,p,f,] <- I_spft[s,p,f,] / fitq_spft[s,p,f,] 
        else
          scaledIdx_spft[s,p,f,] <- I_spft[s,p,f,] / fitq_spf[s,p,f]

        # Scale by ratio of SB and VB
        scaledIdx_spft[s,p,f,] <- scaledIdx_spft[s,p,f,] * fitSB_spt[s,p,] / fitVB_spft[s,p,f,]
      }

      if( blendIdx )
      {
        scaledIdx_spt[s,p,]   <- I_spft[s,p,5,] / fitqComb_spt[s,p,]
      }
    }


  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- obj$ctlList$opMod$stock
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear
  pT            <- obj$ctlList$opMod$pT

  yrs <- seq( from = fYear, by = 1, length.out = t)
  ppJitter <- seq(from = -.3, to = .3, length.out = nP )

  par(  mfcol = c(nPP,nSS), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )
  for(s in 1:nSS)
    for( p in 1:nPP )
    {
      plot( x = range(yrs),
            y = c(0,max(VB_spt[s,p,],SB_spt[s,p,],na.rm = T) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+1.5, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()

      if( Ct )
      {
        # Plot actual catch
        rect( xleft = yrs - .3, xright = yrs + .3, 
              ytop = C_spt[s,p,], ybottom = 0, col = "grey40",
              border = NA )

      }

      lines( x = yrs, y = SB_spt[s,p,], col = "red", lwd = 3 )
      lines( x = yrs, y = fitSB_spt[s,p,], col = "grey60", lwd = 1 )

      if( !blendIdx )
        for( f in 1:nF )
          points( x = yrs, y = scaledIdx_spft[s,p,f,],
                  pch = fleetPCH[f], bg = fleetBG[f],
                  cex = 1.3, col = stockCol[p] )
        
      if( blendIdx )
        points( x = yrs, y = scaledIdx_spt[s,p,],
                pch = 16, col = "grey40" )
  
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }
    # legend( x= "topright",
    #         legend = fleetNames,
    #         pch = fleetPCH,
    #         pt.bg = fleetBG )

} # END plotScaledIndices()

# plotAMIdxResids()
# Plots of standardised AM residuals
plotAMIdxResids <- function(  obj = blob, 
                              iRep = 1, 
                              Ct = TRUE,
                              t = blob$om$tMP )
{
  # Get model dims
  tMP <- obj$om$tMP
  nT  <- obj$om$nT
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nF  <- obj$om$nF

  blendIdx      <- obj$ctlList$opMod$blendIdx

  # Calculate the projection time step
  projt <- t - tMP + 1

  SB_spt        <- array(NA, dim = c(nS,nP,t))
  VB_spt        <- array(NA, dim = c(nS,nP,t))
  fitSB_spt     <- array(NA, dim = c(nS,nP,t))
  fitVB_spft    <- array(NA, dim = c(nS,nP,nF,t))
  fitq_spf      <- array(NA, dim = c(nS,nP,nF))
  fitq_spft     <- array(NA, dim = c(nS,nP,nF,t))
  I_spft        <- array(NA, dim = c(nS,nP,nF,t))
  tauObs_spf    <- array(NA, dim = c(nS,nP,nF))

  # Get biomass arrays
  SB_spt[1:nS,,1:t]     <- obj$om$SB_ispt[iRep,,,1:t]
  VB_spt[1:nS,,1:t]     <- obj$om$vB_ispft[iRep,,,2,1:t]
  fitSB_spt[1:nS,,1:t]  <- obj$mp$assess$retroSB_itspt[iRep,projt,,,1:t]
  fitVB_spft[1:nS,,,1:t]<- obj$mp$assess$retroVB_itspft[iRep,projt,,,,1:t]
  fitq_spf[1:nS,,]      <- obj$mp$assess$retroq_itspf[iRep,projt,,,]
  fitq_spft[1:nS,,,1:t] <- obj$mp$assess$retroq_itspft[iRep,projt,,,,1:t]
  
  tauObs_spf[1:nS,,]    <- obj$mp$assess$retrotauObs_itspf[iRep,projt,,,]

  if( all( is.na(tauObs_spf)) )
    tauObs_spf[1:nS,,] <- 1

  tauObs_spf[is.na(tauObs_spf)] <- 0

  ctlList <- obj$ctlList

  fitSB_spt[fitSB_spt < 0] <- NA

  # idxFleets


  spTVqFleets <- ctlList$mp$assess$spTVqFleets
  spFleets    <- ctlList$mp$assess$spFleets

  # Model dims
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP 
  nT      <- obj$om$nT

  # Now pull indices
  I_spft[1:nS,1:nP,,1:t] <- obj$mp$data$I_ispft[iRep,1:nS,1:nP,,1:t]
  I_spft[I_spft <= 0] <- NA

  nSS <- nS
  nPP <- nP

  stdResids_spft  <- array( NA, dim = c(nSS, nPP, nF, t) )
  stdResids_spt   <- array( NA, dim = c(nSS, nPP, t) )

  for( s in 1:nSS )
    for( p in 1:nPP )
    {
      for( f in spFleets )
      {
        if( tauObs_spf[s,p,f] > 0 )
        {

          if( f %in% spTVqFleets )
            stdResids_spft[s,p,f,] <- (-log( I_spft[s,p,f,]/fitq_spft[s,p,f,] ) + log(fitVB_spft[s,p,f,]))/tauObs_spf[s,p,f]
          if( !f %in% spTVqFleets ) 
            stdResids_spft[s,p,f,] <- (-log( I_spft[s,p,f,]/fitq_spf[s,p,f] ) + log(fitVB_spft[s,p,f,]))/tauObs_spf[s,p,f]


        }
      }
    }

  
  fleetPCH <- 20 + 1:nF
  fleetBG <- RColorBrewer::brewer.pal(nF, "Set1")



  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- obj$ctlList$opMod$stock
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear
  pT            <- obj$ctlList$opMod$pT

  # HACK - need to import fleetnames
  fleetNames <- c("reduction","seineRoe","gillnet","surf","dive","SOK")

  if( nSS == 1 )
    speciesNames <- "Data Pooled"

  yrs <- seq( from = fYear, by = 1, length.out = t)

  ppJitter <- seq( from = -.3, to = .3, length.out = nP )

  par(  mfcol = c(nPP,nSS), 
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )
  for(s in 1:nSS)
    for( p in 1:nPP )
    {
      if( any(!is.na(stdResids_spft[s,p,,])))
        maxResid <- max(abs(stdResids_spft[s,p,,]),na.rm = T)
      else maxResid <- 1

      plot( x = range(yrs),
            y = range(-maxResid,maxResid),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+1.5, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()

      for( f in 1:nF )
      {
        points( x = yrs, y = stdResids_spft[s,p,f,],
                pch = fleetPCH[f], bg = fleetBG[f],
                cex = 1.3 )
        nonNA <- which(!is.na(stdResids_spft[s,p,f,]))
        if( length(nonNA) > 0 )
        {
          yVal <- stdResids_spft[s,p,f,nonNA]
          xVal <- yrs[nonNA]



          dat <- data.frame(x = xVal, y = yVal )
          regLine <- lm( y~x, data = dat )

          pVal <- round(summary(regLine)$coefficients[2,4],3)

          pLabel <- paste("p = ", pVal, sep = "")

          dat <- dat %>%
                  mutate( regLine = predict.lm(regLine, newdata = dat) )


          lines( x = dat$x, y = dat$regLine, col = fleetBG[f], lwd = 2 )
          text( x = dat$x[1], y = 1.2*dat$y[1], label = pLabel, col = fleetBG[f], font = 2 )

        }
      }
      abline( h = 0, lty = 2)
      
      
    }
    
    legend( x= "bottomleft", bty = "n",
            legend = c(fleetNames),
            pch = fleetPCH,
            pt.bg = fleetBG,
            cex = 1.3 )



} # END plotAMIdxResids()


.plotDiagCondition <- function( obj = blob, iRep = 1)
{
  repObj <- obj$ctlList$opMod$histRpt
  posts <- obj$ctlList$opMod$posts

  iPost <- ifelse( obj$ctlList$opMod$posteriorSamples, obj$ctlList$opMod$postDraws_i[iRep], 0 )

  diagCondition( repObj, obj, posts, iRep, iPost )

}

# diagCondition()
# Plots that help diagnose issues with conditioning
# between the fitted OM report, and the conditioned
# ms3R operating model.
diagCondition <- function(  repObj  = blob$ctlList$opMod$histRpt,
                            ms3Obj  = blob,
                            posts   = blob$ctlList$opMod$posts,
                            iRep    = 1,
                            iPost   = 0 )
{

  nP <- repObj$nP
  nT <- repObj$nT
  nF <- repObj$nG
  stockNames    <- dimnames(ms3Obj$ctlList$opMod$histRpt$I_pgt)[[1]]

  par(mfrow = c(3,2), mar = c(1,2,1.5,1), oma = c(3,5,3,3) )

  # Biomass RE
  plotRE_spt( repObj = repObj,
              omObj = ms3Obj$om, 
              posts = posts,
              AMseries = ifelse(iPost,"SB_ipt","SB_pt"), 
              OMseries = "SB_ispt", 
              iRep = iRep,
              iPost = iPost )
  mtext( side = 3, text = "SB_spt", line = 0, font = 2)

  # Recruitment
  plotRE_spt( repObj = repObj, 
              omObj = ms3Obj$om, 
              posts = posts,
              AMseries = ifelse(iPost,"R_ipt","R_pt"), 
              OMseries = "R_ispt", 
              iRep = iRep,
              iPost = iPost )
  mtext( side = 3, text = "R_spt", line = 0, font = 2)  

  # Recruitment errors
  plotRE_spt( repObj = repObj, 
              omObj = ms3Obj$om$errors, 
              nS = 1,
              posts = posts,
              AMseries = ifelse(iPost,"SRdevs_ipt","SRdevs_pt"),
              OMseries = "omegaR_ispt", 
              iRep = iRep,
              iPost = iPost,
              yRange = c(-1,1) )
  mtext( side = 3, text = expression(omega[R]), line = 0, font = 2)  

  # Catch
  plotRE_spft(  repObj = repObj, 
                omObj = ms3Obj$om, 
                AMseries = "C_pgt",
                OMseries = "C_ispft", 
                iRep = iRep,
                iPost = 0 )
  mtext( side = 3, text = expression(C[spft]), line = 0, font = 2)    

  # F
  plotRE_spft(  repObj = repObj, 
                omObj = ms3Obj$om, 
                posts = posts,
                AMseries = ifelse(iPost,"U_ipgt","U_pgt"),
                OMseries = "F_ispft", 
                iRep = iRep,
                iPost = iPost,
                legendON=TRUE, stocks=stockNames )
  mtext( side = 3, text = expression(F[spft]), line = 0, font = 2) 

  # M
  plotRE_at(  repObj = repObj,
                omObj = ms3Obj$om,
                AMseries = ifelse(iPost,"M_iapt","M_apt"),
                OMseries = "M_iaxspt",
                iRep = iRep,
                posts = posts,
                iPost = iPost,
                stocks = stockNames )

  mtext( side = 3, text = expression(M[at]), line = 0, font = 2) 

  mtext( side = 1, outer = TRUE, text = "Time Step", font = 2, line = 2)
  mtext( side = 2, outer = TRUE, text = "Relative error in conditioning", 
          font = 2, line = 3 )



  # # Numbers at age
  # plotRE_axspt( repObj = repObj, omObj = ms3Obj$om, series = "N_iaxspt" )
  # mtext( side = 3, text = expression(N[axspt]), line = 0, font = 2)    

  # # Biomass at age
  # plotRE_axspt( repObj = repObj, omObj = ms3Obj$om, series = "N_axspt" )
  # mtext( side = 3, text = expression(B[axspt]), line = 0, font = 2)    
}

calcREdist <- function( est, true, marg, 
                        qProbs = c(0.025, 0.5, 0.975) )
{
  # Calculate REs
  re <- (est - true)/true

  # calculate distribution over margin
  reQuants <-  apply( X = re, FUN = quantile, 
                      MARGIN = marg, probs = qProbs,
                      na.rm = TRUE )

  reQuants[!is.finite(reQuants)] <- 0

  return(reQuants)

}

plotRE_spt <- function( repObj, omObj, posts=NULL, nS = 1,
                        AMseries = "SB_spt",
                        OMseries = "SB_ispt", iRep = 1,
                        iPost = 0,
                        yRange = NULL )
{

  nP <- repObj$nP
  nT <- repObj$nT

  true_spt  <- array(NA, dim = c(nS,nP,nT))
  est_spt   <- array(NA, dim = c(nS,nP,nT))

  if( iPost ) 
    true_spt[1:nS,,] <- posts[[AMseries]][iPost, ,1:nT,drop = FALSE]
  else
    true_spt[1:nS,,] <- repObj[[AMseries]][,1:nT,drop = FALSE]

  est_spt[1:nS,1:nP,1:nT] <- omObj[[OMseries]][iRep,1:nS,1:nP,1:nT]

  if( OMseries == "omegaR_ispt" )
  {
    
    if(is.null(repObj$avgRcode_p))
      repObj$avgRcode_p <- rep(0,nP)

    for (p in 1:nP)
    {
      
#      if(repObj$avgRcode_p[p]==1)
#        true_spt[1:nS,p,]   <- repObj$SRdevs_pt[p,1:nT]
#      else
#        true_spt[1:nS,p,]   <- repObj$omegaR_pt[p,1:nT]

    }
    sigmaR <- repObj$sigmaR
    est_spt[1:nS,1:nP,1:nT] <- est_spt[1:nS,1:nP,1:nT] - 0.5 * sigmaR
  }

  re_qspt  <- calcREdist( true = true_spt,
                          est  = est_spt,
                          marg = c(1,2,3) )

  stockCols <- RColorBrewer::brewer.pal(n = nP, "Dark2")
  stockLty <- rep(1,nP)

  if( is.null(yRange) )
  {
    yRange <- range(re_qspt, na.rm = TRUE)
    yRange[2] <- max(.01,yRange[2])
    yRange[1] <- min(-.01,yRange[1])
  }

  plot( x = c(1,nT), y = yRange,
        axes = FALSE, type = "n" )
    mfg <- par("mfg")
    if(mfg[1] == mfg[3])
      axis( side = 1 )
    axis( side = 2, las = 1 )
    box()
    grid()
    for( s in 1:nS )
      for( p in 1:nP )
      {
        # polyY <- c(re_qspt[1,s,p,],rev(re_qspt[3,s,p,]))
        # polygon( x = c(1:nT,nT:1), y = polyY,
        #           col = scales::alpha(stockCols[p],alpha = .3),
        #           lty = stockLty[p] )
        lines( x = 1:nT, y = re_qspt[2,s,p,], 
                col = stockCols[p], lty = stockLty[p] )

      }
    abline( h = 0, lty = 3, lwd = .8 )
}

plotRE_spft <- function(  repObj, omObj, posts=NULL,
                          OMseries = "C_spft",
                          AMseries = "C_pgt",
                          iRep = 1, iPost=0, nS = 1,
                          legendON=FALSE, stocks=NULL )
{
  nP <- repObj$nP
  nT <- repObj$nT
  nF <- repObj$nG

  # if(OMseries == "F_ispft")
  #   browser()

  fleetType_f <- omObj$fleetType_f
  sokFleets <- which(fleetType_f >= 2)

  true_spft  <- array(NA, dim = c(nS,nP,nF,nT))
  est_spft   <- array(NA, dim = c(nS,nP,nF,nT))

  if( iPost ) 
    true_spft[1:nS,,,] <- posts[[AMseries]][iPost,,,1:nT,drop = FALSE]  
  else
    true_spft[1:nS,,,] <- repObj[[AMseries]][,,1:nT,drop = FALSE]  


  est_spft[1:nS,1:nP,1:nF,1:nT] <- omObj[[OMseries]][iRep,1:nS,1:nP,1:nF,1:nT]


  # repObj C_pgt contains SOK product for SOK fleets
  # omObj C_spft contains ponded fish for SOK fleets

  # # convert SOK product which is what is reported in AMseries
  # if( OMseries == "C_ispft" & length(sokFleets) > 0 )
  # {
  #   psi_pgt <- repObj$psi_pgt
  #   psi_gt  <- repObj$psi_gt

  #   for (p in 1:nP)
  #     for( f in sokFleets)
  #     {
  #       if(any(psi_pgt[p,f,] != 0))
  #         true_spft[1:nS,p,f,] <- true_spft[1:nS,p,f,]/psi_pgt[p,f,]

  #       if(any(psi_gt[f,] != 0))
  #         true_spft[1:nS,p,f,] <- true_spft[1:nS,p,f,]/psi_gt[f,]
  #     }
  # }

  re_qspt  <- calcREdist( true = true_spft,
                            est  = est_spft,
                            marg = c(1,2,4) )

  re_qft  <- calcREdist(  true = true_spft,
                          est  = est_spft,
                          marg = c(3:4) )



  stockCols <- RColorBrewer::brewer.pal(n = nP, "Dark2")
  stockLty <- rep(1,nP)

  yRange <- range(re_qspt, na.rm = TRUE)
  yRange[2] <- max(.01,yRange[2])
  yRange[1] <- min(-.01,yRange[1])

  col <- RColorBrewer::brewer.pal(n = nF, "Paired")

  plot( x = c(1,nT), y = yRange,
        axes = FALSE, type = "n" )
    mfg <- par("mfg")
    if(mfg[1] == mfg[3])
      axis( side = 1 )
    axis( side = 2, las = 1 )
    box()
    grid()

    for( s in 1:nS )
      for( p in 1:nP )
      {
        polyY <- c(re_qspt[1,s,p,],rev(re_qspt[3,s,p,]))
        polygon( x = c(1:nT,nT:1), y = polyY,
                  col = scales::alpha(stockCols[p],alpha = .3),
                  lty = stockLty[p] )
        lines( x = 1:nT, y = re_qspt[2,s,p,], 
                col = stockCols[p], lty = stockLty[p] )

      }
    abline( h = 0, lty = 3, lwd = .8 )

    for(f in 1:nF )
      lines(x = 1:nT, y=re_qft[2,f,], col = col[f])

  if(legendON)
    legend('topright', bty='n',
           legend=1:nF, lty=1, col=col)

}

plotRE_at <- function(  repObj, omObj, 
                          OMseries = "M_iaxspt",
                          AMseries = "M_apt",
                          posts = NULL,
                          iRep = 1, iPost=0, nS = 1,
                          legendON=FALSE, stocks=NULL )
{
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT

  if( iPost ) 
    true_at <- posts[[AMseries]][iPost,,1,1:nT]
  else
    true_at <- repObj[[AMseries]][,1,1:nT]

  est_at  <- omObj[[OMseries]][iRep,,1,1,1,1:nT]

  re_qt  <- calcREdist( true = true_at[,1:nT],
                        est  = est_at[,1:nT],
                        marg = c(2) )

  stockLty <- rep(1,nP)

  stockCols <- RColorBrewer::brewer.pal(n = nP, "Dark2")

  yRange <- range(re_qt, na.rm = TRUE)
  yRange[2] <- max(.01,yRange[2])
  yRange[1] <- min(-.01,yRange[1])

  plot( x = c(1,nT), y = yRange,
        axes = FALSE, type = "n" )
    mfg <- par("mfg")
    if(mfg[1] == mfg[3])
      axis( side = 1 )
    axis( side = 2, las = 1 )
    box()
    grid()
    polyY <- c(re_qt[1,],rev(re_qt[3,]))
    polygon( x = c(1:nT,nT:1), y = polyY,
              col = scales::alpha(stockCols[1],alpha = .3),
              lty = stockLty[1] )
    lines( x = 1:nT, y = re_qt[2,], 
            col = stockCols[1], lty = stockLty[1] )


    abline( h = 0, lty = 3, lwd = .8 )
}

# So, x-axis is m1, colour is Mb, and pch is sigmaM
assignPCH <- function( Mb )
{
  pchList <- c(21,23,24)
  names(pchList) <- c("0.2","0.4","0.6")

  pch <- pchList[as.character(Mb)]
  pch
}

assignCol <- function( sigmaM, nameVec = c("0","0.1","0.3") )
{
  cols <- brewer.pal("Dark2", n = length(nameVec))[1:length(nameVec)]
  names(cols) <- nameVec

  col <- cols[as.character(sigmaM)]

  col
}

addJitter <- function( Mb, jitt = c(-.05,0,.05), nameVec = c("0.2","0.4","0.6") )
{
  names(jitt) <- nameVec

  xJitt <- jitt[as.character(Mb)]

  xJitt
}


plotTable <- function(  medREsTable = HG_MREs, ylab = "MRE",
                        columns = c("m1","Mb","h","B0","totB0"))
{

  medREsTable <- medREsTable |>
                  mutate( m1sim = as.numeric(m1sim),
                          Mbsim = as.numeric(Mbsim),
                          sigmaMsim = as.numeric(sigmaMsim),
                          pch = sapply( X = Mbsim, FUN = assignPCH),
                          col = sapply( X = sigmaMsim, FUN = assignCol),
                          MbJitter = sapply(X = Mbsim, FUN = addJitter),
                          sigmaMJitter = sapply(  X = sigmaMsim, FUN = addJitter, 
                                                  nameVec = c("0","0.1","0.3"), jitt = c(-.2,0.,.2)), 
                          xJitter = m1sim + MbJitter + sigmaMJitter  )

  nReps <- as.integer(medREsTable$nReps)
  medREsTable$bg <- medREsTable$col
  medREsTable$bg[nReps < 95] <- "white"
  medREsTable$pch[nReps < 75] <- NA


  nPlots <- length(columns)
  par(mfrow = c(nPlots,1),
      mar = c(.1,2,.1,2),
      oma = c(4,4,2,2) )

  m1Range <- c(0.5,7)

  m1s <- unique(medREsTable$m1sim)

  for( cIdx in 1:nPlots )
  {
    
    yRange <- abs(range(as.numeric(medREsTable[nReps >= 75,columns[cIdx]])))

    yRange <- c(-max(yRange),max(yRange))
    plot( x = m1Range, y = yRange,
          axes = FALSE, type = "n" )
      mfg <- par("mfg")
      axis( side = 2, las = 1)
      if( mfg[1] == mfg[3] )
        axis( side = 1, labels = m1s, at = m1s )
      grid()
      box()

      mtext(side = 2, text = paste0(ylab," (", columns[cIdx] ,")"),
            line = 3)

      # Plot points now, see about a jitter later
      points( x = medREsTable$xJitter, y = medREsTable[,columns[cIdx]],
              col = medREsTable$col, pch = medREsTable$pch,
              bg = medREsTable$bg,
              cex = 1.5 )
      abline(h = 0, lty = 2, col = "grey60")

      if(mfg[1] == 1 )
        legend( x = "topright", bty = "n",
                legend = c( "sigmaM = 0",
                            "sigmaM = 0.1",
                            "sigmaM = 0.3",
                            "Mb = 0.2",
                            "Mb = 0.4",
                            "Mb = 0.6"),
                col = c(brewer.pal("Dark2", n = 3),
                        rep("grey40",3)),
                bg = c(brewer.pal("Dark2", n = 3),
                        rep("white",3)),
                pch = c(15,15,15,21,23,24),
                cex = 1.2 )

  }
}

# plotSimEstDepM()
# Plots the true OM depM model,
# and the tulip of estimated depM
# models from a sim-est experiment
plotSimEstDepM <- function( obj = blob, nB = 100,
                            nTrace = 3,
                            sIdx = 1, pIdx = 1, pt = 1)
{
  # First, pull out OM values
  totB0     <- obj$rp[[1]]$totB0_sp[sIdx,pIdx]
  Mb        <- obj$rp[[1]]$Mb_sp[sIdx,pIdx]
  m1        <- obj$rp[[1]]$m1_sp[sIdx,pIdx]

  # then pull retros
  retroMb_i <- obj$mp$assess$retroM_itsp[,pt,sIdx,pIdx]
  retrom1_i <- obj$mp$assess$retrom1_itsp[,pt,sIdx,pIdx]

  nReps <- sum(obj$goodReps)
  goodReps <- which(obj$goodReps)

  Dseq <- seq(from = 0, to = 3, length.out = nB)

  depMcurve_id <- array(NA, dim = c(nReps,nB))
  for( i in 1:length(goodReps) )
  {
    idx <- goodReps[i]

    depMcurve_id[i,] <- retroMb_i[idx] + exp( -retrom1_i[idx] * Dseq )
  }

  traceIdx <- sample(1:nReps,nTrace)

  trueDepM <- Mb + exp(-m1 * Dseq )

  depMcurve_qd <- apply(  X = depMcurve_id, FUN  = quantile, 
                          MARGIN = 2, probs = c(0.025, 0.5, 0.975) )


  maxM <- max(c(Mb + 1, retroMb_i + 1), na.rm = T)

  plot( x = c(0,3), y = c(0,maxM), type = "n",
        las = 1 )
    polygon(  x = c(Dseq,rev(Dseq)),
              y = c(depMcurve_qd[1,],rev(depMcurve_qd[3,])),
              border = NA, col = "grey70" )
    lines( x = Dseq, y = depMcurve_qd[2,], lwd = 3, col = "black", lty = 2 )
    lines( x = Dseq, y = trueDepM, col = "red", lwd = 3 )
    for( tIdx in traceIdx )
      lines( x = Dseq, y = depMcurve_id[tIdx,], lwd = .8)

}


# plotGridTulipBtCt()
# Plots a grid of simulation envelopes
# for biomass and legal catch (TAC)
plotGridTulipBtCtUt <- function(  blobList = mpBlobList,
                                  labels = NULL,
                                  folder = "",
                                  dep   = TRUE,
                                  yLimB   = c(0,2),
                                  yLimC   = c(0,15),
                                  yLimHR  = c(0,.3),
                                  traces = 3,
                                  traceSeed = 1234,
                                  proj  = TRUE,
                                  maxProjT = NULL,
                                  goodReps = NULL,
                                  yrRange = -(150:50),
                                  empRefCurves = NULL)
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

  if(is.null(goodReps))
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
    if(!is.null(empRefCurves))
    {
      eqListQuant <- solveSimEqbriaQuant( empRefCurves = empRefCurves,
                                          maxXspline = maxXspline, 
                                          yrRange = yrRange,
                                          USR_q = Busr_q, reps = 1:nReps,
                                          qProbs = c(0.25, 0.5, 0.75) )  


      B0_isp[,1,1]     <- eqListQuant$B0_q[2]
      Bcrash_isp[,1,1] <- eqListQuant$Bcrash_q[2]
      Ucrash_isp[,1,1] <- eqListQuant$Ucrash_q[2]
      MSY_isp[,1,1]    <- eqListQuant$MSY_q[2]
      Bmsy_isp[,1,1]   <- eqListQuant$Bmsy_q[2]
    } else {
      B0_isp    <- blob$om$B0_isp[1:nReps,,,drop = FALSE]
      Bmsy_isp  <- blob$om$Bmsy_isp[1:nReps,,,drop = FALSE]
      Umsy_isp  <- blob$om$Umsy_isp[1:nReps,,,drop = FALSE]
    }



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
}


# plotGridTulipBtCt()
# Plots a grid of simulation envelopes
# for biomass and legal catch (TAC)
plotGridTulipBtRtMt <- function(  blobList  = sims_CC,
                                  labels    = NULL,
                                  folder    = "",
                                  dep       = TRUE,
                                  yLimB     = c(0,2),
                                  yLimR     = c(0,2000),
                                  yLimM     = c(0,1.4),
                                  yLimSBPR  = c(0,.8),
                                  traces = 3,
                                  traceSeed = 1234,
                                  proj  = TRUE,
                                  maxProjT = NULL,
                                  empRefCurves = NULL)
{
  # Load blobs
  nSims <- length(blobList)
  simLabs <- c()
  nTs     <- c()
  for(simIdx in 1:nSims)
  {
    simLabs[simIdx] <- paste0(blobList[[simIdx]]$ctlList$ctl$mpName,":",blobList[[simIdx]]$ctlList$ctl$scenarioName)
    names(blobList)[simIdx] <- simLabs[simIdx]
    nTs[simIdx] <- blobList[[simIdx]]$om$nT
  }

  if(!is.null(labels))
    plotLabs <- labels
  else plotLabs <- simLabs


  # Now set up plotting area
  par(  mfcol = c(4,nSims),
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


  nReps <- max(which(blobList[[1]]$goodReps))
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
    SB_ispt       <- blobList[[simIdx]]$om$SB_ispt[1:nReps,,,,drop = FALSE]
    R_ispt        <- blobList[[simIdx]]$om$R_ispt[1:nReps,,,,drop = FALSE]
    M_it          <- blobList[[simIdx]]$om$M_iaxspt[1:nReps,2,1,1,1,]

    SSBpr_ispt    <- SB_ispt/R_ispt
    
    
    # Load ref pts
    B0_isp     <- array(0,dim = c(nReps,nS,nP))
    R0_isp     <- array(0,dim = c(nReps,nS,nP))
    bcR0_isp   <- array(0,dim = c(nReps,nS,nP))
    bcM0_isp   <- array(0,dim = c(nReps,nS,nP))
    phi0_isp   <- array(0,dim = c(nReps,nS,nP))
    M0_isp     <- array(0,dim = c(nReps,nS,nP))
    Bmsy_isp   <- array(0,dim = c(nReps,nS,nP))
    USR_isp    <- array(0,dim = c(nReps,nS,nP))
    Bcrash_isp <- array(0,dim = c(nReps,nS,nP))
    Ucrash_isp <- array(0,dim = c(nReps,nS,nP))
    MSY_isp    <- array(0,dim = c(nReps,nS,nP))
    # Umsy_isp  <- array(0,dim = c(nReps,nS,nP))
    # MSY_isp   <- array(0,dim = c(nReps,nS,nP))

    for( i in 1:nReps)
    {
      # Pull out Dstar and Dscalar
      tdxUSR        <- blob$ctlList$opMod$histCtl$hypo$Dstar - 1950
      scalarUSR     <- blob$ctlList$opMod$histCtl$hypo$Dscalar
      
      # First, calculate the USR
      USR_isp[i,,]  <- apply(X = scalarUSR * SB_ispt[i,,,tdxUSR,drop = FALSE], FUN = mean, MARGIN = 2:3 )
      
      B0_isp[i,,]   <- blob$rp[[i]]$B0_sp[,,1]
      R0_isp[i,,]   <- blob$rp[[i]]$R0_sp
      M0_isp[i,,]   <- blob$rp[[i]]$M_xsp[1,,]


      bcR0_isp[i,,] <-  R0_isp[i,,] * exp(-0.5*blob$om$sigmaR_isp[i,,]^2)
      bcM0_isp[i,,] <-  M0_isp[i,,] * exp(-0.5*blob$om$sigmaM_isp[i,,]^2)

    
      # Umsy_isp[i,,] <- blob$rp[[i]]$FmsyRefPts$lUmsy_sp
      # MSY_isp[i,,]  <- blob$rp[[i]]$FmsyRefPts$lYeqFmsy_sp
    }

    if(!is.null(empRefCurves))
    {

      Busr_q <- quantile(USR_isp[,1,1], probs = c(0.25, 0.5, 0.75))
      eqListQuant <- solveSimEqbriaQuant( empRefCurves = empRefCurves,
                                          maxXspline = .3, 
                                          yrRange = -c(150:50),
                                          USR_q = Busr_q, reps = 1:nReps,
                                          qProbs = c(0.25, 0.5, 0.75) )  


      B0_isp[,1,1]     <- eqListQuant$B0_q[2]
      Bcrash_isp[,1,1] <- eqListQuant$Bcrash_q[2]
      Ucrash_isp[,1,1] <- eqListQuant$Ucrash_q[2]
      MSY_isp[,1,1]    <- eqListQuant$MSY_q[2]
      Bmsy_isp[,1,1]   <- eqListQuant$Bmsy_q[2]

      # simEqList <- solveSimEqbriaQuant( empRefCurves = empRefCurves,
      #                                   maxXspline = 0.5, USR = USR_isp[i,1,1], 
      #                                   iRep = i)
      # browser()
      # Bcrash_isp[i,,] <- simEqList$Bcrash
      # Ucrash_isp[i,,] <- simEqList$Ucrash
      # Bmsy_isp[i,,]   <- simEqList$Bmsy
      # MSY_isp[i,,]    <- simEqList$MSY
    }

    LRP_isp <- 0.3*B0_isp
    phi0_isp <- B0_isp/R0_isp
    # USR_isp <- 0.8*Bmsy_isp
    


    sbylab <- "Spawning Biomass (kt)"

    if(dep)
    {
      sbylab <- expression(SB[t]/SB[0])
      for( t in 1:nTs[simIdx] )
      {
        SB_ispt[,,,t] <- SB_ispt[,,,t] / B0_isp
      }
      LRP_isp <- LRP_isp/B0_isp
      Bcrash_isp <- Bcrash_isp/B0_isp
      Bmsy_isp <-  Bmsy_isp/B0_isp
      USR_isp <- USR_isp/B0_isp
      B0_isp  <- B0_isp/B0_isp
      # USR_isp <- USR_isp/Bmsy_isp
      # TRP_isp <- TRP_isp/Bmsy_isp
    }
  
    B0 <- median(B0_isp)
    LRP <- median(LRP_isp)
    Bcrash <- median(Bcrash_isp)
    Ucrash <- median(Ucrash_isp)
    Bmsy <- median(Bmsy_isp)
    USR <- median(USR_isp)
    MSY <- median(MSY_isp)
    # TRP <- mean(TRP_isp)

    SB_qspt   <- apply(X = SB_ispt, FUN = quantile, MARGIN = 2:4,
                        probs = c(0.025, 0.5, .975 ) )

    mnSB_spt   <- apply(X = SB_ispt, FUN = mean, MARGIN = 2:4)

    R_qspt   <- apply(X = R_ispt, FUN = quantile, MARGIN = 2:4,
                        probs = c(0.025, 0.5, .975 ) )

    mnR_spt   <- apply(X = R_ispt, FUN = mean, MARGIN = 2:4)

    M_qt   <- apply(X = M_it, FUN = quantile, MARGIN = 2,
                        probs = c(0.025, 0.5, .975 ) )

    mnM_t   <- apply(X = M_it, FUN = mean, MARGIN = 2)

    SSBpr_qspt   <- apply(X = SSBpr_ispt, FUN = quantile, MARGIN = 2:4,
                        probs = c(0.025, 0.5, .975 ) )

    if(proj)
      xLim <- range(yrs[(tMP-3):nTs[simIdx]])
    else xLim <- range(yrs[1:nTs[simIdx]])


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
      lines( x = yrs, y = mnSB_spt[1,1,], lwd = 3, lty = 2 )
      
      abline(v = yrs[tMP] - 0.5, lty = 2, lwd =.9 )
      abline(v = 2052, lty = 2, lwd =1 )
      
      abline(h = B0, lty = 2, col = "grey30" )
      # abline(h = TRP, lty = 2, col = "darkgreen" )
      # abline(h = USR, lty = 2, col = "darkgreen", lwd=2)
      # abline(h = Bcrash, lty = 5, col = "red", lwd=2)
      # abline(h = Bmsy, lty = 5, col = "purple", lwd=2)
      abline(h = LRP, lty = 2, col = "orange", lwd=2 )
      if(traces > 0)
        for( t in traceIdx )
          lines( x = yrs, y = SB_ispt[t,1,1,], lwd = 1, col = "black")

      if(mfg[2] == mfg[4]){
        axis(side = 4, las = 1)
        legend("topright", bty = "n", lty = c(2,2,5,2,5), col=c("grey30","darkgreen","purple","orange","red"), 
               c("B0","USR","Bmsy","LRP","Bcrash"), lwd = c(1,2,2,2,2))
      }
      
      # abline(h = 2.6, lty = 2, col = "red" )

    plot( x = xLim, y = yLimR,
          type = "n", axes = FALSE )
      mfg <- par("mfg")      

      if( mfg[2] == 1 )
      {
        axis(side = 2, las = 1)
        mtext( side = 2, text = "Rt (1e6)", line = 3 )
      }

      grid()
      box()

      polygon( x = c(yrs,rev(yrs)),
               y = c(R_qspt[1,1,1,],rev(R_qspt[3,1,1,])),
               col = "grey60", border = NA )
      lines( x = yrs, y = R_qspt[2,1,1,], lwd = 3 )      
      lines( x = yrs, y = mnR_spt[1,1,], lwd = 3, lty = 2 )
      
      abline(v = yrs[tMP] - 0.5, lty = 2, lwd = .9)


      if(traces > 0)
        for( t in traceIdx )
          lines( x = yrs, y = R_ispt[t,1,1,], lwd = 1, col = "black")

      abline(h = median(R0_isp), col = "purple", lty = 5, lwd=2)
      abline(h = median(bcR0_isp), col = "green", lty = 5, lwd=2)

      if( mfg[2] == mfg[4] ){
        legend("topright", 
                bty = "n",
                lty = 5, col=c("purple","green"), 
                legend = c("R0","bias Corrected R0"), 
                lwd=2)
        axis(side = 4, las = 1)
      }


    plot( x = xLim, y = yLimM,
          type = "n", axes = FALSE )
      mfg <- par("mfg")      

      if( mfg[2] == 1 )
      {
        axis(side = 2, las = 1)
        mtext( side = 2, text = "Natural Mortality (/yr)", line = 3 )
      }

     if( mfg[2] == mfg[4] ){
        legend("topright", 
                bty = "n",
                lty = 5, col=c("purple","green"), 
                legend = c("M0","bias Corrected M0"), 
                lwd=2)
        axis(side = 4, las = 1)
      }


      grid()
      box()
      polygon( x = c(yrs,rev(yrs)),
               y = c(M_qt[1,],rev(M_qt[3,])),
               col = "grey60", border = NA )
      lines( x = yrs, y = M_qt[2,], lwd = 3 )      
      lines( x = yrs, y = mnM_t, lwd = 3, lty = 2 )      
      if(simLabs[simIdx] != "NoFish"){
        abline( h = blobList[[simIdx]]$ctlList$mp$hcr$Uref_p[1], lty = 2, col = "grey30")
        abline(h = Ucrash, lty = 5, col = "red", lwd=2)
      }
      abline(v = yrs[tMP] - 0.5, lty = 2, lwd =.9 )

      if(traces > 0)
        for( t in traceIdx )
          lines( x = yrs, y = M_it[t,], lwd = 1, col = "black")


      abline(h = median(M0_isp), lty = 5, col = "purple", lwd=2)
      abline(h = median(bcM0_isp), col = "green", lty = 5, lwd=2)
      

    plot( x = xLim, y = yLimSBPR,
          type = "n", axes = FALSE )
      mfg <- par("mfg")      

      if( mfg[2] == 1 )
      {
        axis(side = 2, las = 1)
        mtext( side = 2, text = "SSB per recruit (kg)", line = 3 )
      }

      if( mfg[2] == mfg[4] ){
        legend("topright", 
                bty = "n",
                lty = 5, col=c("purple","green"), 
                legend = c("med phi0","mean phi0"), 
                lwd=2)
        axis(side = 4, las = 1)
      }
      axis(side = 1 )


      grid()
      box()
      polygon( x = c(yrs,rev(yrs)),
               y = c(SSBpr_qspt[1,1,1,],rev(SSBpr_qspt[3,1,1,])),
               col = "grey60", border = NA )
      lines( x = yrs, y = SSBpr_qspt[2,1,1,], lwd = 3 )      
      
      abline(v = yrs[tMP] - 0.5, lty = 2, lwd =.9 )
      

      if(traces > 0)
        for( t in traceIdx )
          lines( x = yrs, y = SSBpr_ispt[t,1,1,], lwd = 1, col = "black")

      abline(h =median(phi0_isp[,1,1]), lty = 2, lwd =2, col = "purple" )
      abline(h =mean(phi0_isp[,1,1]), lty = 2, lwd =2, col = "green" )
  }

  mtext(side = 1, text=  "Year", outer = TRUE, line = 2)
}
