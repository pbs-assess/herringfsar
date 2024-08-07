# Plotting functions for SISCAH TMB MP
# assessCAplots.R

# plotRule()
# Plot of Harvest control rule and recommended
# TAC based on SSPM estimate of biomass
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
  nT    <- length(C_t)
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


plotData <- function( ctlFile = "fitCtlFile.txt")
{
  # Read control file
  ctlPath <- file.path("./controlFiles",ctlFile)
  # read in control file
  ctlTable    <- .readParFile ( ctlPath )
  ctlList     <- .createList  ( ctlTable )

  # Read data
  dataCtl <- ctlList$data
  # Get model dimensions
  nA        <- dataCtl$nA
  minAge    <- dataCtl$minAge
  minAge_g  <- dataCtl$minAge_g
  fYear     <- dataCtl$fYearAssess
  lYear     <- dataCtl$lYearAssess
  fYearDat  <- dataCtl$fYearData
  lYearDat  <- dataCtl$lYearData
  yearsDat  <- fYearDat:lYearDat

  yearsAss  <- fYear:lYear
  tAssess   <- which(yearsDat %in% yearsAss)


  # save plots of data to data folder
  # Now load the fleetIDs
  loadStockSpecNameLists()


  dataList <- loadData( dataCtl$dataFolder,
                        ext = ".csv",
                        ctl = dataCtl )



  I_pgt     <- dataList$dataArrays$I_pgt[,,tAssess,drop = FALSE]
  C_pgt     <- dataList$dataArrays$C_pgt[,,tAssess,drop = FALSE]
  A_apgt    <- dataList$dataArrays$A_apgt[,,,tAssess,drop = FALSE]
  W_apgt    <- dataList$dataArrays$W_apgt[,,,tAssess,drop = FALSE]
  W_apt     <- dataList$dataArrays$W_apt[,,tAssess,drop = FALSE]
  mI_gt     <- dataList$dataArrays$mI_gt[,tAssess,drop = FALSE]
  mC_gt     <- dataList$dataArrays$mC_gt[,tAssess,drop = FALSE]
  mA_agt    <- dataList$dataArrays$mA_agt[,,tAssess,drop = FALSE]
  rI_pgt    <- dataList$dataArrays$rI_pgt[,,tAssess,drop = FALSE]
  combI_pt  <- dataList$dataArrays$combI_pt[,tAssess,drop = FALSE]

  nG <- dim(C_pgt)[2]
  nP <- dim(C_pgt)[1]

  gearLabs <- dimnames(C_pgt)[[2]]
  stockLabs <- dimnames(C_pgt)[[1]]

  datFolder <- file.path("./Data",dataCtl$dataFolder)
  plotFolder <- file.path(datFolder,"plots")
  if(!dir.exists(plotFolder))
    dir.create(plotFolder)
  graphics.off()

  # Want to plot
  # 1. Data summary
  summPlotFile <- file.path(plotFolder,"plotDataSummary.pdf")

  pdf(summPlotFile, width = 8.5, height = 11)
  plotStockDataSummary(dataList$dataArrays)
  dev.off()

  # 2. Index/catch series
  catIdxPlotFile <- file.path(plotFolder,"plotCatchIdx.pdf")
  pdf(catIdxPlotFile, width = 11, height = 8.5)
  plotDatIdxCatch(dataList$dataArrays)
  dev.off()

  # 3. Age data by year
  for( g in 1:nG )
    for( p in 1:nP )
    {
      if(any(A_apgt[2,p,g,] >= 0))
      {
        ageCompPlotFile <- file.path(plotFolder,paste0("plotAgeComps_",stockLabs[p],"_",gearLabs[g],".pdf"))
        pdf(ageCompPlotFile, width = 11, height = 8.5)
        plotDataCompYrs(  dataList = dataList$dataArrays,
                          pIdx = p,
                          gIdx = g,
                          fYear = 1951 )  
        dev.off()
      }
    }
  
  message("Data plots complete!\n")

}


# plotDatIdxCatch()
# Functions for plotting catch/index data
# on same axes, with smoothed indices.
plotDatIdxCatch <- function( dataList, fYear = 1951 )
{
  C_pgt     <- dataList$C_pgt
  I_pgt     <- dataList$I_pgt
  combI_pt  <- dataList$combI_pt

  # Get model dims
  nP <- dim(C_pgt)[1]
  nG <- dim(C_pgt)[2]
  nT <- dim(C_pgt)[3]

  yrs <- seq(from = fYear, by = 1, length.out = nT)

  I_pgt[I_pgt < 0] <- NA
  C_pt <- apply(X = C_pgt, FUN = sum, MARGIN = c(1,3) )
  
  # Fleet colours
  fleetCols <- RColorBrewer::brewer.pal(nG, "Dark2")
  

  par( mfrow = c(nP,1), mar = c(.1,.1,.1,.1), oma = c(3,3,1,1) )
  for( p in 1:nP )
  {
    maxY <- max(I_pgt[p,,], C_pt[p,], combI_pt[p,], na.rm = T)
    
    plot(x = range(yrs), y = c(0,maxY), type = "n", axes = FALSE )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis(side = 1)
      axis(side = 2, las = 1)
      grid()
      box()

      rect( xleft = yrs - 0.3, xright = yrs + 0.3,
            ybottom = 0, ytop = C_pt[p,], col = "salmon",
            border = NA )

      if( any(combI_pt > 0) )
      {

        # Plot blended index
        points(x = yrs, y = combI_pt[p,], pch = 16,
                cex = 1.2, bg = "grey60" )

        lines(lowess(y = combI_pt[p,], x = yrs, f = 0.1), 
          col = "grey40", lwd = 4)
      } # END blended spawn index

      if( all(combI_pt < 0))
      {
        # Plot split index
        for( g in 1:nG)
        {
          points(x = yrs, y = I_pgt[p,g,], pch = 16, col=fleetCols[g])
          if( any(!is.na(I_pgt[p,g,])))
          {
            whichYrs <- which(!is.na(I_pgt[p,g,]))
            lines(lowess(x = yrs[whichYrs],y = I_pgt[p,g,whichYrs], f = .3),
                  col = fleetCols[g], lwd = 4 )

          }
        }
        
      } # END split survey design
  } # END p loop
  
} # END plotDatIdxCatch()

# Plot comp fits
plotDataCompYrs <- function(  dataList = dataList$dataArrays,
                              comps = "age",
                              pIdx = 1,
                              gIdx = 1,
                              fYear = 1951 )
{
  
  obs_at    <- dataList$A_apgt[,pIdx,gIdx,]
  gearLabs  <- dimnames(dataList$A_apgt)[[3]]
  stockLabs <- dimnames(dataList$A_apgt)[[2]]

  # Pull model dims
  nT      <- dim(obs_at)[2]  
  nA      <- dim(obs_at)[1]
  nG      <- length(gearLabs)

  # Make colours vector
  cols    <- brewer.pal( n = nG, "Dark2" )

  # Make years vector
  years   <- seq(fYear, length.out = nT, by = 1)

  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  gearTimes <- which(obs_at[1,] >= 0)

  
  # Count the number of age observations
  # there are, and make the plotting window
  nObs <- length(gearTimes)
  nCols <- round(sqrt(nObs))
  nRows <- ceiling(nObs/nCols)

  par(  mfcol = c(nRows,nCols), 
        mar = c(0,0,0,0),
        oma = c(3,3,3,3) )

  obsProp_at <- obs_at
  nObs_t <- apply(X = obs_at, FUN = sum, MARGIN = 2)
  # Loop over times, normalise to proportions
  for(tIdx in gearTimes)
    obsProp_at[,tIdx] <- obs_at[,tIdx]/sum(obs_at[,tIdx])

  yMax <- max(obsProp_at, na.rm = T ) 

  for( tIdx in gearTimes )
  { 

    obsProp_a <- obsProp_at[,tIdx]  

    plot( x = c(1,nA), y = c(0,yMax ),
          xlab = "", ylab = "", type = "n", las = 1,
          axes = FALSE )
      box()
      grid()
      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )

      if( mfg[2] == 1 )
        axis( side = 2, las = 1 )

      rect( xleft = 1:nA - .3, xright = 1:nA + .3,
            ybottom = 0, ytop = obsProp_a[1:nA],
            col = "grey40", border = NA )
      panLab( x=.5, y = .95, txt = years[tIdx], font = 2 )
      legend( x = "topright",
              legend = paste("N = ", nObs_t[tIdx], sep = "" ),
              bty = "n" )

  }
    mtext( side = 3, outer = T, text = paste(stockLabs[pIdx], " ", gearLabs[gIdx]), line = 2, font = 2)
    mtext(  side = 1, outer = T, text = "Age", line = 2 )
    mtext(  side = 2, outer = T, text = "Proportion-at-age", 
            line = 2 )

} # END plotCompFitYrs()



# Plot multi-panel time series plots
plotMulti <- function(  ts = c("Nt","Bt","Ft"),
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
}



# Plot biomass
plotSBt <- function(  repList = reports,  
                      noPar = FALSE, 
                      pIdx = 1,
                      labcex = 1,
                      heading = NULL )
{
  report    <- repList$repOpt
  initYear  <- repList$fYear


  # Pull stuff from report
  SBt       <- report$SB_pt[pIdx,]
  Cgt       <- report$totC_pgt[pIdx,,]

  B0        <- signif(report$B0_p[pIdx],3)
  M         <- signif(report$M_p[pIdx],3)

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
    M  <- signif(mean(repList$posts$M_iapt[,1,pIdx,tInitModel_p[pIdx]+ 1]),3)

  } else {
    SB_qt <- 0
  }

  SBt[SBt == 0] <- NA

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

  
  # Get number of time steps
  nT    <- report$nT
  years <- seq(from = initYear, length = nT+1, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  depRange  <- range(SBt/B0, na.rm = T)
  maxDep    <- ceiling(depRange[2])

  depAxis   <- round(seq(0, maxDep, length.out = 5),3)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(5,4,3,1) )

  # Now plot recruitments
  plot( x = range(years), y = c(0,max(SBt,SB_qt, na.rm =T) ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    mfg <- par("mfg")
    # Add SSB axis
    axis(side = 2, las =1 )
    # Add depletion axis
    axis( side = 4, las = 1, labels = depAxis, at = B0 * depAxis)

    if( mfg[1] == mfg[3] )
      axis(side = 1 )

    box()
    grid()
    # rect( xleft = years[1:nT] - .3, xright = years[1:nT] + .3,
    #       ybottom = 0, ytop = Ct, border = NA, col = "grey60" )
    # Plot Biomass
    if( plotCI )
    {
      polyCol <- scales::alpha("red", .5)
      polygon(  x = c(years[1:nT],rev(years[1:nT])),
                y = c(SB_qt[1,1:nT], rev(SB_qt[3,1:nT])),
                col = polyCol, border = NA )
    }
    lines( x = years[1:nT], y = SBt[1:nT], lwd = 2, col = "red" )
    # abline( v = vertLines, lwd = .8, lty = 2, col = "grey80")
    abline( h = B0, lty = 3, lwd = .8, col = "red" )
    # points(x = years[(nT+1)], y = SBt[(nT+1)], col = "black", bg = "red", pch = 21)
    
    if(mfg[1] == mfg[3])
      mtext(side = 1, text = "Year", line = 3)

    if( mfg[2] == 1 )
      mtext(side = 2, text = "Spawning Biomass (kt)", line = 3.5, cex = labcex)
    if( mfg[2] == mfg[4] )
      mtext(side = 4, text = "Depletion", line = 2.5 )
    if( mfg[1] == 1 & !is.null(heading) )
      mtext(side = 3, text = stockLabel, line = 1, font = 2 )

    # panLab( x = 0.7, y = 0.8, txt = paste( "B0 = ",  B0, sep = "") )
    # panLab( x = 0.7, y = 0.7, txt = paste( "M = ",  M, sep = "") )
} # END plotSBt()

# Plot biomass
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

    
    text( x = years[nT - 7], y = c(.4,.3,.2,.1) * yMax,
          labels = c( paste( "B0 = ",  B0, sep = ""),
                      paste( "M0 = ", M0, sep = ""),
                      paste( "qs = ",  round(qg[4],2), sep = ""),
                      paste( "qd = ",  round(qg[5],2), sep = "")),

          cex = .9 )
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

# Plot age fits averaged over time
plotAgeFitAvg <- function(  repList = reports )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear
  gearLabels<- repList$gearLabs

  # Pull predicted and observed ages
  predAge <- report$predPA_apgt
  obsAge  <- report$A_apgt


  minAge_g <- repList$data$minAge_g

  gearLabs  <- dimnames(obsAge)[[3]]
  stockLabs <- dimnames(obsAge)[[2]]


  # replace missing entries with NAs
  predAge[predAge < 0] <- NA
  obsAge[obsAge < 0] <- NA


  # Get time steps
  nT <- dim(predAge)[4]
  nG <- dim(predAge)[3]
  nP <- dim(predAge)[2]

  # Multiply the yearly predicted proportions by the
  # yearly sample sizes
  for( g in 1:nG )
    for( p in 1:nP )
      for( t in 1:nT )
      {
        thisSamp <- sum( obsAge[,p,g,t], na.rm = T )
        if( thisSamp < 0 )
          thisSamp <- 0

        predAge[,p,g,t] <- thisSamp * predAge[,p,g,t]
      }


  
  # Average over time
  predAge <- apply( X = predAge, FUN = sum, MARGIN = c(1,2,3), na.rm = T )
  obsAge <- apply( X = obsAge, FUN = sum, MARGIN = c(1,2,3), na.rm = T )    

  # Pull model dims
  nG      <- report$nG
  nA      <- report$nA
  nT      <- report$nT
  nP      <- report$nP
  minPA   <- report$minPropAge

  # Make years vector
  years   <- seq(initYear, length = nT+1, by = 1)

  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  ageGears <- c()
  gearTimes <- vector(mode = "list", length = nG)
  for( gIdx in 1:nG )
    if( sum(obsAge[,,gIdx]) > 0  )
      ageGears <- c(ageGears,gIdx)

  # Make colours vector
  cols    <- brewer.pal( n = nG, "Dark2" )


  par(  mfcol = c(length(ageGears),nP), 
        mar = c(1,2,1,2),
        oma = c(3,3,3,3) )

  # ok, ageGears are the ones we want to plot,
  for( pIdx in 1:nP)
    for( aIdx in 1:length(ageGears) )
    { 
      gIdx <- ageGears[aIdx]
      # get age obs and preds
      ageObs  <- obsAge[,pIdx,gIdx]
      nSamp   <- sum(ageObs)
      ageObs  <- ageObs / sum(ageObs)
      agePred <- predAge[,pIdx,gIdx]
      agePred <- agePred / sum( agePred )

      plot( x = c(1,nA), y = c(0,1.1*max(ageObs,agePred,na.rm = T) ),
            xlab = "", ylab = "", type = "n", las = 1 )
        legend( x = "topright",
                bty = "n",
                legend = paste("N = ", nSamp, sep = "" ) )
        rect( xleft = minAge_g[gIdx]:nA - .3, xright = minAge_g[gIdx]:nA + .3,
              ybottom = 0, ytop = ageObs[minAge_g[gIdx]:nA],
              col = "grey40", border = NA )
        lines(  x = minAge_g[gIdx]:nA, y = agePred[minAge_g[gIdx]:nA], lwd = 3,
                col = cols[gIdx] )
        points(  x = minAge_g[gIdx]:nA, y = agePred[minAge_g[gIdx]:nA],
                col = cols[gIdx], pch = 16, cex = 1.5 )
        abline( h = minPA, lty = 2, lwd = .8 )
        
        
        # panLab( x=.5, y = .95, txt = gearLabels[gIdx] )
        mfg <- par("mfg")
        if(mfg[1] == 1 )
          mtext( side = 3, text = stockLabs[pIdx], line = 1, font = 2)

        if( mfg[2] == mfg[4] )
        {
          corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          par(xpd = TRUE) #Draw outside plot area
          text( x = corners[2]+0.5, 
                y = mean(corners[3:4]), 
                labels = gearLabs[gIdx], srt = 270,
                font = 2, cex = 1.5 )
          par(xpd = FALSE)
        }
    }
  
  mtext( side = 1, outer = T, text = "Age", line = 2 )
  mtext( side = 2, outer = T, text = "Proportion", line = 2 )

}

# Plot Biomass, catch and indices on the same
# axes
plotItResids <- function( repList =reports,
                          noPar = FALSE,
                          pIdx = 1, labcex = .8 )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear
  gearLabels<- repList$gearLabs
  # Pull stuff from report
  Bt     <- report$SB_pt[pIdx,]
  Igt    <- report$I_pgt[pIdx,,]
  qg     <- report$qhat_pg[pIdx,]
  rIgt   <- repList$data$rI_pgt[pIdx,,]

  # Need to pull vuln biomass/numbers
  vulnBgt <- report$vulnB_pgt[pIdx,,]
  vulnNgt <- report$vulnN_pgt[pIdx,,]

  # Pull probability of positive indices
  probPosIdx_gt <- report$probPosIdx_pgt[pIdx,,]

  # Get survey type and whether
  # idx is relative or absolute
  surveyType  <- report$survType_g      # abs or relative?
  indexType   <- report$indexType_g     # biomass or numbers
  calcIndex   <- report$calcIndex_g     # calc index at all?

  # Model dimensions
  nG          <- report$nG
  nT          <- report$nT
  surveyGears <- which(calcIndex == 1)

  # Get survey names
  surveyNames <- dimnames(vulnBgt)[[1]][surveyGears]
  stockLabel  <- dimnames(report$I_pgt)[[1]]

  if( stockLabel == "Aggregate" )
  {
    # We need to pull the aggregated stock
    # indices in at this point
    mIgt    <- report$mI_gt
    mqg     <- report$qhat_g

    mIgt[mIgt < 0] <- 0

    for( g in surveyGears )
    {
      if( all(Igt[g,] < 0) & any(mIgt[g,] > 0))
      {
        Igt[g,] <- mIgt[g,]
        qg[g]   <- mqg[g]
      }
    }
  }

  Igt[Igt<=0] <- NA

  if(!is.null(repList$posts))
  {
    tInitModel_p <- repList$repOpt$tInitModel_p
    SB_it <- repList$posts$SB_ipt[,pIdx,]

    Bt <- apply(X = SB_it, FUN = mean, MARGIN = c(2))

    q_ig  <- repList$posts$q_ipg[,pIdx,]
    qg   <- apply( X = q_ig, FUN = mean, MARGIN = 2)

    probPosIdx_igt <- repList$posts$probPosIdx_ipgt[,pIdx,,]
    probPosIdx_gt  <- apply( X = probPosIdx_igt, FUN = mean, MARGIN = c(2,3))

  } 


  # Scale indices by q
  IgtScaled <- Igt
  for( g in 1:nG)
  {
    IgtScaled[g,] <- Igt[g,] / qg[g] / probPosIdx_gt[g,]
  }

  IgtScaled[IgtScaled < 0] <- NA

  cols <- brewer.pal(nG, name = "Dark2")

  # Get number of time steps
  years <- seq(from = initYear, length = nT, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  tauObs_g <- report$tauObs_pg[pIdx,]

  resids <- IgtScaled

  AC_g <- rep(NA, nG)

  # Calc residuals
  for( g in 1:nG )
  {
    posIdx <- which(IgtScaled[g,] > 0)

    # Check index type, calcs resids
    if( surveyType[g] == 0 )
      resids[g,posIdx] <- log(IgtScaled[g,posIdx] / vulnBgt[g,posIdx] ) / tauObs_g[g]
    if( surveyType[g] == 1 )
      resids[g,posIdx] <- log(IgtScaled[g,posIdx] / vulnNgt[g,posIdx] ) / tauObs_g[g]
    if( surveyType[g] == 2 )
      resids[g,posIdx] <- log(IgtScaled[g,posIdx] / Bt[posIdx] ) / tauObs_g[g]

    if(length(posIdx))
      AC_g[g] <- acf(resids[g,posIdx], lag.max = 1, plot = FALSE, na.action = "na.pass")$acf[2]
  }   

  combAC <- c()

  combResids <- rep(NA,nT)
  if( any(repList$data$whichCombIdx_g == 1) )
    combResids <- -1 * report$zComb_pt[pIdx,]

  combAC <- acf(combResids, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2]


  # Now set up plotting window
  if( !noPar )
    par( mfrow = c(1, 1), mar = c(2,3,1,3), oma = c(3,3,1,1) )

  plot( x = range(years), y = range(-2,2,resids,combResids,na.rm=T), xlab = "",
          ylab = "", type = "n", las = 1 )
    mfg <- par("mfg")
    # show 0 line
    abline( h = 0, lty = 3, lwd = .8 )
    abline( v = vertLines, lwd = .8, lty = 2, col = "grey80")
    # Now show the resids
    for(g in surveyGears)
    {
      residTimes <- which(!is.na(resids[g,]))
      if(length(residTimes) > 0)
      {
        points( x = years[residTimes], y = resids[g,residTimes], pch = 16,
                col = cols[g] )
        # Fit a regression
        dat <- data.frame(x = years[residTimes], y = resids[g,residTimes])
        residModel <- lm( y~x, data = dat )
        dat$pred <- predict.lm(object = residModel, newdata = dat)

        pVal <- round(summary(residModel)$coefficients[2,4],2)

        lines( x = dat$x, y = dat$pred,
                col = cols[g], lwd = 2 )
        text( x = dat$x[4], y = min(dat$y) + .3, label = paste("p = ", pVal, sep = ""),
              col = cols[g], font = 2 )

        meanResid <- mean(resids[g,residTimes])
        segments(  x0 = years[residTimes[1]],
                  x1 = years[residTimes[length(residTimes)]], 
                  y0 = meanResid, col = cols[g], lty = 2, lwd = 2 ) 

        surveyNames[g] <- paste0(surveyNames[g],", m = ", round(meanResid,2),", AC = ", round(AC_g[g],2))
      }
    }
    combResids[combResids == 0] <- NA
    combResidTimes <- which(!is.na(combResids))

    if(length(combResidTimes) > 0)
    {
      points( x = years[combResidTimes], y = combResids[combResidTimes], pch = 16,
                col = "grey40" )
        # Fit a regression
        dat <- data.frame(x = years[combResidTimes], y = combResids[combResidTimes])
        residModel <- lm( y~x, data = dat )
        dat$pred <- predict.lm(object = residModel, newdata = dat)

        pVal <- round(summary(residModel)$coefficients[2,4],2)

        lines( x = dat$x, y = dat$pred,
                col = "grey40", lwd = 2 )
        text( x = dat$x[7], y = min(dat$y + .3), label = paste("p = ", pVal, sep = ""),
              col = "grey40", font = 2 )

        meanResid <- mean(combResids[combResidTimes])
        segments(  x0 = years[combResidTimes[1]],
                  x1 = years[combResidTimes[length(combResidTimes)]], 
                  y0 = meanResid, col = "grey40", lty = 2, lwd = 2 )

        surveyNames <- c(surveyNames,paste0("Spawn Index, m = ", round(meanResid,2),", AC = ", round(combAC,2)))
    }
    legend( x = "topleft", bty = "n",
            legend = surveyNames,
            pch = 16,
            col = c(cols[which(calcIndex == 1)],"grey40") )
    if( mfg[2] == 1)
      mtext(  side = 2, outer = F,
              text = paste("Std. log-residuals", sep = ""),
              line = 3, cex = labcex )
}



# plotStockDataSummary()
# A visual summary of the data
# available for a given stock
plotStockDataSummary <- function( dataList,
                                  fYear = 1951 )
{
  # Pull data
  I_pgt     <- dataList$I_pgt
  A_apgt    <- dataList$A_apgt
  C_pgt     <- dataList$C_pgt
  combI_pt  <- dataList$combI_pt
  # Mixed catch
  mC_gt  <- dataList$mC_gt

  # First, get dimension
  nT <- dim(I_pgt)[3]
  nP <- dim(I_pgt)[1]
  yrs <- seq( from = fYear, by = 1, length.out = nT+4)

  # Replace -1s with NA
  A_apgt[A_apgt < 0] <- NA
  I_pgt[I_pgt < 0] <- 0
  combI_pt[combI_pt < 0] <- NA

  gearCols <- RColorBrewer::brewer.pal(6, "Dark2")

  A_pgt <- apply( X = A_apgt, FUN = sum, MARGIN = c(2,3,4),
                  na.rm = TRUE )
  A_gt <- apply( X = A_pgt, FUN = sum, MARGIN = c(2,3),
                  na.rm = TRUE )
  sampSize_pg <- apply( X = A_pgt, FUN = sum, MARGIN = c(1,2) )
  sampSize_g  <- apply( X = A_gt, FUN = sum, MARGIN = c(1) )
  
  pchA_pgt <- A_pgt
  pchA_pgt[A_pgt < 200 ] <- 1
  pchA_pgt[A_pgt >= 200 ] <- 16
  pchA_pgt[A_pgt == 0 ] <- NA

  pchA_gt <- A_gt
  pchA_gt[A_gt < 200 ] <- 1
  pchA_gt[A_gt >= 200 ] <- 16
  pchA_gt[A_gt == 0 ] <- NA

  pchI_pgt <- I_pgt
  pchI_pgt[ I_pgt > 0 ] <- 16
  pchI_pgt[ I_pgt == 0 ] <- NA

  pchC_pgt <- C_pgt
  pchC_pgt[ C_pgt > 0 ] <- 16
  pchC_pgt[ C_pgt == 0 ] <- NA

  pchCombI_pt <- combI_pt
  pchCombI_pt[ combI_pt > 0] <- 16
  pchCombI_pt[ combI_pt == 0] <- 1


  I_gt <- apply( X = I_pgt, FUN = sum, MARGIN = c(2,3))
  pchI_gt <- I_gt
  pchI_gt[ I_gt > 0 ] <- 16
  pchI_gt[ I_gt == 0 ] <- NA  

  C_gt <- apply( X = C_pgt, FUN = sum, MARGIN = c(2,3))
  C_gt <- C_gt + mC_gt
  pchC_gt <- C_gt
  pchC_gt[ C_gt > 0 ] <- 16
  pchC_gt[ C_gt == 0 ] <- NA 

  combI_t <- apply( X = combI_pt, FUN = sum, MARGIN = c(2), na.rm = T)
  pchCombI_t <- combI_t
  pchCombI_t[combI_t > 0] <- 16
  pchCombI_t[combI_t == 0] <- NA


  ageGears <- 1:3
  catGears <- c(1:3)
  idxGears <- 4:5

  # Make a plotting area
  # Update this so it's responsive to data structures
  stockLab <- dimnames(combI_pt)[[1]]

  yLabs <- c("Ages", "Catch", "Spawn Idx" )

  gearLabs <- dimnames(C_pgt)[[2]]

  if( any(!is.na(combI_pt)))
    gearLabs <- c(gearLabs[1:3],"Spawn Index",gearLabs[-(1:5)])

  if(nP == 1)
    nPlots <- 1
  if(nP > 1)
    nPlots <- nP + 1

  par(  mfrow = c(nPlots,1), mar = c(0.5,1,.5,1),
        oma = c(7,6,3,4), xpd = FALSE )
  # Now plot aggregate/mixed data
  plot( x = range(yrs), y = c(0.5,3.5),
          type = "n", axes = FALSE )
    axis( side = 2, las = 1, at = 1:3,
          labels = yLabs)
    mtext( side = 3, text = "SISCAH Data Summary", font = 2)

    box()
    grid()  
    abline(h = c(1:2) + .5, lwd = 2)
    rmtext( txt = "Aggregate", cex = 1.5,
            line = .05, outer = TRUE, font = 2 )
    # Plot age comps
    points( x = yrs[1:nT], y = rep(1,nT)-.2,
            pch = pchA_gt[1,],
            col = gearCols[1] )
    points( x = yrs[1:nT], y = rep(1,nT),
            pch = pchA_gt[2,],
            col = gearCols[2] )
    points( x = yrs[1:nT], y = rep(1,nT)+.2,
            pch = pchA_gt[3,],
            col = gearCols[3] )
    # Plot Catch

    points( x = yrs[1:nT], y = rep(2,nT)-.3,
            pch = pchC_gt[1,],
            col = gearCols[1] )
    points( x = yrs[1:nT], y = rep(2,nT),
            pch = pchC_gt[2,],
            col = gearCols[2] )
    points( x = yrs[1:nT], y = rep(2,nT)+.3,
            pch = pchC_gt[3,],
            col = gearCols[3] )
    
    # Plot Indices
    points( x = yrs[1:nT], y = rep(3,nT)-.2,
            pch = pchI_gt[4,],
            col = gearCols[4] )
    points( x = yrs[1:nT], y = rep(3,nT)+.2,
            pch = pchI_gt[5,],
            col = gearCols[5] )
    points( x = yrs[1:nT], y = rep(3,nT),
            pch = pchCombI_t,
            col = "grey40" )



    text( x = yrs[nT]+4, y = 1 + c(-.2,0,.2), 
          labels = sampSize_g[1:3], cex = .8 )

  if(nP > 1)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs), y = c(0.5,3.5),
            type = "n", axes = FALSE )

        axis( side = 2, las = 1, at = 1:3,
              labels = yLabs)
        box()
        grid()
        rmtext( txt = stockLab[p], cex = 1.5,
                font = 2, line = .05, outer = TRUE )

        # Plot age comps
        points( x = yrs[1:nT], y = rep(1,nT)-.2,
                pch = pchA_pgt[p,1,],
                col = gearCols[1] )
        points( x = yrs[1:nT], y = rep(1,nT),
                pch = pchA_pgt[p,2,],
                col = gearCols[2] )
        points( x = yrs[1:nT], y = rep(1,nT) + .2,
                pch = pchA_pgt[p,3,],
                col = gearCols[3] )
        # Plot catch
        points( x = yrs[1:nT], y = rep(2,nT)-.3,
                pch = pchC_pgt[p,1,],
                col = gearCols[1] )
        points( x = yrs[1:nT], y = rep(2,nT)-.1,
                pch = pchC_pgt[p,2,],
                col = gearCols[2] )
        points( x = yrs[1:nT], y = rep(2,nT)+.1,
                pch = pchC_pgt[p,3,],
                col = gearCols[3] )
        points( x = yrs[1:nT], y = rep(2,nT)+.3,
                pch = pchC_pgt[p,6,],
                col = gearCols[6] )

        # Plot indices
        points( x = yrs[1:nT], y = rep(3,nT)-.2,
                pch = pchI_pgt[p,4,],
                col = gearCols[4] )
        points( x = yrs[1:nT], y = rep(3,nT)+.2,
                pch = pchI_pgt[p,5,],
                col = gearCols[5] )
        points( x = yrs[1:nT], y = rep(3,nT),
                pch = pchCombI_pt[p,],
                col = "grey40" )
        abline(h = c(1:2)+.5, lwd = 2)

    
        text( x = yrs[nT]+4, y = 1 + c(-.2,0,.2), 
              labels = sampSize_pg[p,1:3], cex = .8 )

    }
  }
  axis( side = 1)
  mtext( side = 1, text = "Year", font = 2, outer = TRUE,
          line = 2)

    if( any(!is.na(combI_pt)))
      gearCols <- c(gearCols[1:3],"grey40",gearCols[6])

    par(xpd=TRUE, oma = c(0,0,0,0))
    # plot( x = c(0,1), y = c(0,1), add = TRUE)
    legend( x = "bottom", bty = "n",
            # inset = c(0,-1), xpd = TRUE,
            legend = gearLabs, 
            pch = 16, col = gearCols,
            horiz = TRUE, cex = 1.5, pt.cex = 2)
}



# panLab      (Place text labels in plot region)
# Purpose:    Place a text label in the plot region defined by (0,1), (0,1).
#             The ... notation allows all parameters available to "text" to be
#             passed.
# Parameters: x, y are the coordinates of the label
#             txt is the text
# Returns:    NULL (invisibly)
# Source:     A.R. Kronlund
# Revised:    K.Holt; 13-Jan-10 to accomodate axes on log scale
panLab <- function( x, y, txt, ... )
{
  # Allows text to be placed in plot panel at 0<x<1, 0<y<1.
  usr <- par( "usr" )
  
  yLog <- par("ylog")
  xLog <- par("xlog")
  
  # Check for log-transformed axes and adjust usr commands as needed
  # note: when a log scale is in use, 
  #           usr gives limits in the form 10 ^ par("usr")

  # Case 1: neither axis is on the log scale
  if (yLog==FALSE & xLog==FALSE)
  {
    par( usr=c(0,1,0,1) )
  }
  # Case 2: only the y-axis is on log scale
  if (yLog==TRUE & xLog==FALSE) 
  {
    usr[3:4]<-10 ^ par("usr")[3:4]
    par( usr=c(0,1,0,1), ylog=FALSE )
  } 
  # Case 3: only the x-axis is on log scale
  if (yLog==FALSE & yLog==TRUE) 
  {
    usr[1:2]<-10 ^ par("usr")[1:2]
    par( usr=c(0,1,0,1), xlog=FALSE )
  } 
  # Case 4: both axes are on the log scale
  if (yLog==TRUE & xLog==TRUE) 
  {
    usr[1:4]<-10 ^ par("usr")[1:4]
    par( usr=c(0,1,0,1), xlog=FALSE, ylog=FALSE )
  } 
  text( x, y, txt, ... )
  par( usr=usr )

}