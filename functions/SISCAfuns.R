# -----------------------------------------------------------------------------
#
# SISCAfuns.R
#
# Spatially Integrated Statistical Catch-at-Age assessment
# model fitting functions.
# 
# 
# -----------------------------------------------------------------------------

Sys.unsetenv("PKG_CXXFLAGS")

# fitSISCA()
# Wrapper function to fit the assessCA model under a data scenario 
# and model hypothesis, both defined in the ctlFile. If fit is succesful
# then model outputs are saved to a folder in the ./Outputs/fits/ directory
# inputs:   ctlFile=character with the name/path of the control file
#           folder=optional character name of output folder 
#                   ie saves to ./Outputs/<folder>
# ouputs:   NULL
# usage:    from the console to run the procedure
fitSISCA <- function (  ctlFile = "fitCtlFile.txt", 
                        ctlDir  = "controlFiles",
                        folder=NULL, quiet=FALSE )
{ 
  ctlPath     <- file.path(".",ctlDir,ctlFile)
  baseCtlPath <- file.path(".",ctlDir,"baseCtlFile.txt")
  
  # read in base control file
  baseCtlTable <- .readParFile( baseCtlPath )
  baseCtlList  <- .createList( baseCtlTable )

  # read in control file
  ctlTable    <- .readParFile( ctlPath )
  ctlList     <- .createList( ctlTable )
  
  # Now can we overwrite elements... just the data elements
  for(i in 1:length(ctlList$data))
  {
    parName                     <- names(ctlList$data)[i]
    baseCtlList$data[[parName]] <- ctlList$data[[parName]]
  }

  for(i in 1:length(ctlList$ctrl))
  {
    parName                     <- names(ctlList$ctrl)[i]
    baseCtlList$ctrl[[parName]] <- ctlList$ctrl[[parName]]
  }

  baseCtlList$mp <- ctlList$mp

  # Run simEst Procedure
  reports <- .runSISCA( obj = baseCtlList )

  # save output to project folder
  # First, if a folder name isn't nominated, create a default sim folder
  if ( is.null(folder) )
  {
    stamp <- paste( format(Sys.time(),format="%Y%m%d%H%M%S" ),sep="" )
    folderName <- paste ( "fit_",stamp, sep = "" )
    folder <- stamp
  } else folderName <- paste ("fit_", folder, sep = "")

  # Now paste together the path to the folder and create it
  path <- file.path (getwd(),"Outputs",folderName)
  dir.create ( path )
  message( "\nMSG (saveSim) Created assessment fit folder ./Outputs/", folderName, "\n" )
  # Save blob
  # Add folder to path
  reports$folder <- path
  reports$ctlList <- baseCtlList

  fitTableList <- makeParFitTableMLE( reports )

  reports$parTable      <- fitTableList$parEstTable


  # Save report
  saveRDS(reports, file = file.path(path,paste(folderName,".rds",sep="")))
  
  # Save grad table
  write.csv(reports$grad, file = file.path(path,"gradRpt.csv"))

  # Copy cpp to sim folder for posterity
  file.copy(from="SISCA.cpp",to=file.path(path,"SISCA.cpp"))

  # Copy control file to sim folder for posterity
  cat(  "# fitCtlFile.txt, written to ", folder, " on ", as.character(Sys.Date()),"\n", sep = "", 
        file = file.path(path,"fitCtlFile.txt"))
  write.table(  ctlTable,
                file = file.path(path,"fitCtlFile.txt"),
                row.names = FALSE,
                quote = FALSE, qmethod = "double",
                append = TRUE )
  cat(  "# <End File>", sep = "", 
        file = file.path(path,"fitCtlFile.txt"),
        append = TRUE)
  
  makeMPReport( fitID = folder )

  beepr::beep()
  # Done
} # END fitSISCA()



# .runSISCA()
# Procedure to create data, parameter and map 
# lists, and fit the assessCA TMB model. There are
# some path specific features here that relate to YE
# assessment stuff, so take care if copying to other
# locations
.runSISCA <- function( obj = controlList )
{
  # Get data scenario and model hypothesis control lists
  dataCtl <- obj$data
  hypoCtl <- obj$hypo
  ctrlObj <- obj$ctrl
  phases  <- obj$phases

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

  yrDatChar   <- as.character(yearsDat)
  yrAssChar   <- as.character(yearsAss)

  stockNames  <- dataCtl$stock
  gearNames   <- dataCtl$fleetnames_g
  ageNames    <- as.character(1:nA)

  nP          <- length( stockNames )
  nG          <- length( gearNames )
  nT          <- length( yearsAss )

  initModelYear_p <- hypoCtl$initModelYear_p

  # Now load the fleetIDs
  loadStockSpecNameLists()

  dataList <- loadData( dataCtl$dataFolder,
                        ext = ".csv",
                        ctl = dataCtl )

  I_pgt     <- dataList$dataArrays$I_pgt[,,tAssess,drop = FALSE]
  E_pgt     <- dataList$dataArrays$E_pgt[,,tAssess,drop = FALSE]
  C_pgt     <- dataList$dataArrays$C_pgt[,,tAssess,drop = FALSE]
  A_apgt    <- dataList$dataArrays$A_apgt[,,,tAssess,drop = FALSE]
  W_apgt    <- dataList$dataArrays$W_apgt[,,,tAssess,drop = FALSE]
  W_apt     <- dataList$dataArrays$W_apt[,,tAssess,drop = FALSE]
  mI_gt     <- dataList$dataArrays$mI_gt[,tAssess,drop = FALSE]
  mC_gt     <- dataList$dataArrays$mC_gt[,tAssess,drop = FALSE]
  mA_agt    <- dataList$dataArrays$mA_agt[,,tAssess,drop = FALSE]
  rI_pgt    <- dataList$dataArrays$rI_pgt[,,tAssess,drop = FALSE]
  combI_pt  <- dataList$dataArrays$combI_pt[,tAssess,drop = FALSE]


  whichCombIdx_g <- rep(0,nG)
  if( dataCtl$combSpawnIdx )
  {
    whichCombIdx_g[4:5] <- 1
  }

  if( hypoCtl$backSplitSOK )
  {
    if( nP == 1 )
      meanSOKprop_p <- C_pgt[1,6, ]
    else
      meanSOKprop_p <- apply(X = C_pgt[ ,6, ], FUN = sum, MARGIN = 1 )
    meanSOKprop_p <- meanSOKprop_p / sum(meanSOKprop_p)

    for( p in 1:nP )
      C_pgt[p,6,26:44] <- mC_gt[6,26:44] * meanSOKprop_p[p]

    mC_gt[6,] <- 0   

  }


  # Calculate first time step from first year
  tInitModelYear_p  <- initModelYear_p - fYear
  tInitModelYear_p[tInitModelYear_p < 0] <- 0

  # remove any data for populations that have
  # initModelYear after fYear
  for( p in 1:nP )
  {
    if( initModelYear_p[p] > fYear)
    {
      C_pgt[ p,,1:(tInitModelYear_p[p] - 1)] <- 0
      I_pgt[ p,,1:(tInitModelYear_p[p] - 1)] <- -1
      A_apgt[ ,p,,1:(tInitModelYear_p[p] - 1)] <- -1
      W_apgt[ ,p,,1:(tInitModelYear_p[p] - 1)] <- -1
    }
  } 

  # Sum catch for an initial B0 value
  sumCat <- apply( X = C_pgt, FUN = sum, MARGIN = c(1), na.rm = T )
  

  # Figure out years for propEff estimation
  sokTimes <- c()
  sokFleets <- which(dataCtl$fleetType[gearNames] %in% c(2,3))

  if( length(sokFleets) > 0)
  {
    for( g in sokFleets)
    {
      if( g > nG)
        next
      # First do stock specific
      for( p in 1:nP )
      {
        whichTimes <- which( C_pgt[p,g,,drop = FALSE] > 0 )
        sokTimes <- union(sokTimes, whichTimes)
      }
      # Then do mixed
      whichTimes <- which( mC_gt[g,] > 0)
      sokTimes <- union(sokTimes,whichTimes)
    }
  }

  sokOn <- rep(0,nT)
  sokOn[sokTimes] <- 1

  # Calculate the initial steepness value
  initSteep <- hypoCtl$initSteep

  # Then convert to logit scale (with adjustments 
  # for bounds) for the inital value - good for fixing in 
  # estimation
  initLogitSteep <- -log( 0.78 / (initSteep - 0.2) - 1 )



  # Generate the model switch settings we're using
  fleetSwitches <- setFleetSwitches()
  
  calcIndex <- integer(length = nG)
  names(calcIndex) <- gearNames
  for( g in 1:nG )
    if( any(I_pgt[,gearNames[g],] > 0) | 
        any(mI_gt[gearNames[g],] > 0) )
      calcIndex[g] <- 1

  # Initialised in a fished state (0 == unfished, 1 == fished)
  initFished    <- hypoCtl$initFished_p[1:nP]

  if( hypoCtl$fishedInitMethod == "surv")
    nFishedInitAges <- nA - 1

  if( hypoCtl$fishedInitMethod == "nums")
    nFishedInitAges <- nA 

  # Recruitment deviations - need to add 
  # different initial and terminal rec devs
  # for each stock
  yFirstRecDev_p  <- hypoCtl$yFirstRecDev_p
  yLastRecDev_p   <- hypoCtl$yLastRecDev_p

  # Initial conditions
  initFcode_p     <- hypoCtl$initFcode_p
  initRcode_p     <- hypoCtl$initRcode_p
  initDevCode_p   <- hypoCtl$initDevCode_p
  avgRcode_p      <- hypoCtl$avgRcode_p
  initMdev_p      <- hypoCtl$initMdev_p

  if( all(initFished == 0) )
    phases$log_initN_ap <- -1

  nInitDevs <- sum(initFished) * nA 

  # Generate tFirstRecDev from yFirstRecDev and fYear
  tFirstRecDev_p    <- yFirstRecDev_p - fYear
  tLastRecDev_p     <- yLastRecDev_p - fYear

  # Adjust so that rec devs are inside 1:(nT - minAge - 1 )
  tFirstRecDev_p[ tFirstRecDev_p < 1]       <- 1
  tLastRecDev_p[  tLastRecDev_p > (nT-2) ] <- nT-2

  nRecDevs_p      <- tLastRecDev_p - tFirstRecDev_p

  # Turn off multi-stock parameters if nP == 1
  if( nP == 1 )
  {
    phases$epsSelAlpha_pg   <- -1
    phases$epsSelBeta_pg    <- -1
    phases$epsM_p           <- -1
    phases$epsSteep_p       <- -1
  }

  # Now make a vector to switch time varying selectivity
  # on and off
  tvSelFleets <- rep(0,nG)
  names(tvSelFleets) <- gearNames
  tvSelFleets[ hypoCtl$tvSelFleets == 1 ] <- 1 

  for( g in 1:nG )
    for( p in 1:nP )
      for( t in 1:nT )
      {
        if( A_apgt[1,p,g,t] >= 0)
          if( sum(A_apgt[,p,g,t]) < dataCtl$minCompSampSize )
            A_apgt[,p,g,t] <- -1
      }

  # Calculate mean sample sizes for compositional data
  # First, stock specific
  meanA_apgt <- A_apgt
  meanA_apgt[meanA_apgt < 0] <- NA
  meanA_pgt <- array(NA, dim = c(nP,nG,nT))
  meanA_pgt[1:nP,,] <- apply( X = meanA_apgt, FUN = sum, MARGIN = c(2,3,4), na.rm = T)
  meanA_pgt[meanA_pgt == 0] <- NA
  meanA_pg <- array(NA, dim = c(nP,nG))
  meanA_pg[1:nP,] <- apply( X = meanA_pgt, FUN = mean, MARGIN = c(1,2), na.rm = T )

  meanA_pg[is.nan(meanA_pg)] <- 0

  # Now mixed
  meanA_agt <- mA_agt
  meanA_agt[meanA_agt < 0] <- NA
  meanA_gt <- apply( X = meanA_agt, FUN = sum, MARGIN = c(2,3), na.rm = T)
  meanA_gt[meanA_gt == 0] <- NA
  meanA_g <- apply( X = meanA_gt, FUN = mean, MARGIN = c(1), na.rm = T ) 
  

  meanA_g[is.nan(meanA_g)] <- 0

  # Juvenile M source
  if( hypoCtl$juveMsource == "est" )
    juveMsource <- 1

  if( hypoCtl$juveMsource == "mean" )
    juveMsource <- 0

  # Calculate the number of selectivity deviations
  nSelDevs <- 0
  for( gIdx in which(hypoCtl$tvSelFleets == 1) )
    for( pIdx in 1:nP)
      if( any(A_apgt[1,pIdx,gIdx,yrAssChar] >= 0) )
        nSelDevs <- nSelDevs + length(which(A_apgt[1,pIdx,gIdx,yrAssChar] >= 0))


  idxWithZeroes_pg <- array(0,dim = c(nP,nG))
  deltaIdx_g <- dataCtl$deltaIdx[1:nG]

  if( dataCtl$insertZeroes )
  {
    
    # Loop through gears and stocks and replace missing data
    # with zeroes
    if( nP > 1)
    { 
      # Surface survey
      

      # Dive survey
      diveIdx <- (1988 - fYear + 1):nT
      for( p in 1:nP )
      {
        whichMissing <- which(I_pgt[p,5,diveIdx] < 0)
        I_pgt[p,5,diveIdx[whichMissing]] <- 0
        idxWithZeroes_pg[p,5] <- 1
      }

      
    }
  } else {
    combI_pt[combI_pt == 0] <- -1
  }


  deltaIdx_pg <- idxWithZeroes_pg
  deltaIdx_pg[,deltaIdx_g <= 0 ] <- 0

  deltaIndVar <- 0
  if( hypoCtl$deltaIndVar == "Biomass" )
    deltaIndVar <- 1
  if( hypoCtl$deltaIndVar == "Depletion" )
    deltaIndVar <- 2

  SRindVar <- 2
  if( hypoCtl$SRindVar == "Biomass" )
    SRindVar <- 1
  if( hypoCtl$SRindVar == "Eggs" )
    SRindVar <- 2  

  scaleSel_gt <- array(1, dim = c(nG,nT), dimnames = list(  gear = gearNames, 
                                                            year = yearsAss))
 
  diluteN_apt       <- array(0,dim =c(nA,nP,nT))  
  minAgeDilute      <- nA-1
  alphaPropDilute   <- 0
  betaPropOverlap_g <- rep(0,nG)
  names(betaPropOverlap_g) <- gearNames

  C_t <- apply(X = C_pgt, FUN = sum, MARGIN = 3)
  posCatIdx <- which(C_t > 0)
  catAllocYrs <- c(max(posCatIdx)-9,max(posCatIdx))

  omegaRradius <- 10
  if(!is.null(hypoCtl$omegaRradius))
    omegaRradius <- hypoCtl$omegaRradius
  
  # Generate the data list
  data <- list( # Input data
                I_pgt           = I_pgt,
                C_pgt           = C_pgt,
                E_pgt           = E_pgt,
                A_apgt          = A_apgt,
                W_apgt          = W_apgt,
                W_apt           = W_apt,
                mI_gt           = mI_gt,
                mC_gt           = mC_gt,
                mA_agt          = mA_agt,
                rI_pgt          = rI_pgt,
                combI_pt        = combI_pt,
                diluteN_apt     = diluteN_apt,
                # Model switches
                survType_g      = dataCtl$survType,
                indexType_g     = dataCtl$idxType,
                deltaIdx_pg     = deltaIdx_pg,
                deltaIndVar     = deltaIndVar,
                qPrior_g        = hypoCtl$qPrior_g[gearNames],
                calcIndex_g     = calcIndex[gearNames],
                selType_g       = hypoCtl$selFun[gearNames],
                scaleSel_gt     = scaleSel_gt,
                hierSel         = as.integer(hypoCtl$hierSel),
                condMLEtauObs   = hypoCtl$condMLEtauObs,
                ageCompWeight_g = hypoCtl$ageCompWt[gearNames],
                idxLikeWeight_g = hypoCtl$idxLikeWt[gearNames],
                catLikeWeight_g = hypoCtl$catLikeWt[gearNames],
                fleetTiming_g   = dataCtl$fleetTiming_g[gearNames],
                fleetType_g     = as.integer(dataCtl$fleetType[gearNames]),
                mortType_g      = as.integer(dataCtl$mortType[gearNames]),
                catSeriesType_g = as.integer(dataCtl$catSeriesType[gearNames]),
                # Dilute catches with another mixing stock
                minAgeDilute    = minAgeDilute,
                alphaDilute     = alphaPropDilute,
                betaDilute_g    = betaPropOverlap_g[gearNames],
                # Selectivity
                selX_g          = hypoCtl$fleetSelX[gearNames],
                tvSelFleets     = tvSelFleets,
                # Initial conditions
                initCode_p      = c(initFished),
                initFcode_p     = initFcode_p,
                initRcode_p     = initRcode_p,
                initMethod      = hypoCtl$fishedInitMethod,
                tInitModel_p    = tInitModelYear_p,
                # posfun
                posPenFactor    = c(dataCtl$posPenFactor),
                # Rec devs
                firstRecDev_p   = tFirstRecDev_p,
                lastRecDev_p    = tLastRecDev_p,
                omegaRradius    = omegaRradius,
                # Age likelihood controls
                minPropAge      = dataCtl$minPropAge,
                minAge_g        = as.integer(dataCtl$minAge_g),
                # year range for averaging M/weight for projections
                nYearsProjAve   = as.integer(5),
                # Timing
                spawnTiming     = hypoCtl$spawnTiming,
                moveTiming      = hypoCtl$moveTiming,
                fec             = hypoCtl$SOKfec,
                gamma_g         = hypoCtl$SOKgamma_g,
                pFem            = hypoCtl$SOKpFem,
                postPondM_g     = hypoCtl$SOKpostPondM_g,
                sokInitF        = hypoCtl$SOKInitF,
                useMovement     = hypoCtl$useMovement,
                juveMage        = hypoCtl$juveMage,
                juveMsource     = juveMsource,
                calcPropEff_t   = sokOn,
                jeffWtB0        = hypoCtl$jeffWtB0,
                lnm1PriorWt     = hypoCtl$lnm1PriorWt,
                meanRdevWt      = hypoCtl$meanRdevWt,
                compLikeFun     = hypoCtl$ageCompLik,
                meanSampSize_pg = meanA_pg,
                meanSampSize_g  = meanA_g,
                avgRcode_p      = avgRcode_p,
                SRindVar        = SRindVar,
                whichCombIdx_g  = whichCombIdx_g,
                densityDepM     = as.integer(hypoCtl$densityDepM),
                corrMdevs       = as.integer(hypoCtl$corrMortDevs),
                corrRdevs       = as.integer(hypoCtl$corrRecDevs),
                corrParWeight   = hypoCtl$corrParWeight,
                mixComps_g      = rep(0,nG),
                catAllocYrs     = catAllocYrs,
                C_v             = 0,
                dataLikeWt      = hypoCtl$dataLikeWt,
                priorDensWt     = hypoCtl$priorDensWt )


  # pull some pars from the fitCtlFile to
  # calculate some initial parameter values
  initSelAlpha  <- hypoCtl$fleetSelAlpha[gearNames]
  initSelBeta   <- hypoCtl$fleetSelBeta[gearNames]

  initdSelAlpha <- hypoCtl$fleetdSelAlpha[gearNames]
  initdSelBeta  <- hypoCtl$fleetdSelBeta[gearNames]

  aMat50 <- hypoCtl$matPars[1]
  aMat95 <- hypoCtl$matPars[2]

  initTau2_g <- hypoCtl$tau2ObsPriorMode[1:nG]

  initPhi1  <- rep(0,nG)
  initPsi   <- rep(0,nG)

  if( hypoCtl$ageCompLik == 2 )
  {
    initPhi1  <- rep(0, nG)
    initPsi   <- rep(-5, nG)
  }


  #### MOVEMENT MATRIX ####
  # Really currently a placeholder.
  mov_ppa <- array( 0, dim = c(nP,nP,nA) )

  if( hypoCtl$moveMatType == "Identity" )
  {
    for( a in 1:nA )
      mov_ppa[,,a] <- diag(1,nP)
  }


  # Make initial selectivity devs array
  initepsSelAlpha_pg <- array(0, dim = c(nP,nG))
  initepsSelBeta_pg  <- array(0, dim = c(nP,nG))

  initepsdSelAlpha_pg <- array(0, dim = c(nP,nG))
  initepsdSelBeta_pg  <- array(0, dim = c(nP,nG))

  if( hypoCtl$hierSel == 0 )
  {
    for( p in 1:nP )
    {
      initepsSelAlpha_pg[p,]  <- log(initSelAlpha)
      initepsSelBeta_pg[p,]   <- log(initSelBeta)

      initepsdSelAlpha_pg[p,]  <- log(initdSelAlpha)
      initepsdSelBeta_pg[p,]   <- log(initdSelBeta)
    }
  }


  initTau2_pg <- array( 0, dim =c(nP,nG) )

  for( p in 1:nP )
    initTau2_pg[ p, ] <- initTau2_g



  # Update sdq if using shrinkage on any survey
  sdq <- hypoCtl$sdq
  sdq[hypoCtl$qPrior_g == 2] <- .5
  if(hypoCtl$qPrior_g[5] == 2)
    sdq[5] <- .01

  lntauObsComb_pg  <- array(0, dim = c(nP,nG))

  lntauObsComb_pg[,4] <- log(sqrt(1.166*initTau2_g[5]))  
  lntauObsComb_pg[,5] <- log(sqrt(initTau2_g[5]))  

  initM <- hypoCtl$Mprior[1]
  if( hypoCtl$densityDepM )
  {
    initM <- hypoCtl$ddBaseM
  }

  # Get initial lnqF_g parameter
  lnqF_g  <- hypoCtl$mlnqF_g[gearNames]
  lnqFinit_pg <- array(0,dim = c(nP,nG))
  
  for(p in 1:nP )
  {
    lnqFinit_pg[p,] <- hypoCtl$mlnqF_g
  }
    

  # Generate parameter list
  pars <- list( lnB0_p                = log(5*sumCat),
                lnRinit_p             = rep(10,nP),
                lnRbar_p              = rep(10,nP),
                logit_ySteepness      = initLogitSteep,
                lnM                   = log(initM),
                lnMjuve               = log(hypoCtl$Mjuve),
                lnm1                  = log(hypoCtl$m1Prior[1]),
                epslnm1_p             = rep(0,nP),
                fDevs_ap              = matrix( data=0, nrow=nFishedInitAges, ncol=nP ),
                lnFinit_p             = rep(-4,nP),
                epsM_p                = rep(0,nP),
                lnsigmaStockM         = log(hypoCtl$sigmaMStock),
                epsSteep_p            = rep(0,nP),
                lnsigmaStockSteep     = log(hypoCtl$sigmaSteepStock),
                
                # Selectivity (up)
                lnSelAlpha_g          = log(initSelAlpha),
                lnSelBeta_g           = log(initSelBeta),
                epsSelAlpha_pg        = initepsSelAlpha_pg,
                epsSelBeta_pg         = initepsSelBeta_pg,
                epsSelAlpha_vec       = rep(0,nSelDevs),
                epsSelBeta_vec        = rep(0,nSelDevs),
                
                # Selectivity (dn)
                lndSelAlpha_g         = log(initdSelAlpha),
                lndSelBeta_g          = log(initdSelBeta),
                epsdSelAlpha_pg       = initepsdSelAlpha_pg,
                epsdSelBeta_pg        = initepsdSelBeta_pg,
                epsdSelAlpha_vec      = rep(0,nSelDevs),
                epsdSelBeta_vec       = rep(0,nSelDevs),
                
                # Selectivity dev SDs
                lnsigmaSelAlpha_g     = log(hypoCtl$sigmaSelStock_g),
                lnsigmaSelBeta_g      = log(hypoCtl$sigmaSelStock_g),
                lnsigmaTVsel          = log(hypoCtl$tvSelSD),
                
                # FPUE/per cap Pred Mort
                lnqFinit_pg           = lnqFinit_pg,
                deltaqF_pgt           = array(0,dim = c(nP,nG,nT)),
                resCatCV_g            = hypoCtl$catchCV[gearNames],
                
                # obs error
                lntau2Obs_pg          = log(initTau2_pg),
                lntau2Obs_g           = log(initTau2_g),
                lntauObsComb_pg       = lntauObsComb_pg,
                # Rec devs
                recDevs_pt            = array(0, dim = c(nP,nT-1) ),
                lnsigmaR              = log(hypoCtl$sigmaR),
                priorSigR             = hypoCtl$priorSigR,
                # M devs
                omegaM_pt             = array(0, dim = c(nP,nT)),
                lnsigmaM              = log(hypoCtl$sigmaM),
                # Correlation in proc errors
                off_diag_R            = rep(0, nP * (nP-1) / 2),
                off_diag_M            = rep(0, nP * (nP-1) / 2),
                # Priors
                obstau2IGa            = rep(hypoCtl$tau2ObsIGa[1],nG),
                obstau2IGb            = (hypoCtl$tau2ObsIGa[1]+1)*hypoCtl$tau2ObsPriorMode,
                sig2RPrior            = c(1,2),
                sig2MPrior            = c(1,0.04),
                rSteepBetaPrior       = hypoCtl$steepnessPrior,
                Mprior                = hypoCtl$Mprior,
                m1Prior               = hypoCtl$m1Prior,
                mlnq_g                = rep(0,nG),
                sdlnq_g               = hypoCtl$sdq,
                mq                    = hypoCtl$mq,
                sdq                   = sdq,
                lnqComb_pg            = array(0, dim = c(nP,nG)),
                mlnqF_g               = hypoCtl$mlnqF_g[gearNames],
                sdlnqF_g              = hypoCtl$sdlnqF_g[gearNames],
                sddeltaqF_g           = hypoCtl$sddeltaqF_g[gearNames],
                # revise maturity ages and growth model - most don't matter since
                # the WAA is empirical
                mat_a                 = hypoCtl$matVec,
                fec_a                 = rep(200,nA),
                Linf                  = hypoCtl$vonLinf,
                L1                    = hypoCtl$vonL1,
                vonK                  = hypoCtl$vonK,
                inputL1               = hypoCtl$inputL1,
                lenWt                 = hypoCtl$alloLW,
                mlnSelAlpha_g         = log(hypoCtl$fleetSelAlpha),
                mlnSelBeta_g          = log(hypoCtl$fleetSelBeta),
                sdSel_g               = hypoCtl$sdSel_g,
                selPriorWt_g          = hypoCtl$selPriorWt_g[gearNames],
                # Movement
                mov_ppa               = mov_ppa,
                # SOK model
                propEffBounds         = hypoCtl$sokPropEffBounds,
                logitPropEff_vec      = rep(0,length(sokTimes)),
                mPsi                  = .06,
                sdPsi                 = .1,
                logitphi1_g           = initPhi1,
                logitpsi_g            = initPsi,
                lnSDProbPosIdx_pg     = array(0,dim = c(nP,nG)),
                meanProbPosIdx_pg     = array(0,dim = c(nP,nG)),
                muSDProbPosIdx_g      = rep(hypoCtl$priorDeltaLNsd[1],nG),
                muMeanProbPosIdx_g    = rep(hypoCtl$priorDeltaLNmean[1],nG),
                sigSDProbPosIdx_g     = rep(hypoCtl$priorDeltaLNsd[2],nG),
                sigMeanProbPosIdx_g   = rep(hypoCtl$priorDeltaLNmean[2],nG) )

  # Set phases for correlation matrix pars to negative
  # if not required
  if( hypoCtl$ageCompLik == 0 )
  {
    phases$logitphi1_g  <- -1
    phases$logitpsi_g   <- -1
  }

  # Set phases for correlation matrix pars to negative
  # if not required
  if( hypoCtl$ageCompLik == 1 )
  {
    phases$logitpsi_g   <- -1
  }

  if( hypoCtl$condMLEtauObs )
  {
    phases$lntau2Obs_g  <- -1
    phases$lntau2Obs_pg <- -1
  }

  if( !hypoCtl$hierSel )
  {
    phases$epsSelAlpha_pg   <- phases$lnSelAlpha_g
    phases$epsSelBeta_pg    <- phases$lnSelBeta_g
    phases$lnSelAlpha_g     <-  -1
    phases$lnSelBeta_g      <-  -1
    
  }

  if( !hypoCtl$densityDepM )
  {
    phases$lnm1       <- -1
    phases$epslnm1_p  <- -1
  }

  if( hypoCtl$densityDepM & !hypoCtl$ddMprocError )
  {
    phases$omegaM_pt <- -1
  }


  if( !hypoCtl$corrRecDevs)
    phases$off_diag_R <- -1

  if( !hypoCtl$corrMortDevs | hypoCtl$densityDepM)
    phases$off_diag_M <- -1

  if( !dataCtl$combSpawnIdx )
  {
    phases$lnqComb_pg       <- -1
    phases$lntauObsComb_pg  <- -1
  }

  if( all(hypoCtl$initFcode_p == 0 ) )
    phases$lnFinit_p  <- -1

  if( all(hypoCtl$initRcode_p == 0 ) )
    phases$lnRinit_p  <- -1

  if( juveMsource == "mean" )
    phases$lnMjuve <- -1

  if( all(hypoCtl$avgRcode_p == 0))
    phases$lnRbar_p <- -1

  if( nP == 1 )
  {
    phases$epsM_p         <- -1
    phases$epslnm1_p      <- -1
    phases$epsSelAlpha_pg <- -1
    phases$epsSelBeta_pg  <- -1
    phases$lntau2Obs_pg   <- -1
    phases$epsSteep_p     <- -1

  }

  # # Generate map list - this will have to be updated
  # # so that it's sensitive to useFleetsIdx
  # # Map out selectivity deviations
  # if( hypoCtl$tvSelAlpha )
  #   epsSelAlphaMap <- factor( seq(from = 100, by = 1, length = nSelDevs ))
  # else epsSelAlphaMap <- factor(rep(NA,nSelDevs))
  # if( hypoCtl$tvSelBeta )
  #   epsSelBetaMap <- factor( seq(from = 100 + nSelDevs + 1, by = 1, length = nSelDevs ))
  # else epsSelBetaMap <- factor(rep(NA,nSelDevs))

  # then map out the average selectivity paramaters
  fleetSelEst <- hypoCtl$fleetSelEst[gearNames]
  
  aveSelAlphaMap                  <- 300 + fleetSelEst
  aveSelAlphaMap[fleetSelEst < 0] <- NA
  for( g in 1:(nG) )
  {
    if( g %in% sokFleets)
      next
    
    if( (all(A_apgt[,,g,] < 0) ) & (all(mA_agt[,g,] < 0) ) )
      aveSelAlphaMap[g] <- NA
  }
  names(aveSelAlphaMap)           <- gearNames
  aveSelBetaMap                   <- aveSelAlphaMap + nG + 1

  # stockSelDevMaps
  stockSelDevMapIdx <- seq(600,by=1,length.out = nG*nP)


  # Now make the map array for stock specific devs
  epsSelAlphaMap_pg <- array( 600, 
                              dim = c(nP,nG),
                              dimnames = list(stockNames,gearNames) )


  for( p in 1:nP )
    epsSelAlphaMap_pg[p,] <- epsSelAlphaMap_pg[p,] + (p-1) * (nG) + fleetSelEst[gearNames] 
  
  epsSelAlphaMap_pg[,fleetSelEst[gearNames] < 0] <- NA


  for( g in 1:(nG-1) )  # HACK to avoid SOK
  {
    for( p in 1:nP )
      if( all(A_apgt[,p,g,] < 0 ) )
        epsSelAlphaMap_pg[p,g] <- NA

    if( hypoCtl$hierSel)
    {
      if( sum(is.na(epsSelAlphaMap_pg[,g]) ) >= 2 )
        epsSelAlphaMap_pg[,g] <- NA 

      if( sum(is.na(epsSelAlphaMap_pg[,g]) ) == 1 )
      {
        whichEst <- which(!is.na(epsSelAlphaMap_pg[,g]))[1]
        epsSelAlphaMap_pg[-whichEst,g] <- NA
      }
    }
  }

  epsSelBetaMap_pg                    <- epsSelAlphaMap_pg + nP * nG + 1

  # Now we want to map different fleet obs error variances on
  # and off
  # First make dummy arrays
  obstau2map_pg <- array(NA, dim = c(nP,nG) )
  obstau2map_g  <- array(NA, dim = c(nG) )
  lowMapPar <- 700
  for( g in 1:nG )
  {
    # Now loop and fill
    if( calcIndex[g] == 1 )
    {
      for( p in 1:nP)
        if( any(I_pgt[p,g,] > 0))
        {
          obstau2map_pg[p,g] <- lowMapPar
          lowMapPar <- lowMapPar + 1
        }

      if( any(mI_gt[g,] > 0 ) )
      {
        obstau2map_g[g] <- lowMapPar
        lowMapPar <- lowMapPar + 1
      }
    }
  }

  omegaMmap_pt <- array( data=1:(nP*(nT)), dim=c(nP,nT) )
  recDevMap_pt <- array( data=1:(nP*(nT)), dim=c(nP,nT-1) )

  if( hypoCtl$mapRecDevs )
  {
    for(p in 1:nP )
      recDevMap_pt[p,] <- 1:(nT-1)

    minRecDev <- min(tFirstRecDev_p)
    
    if(minRecDev > 1)
      recDevMap_pt[,1:min(tFirstRecDev_p-1)] <- NA
  }

  # identical M devs
  if( hypoCtl$mapMort )
  {
    if(!hypoCtl$densityDepM)
      for( t in 2:(nT-1))
          omegaMmap_pt[,t] <- t

    if(hypoCtl$densityDepM)
      for( t in 1:(nT-1))
          omegaMmap_pt[,t] <- t    
  }

  # Loop over stocks and fix omegaM_pt to 0 if the stock
  # has a later init year
  for( p in 1:nP )
  {
    if( hypoCtl$densityDepM & tInitModelYear_p[p] >= 1 )
      omegaMmap_pt[p,1:(tInitModelYear_p[p])] <- NA

    if( !hypoCtl$densityDepM )
      omegaMmap_pt[p,1:(tInitModelYear_p[p]+1)] <- NA

    if( tFirstRecDev_p[p] > 1)
      recDevMap_pt[p,1:(tFirstRecDev_p[p]-1)] <- NA
    
    recDevMap_pt[p,tLastRecDev_p[p]:(nT-1)] <- NA
  }

  # Map out qF parameters
  mortType_g      <- dataCtl$mortType_g[gearNames]
  catSeriesType_g <- dataCtl$catSeriesType[gearNames]
  mapqF_pg <- array(1:(nP*nG), dim = c(nP,nG))
  mapqF_pg[,(mortType_g == 0 )] <- NA
  
  mapqFjumps_pgt <- array(1:(nP*nG*nT), dim = c(nP,nG,nT))
  mapqFjumps_pgt[,,1] <- NA
  mapqFjumps_pgt[,(mortType_g==0 ),] <- NA
  # Just in case
  for( g in 1:nG )
  {
    sumCat <- sum(mC_gt[g,]) + sum(C_pgt[,g,])

    if(sumCat == 0)
      mapqF_pg[,g]        <- NA
    
    rwOn <- FALSE
    for(t in 1:nT)
    {
      sumCat <- sum(mC_gt[g,t]) + sum(C_pgt[,g,t])
      if(sumCat == 0 | !rwOn)
      {
        mapqFjumps_pgt[,g,t] <- NA
      }
      if(sumCat > 0)
        rwOn <- TRUE
    }
  }

  # May want to start Mdevs later
  tInitMdev_p <- initMdev_p - fYear
  for( p in 1:nP )
    omegaMmap_pt[p,1:(tInitMdev_p[p])] <- NA

  # Make fDevs_ap map

  fDevsMap_ap <- array( data = seq( from = 1, 
                                    by = 1, 
                                    length.out = (nFishedInitAges)*nP),
                        dim = c(nFishedInitAges,nP) )

  
  fDevsMap_ap[,initFished == 0] <- NA
  fDevsMap_ap[,initDevCode_p == 0] <- NA

  if( hypoCtl$mapInitDevs )
  {
    for( p in 1:nP )
      fDevsMap_ap[,p] <- 1:(nFishedInitAges-1)
  }


  maplogitphi1_g  <- seq(1e4, by = 1, length.out = nG)
  maplogitphi1_g[is.na(aveSelAlphaMap[gearNames])] <- NA
  if(nG > 5)
    maplogitphi1_g[6:nG] <- NA

  if( hypoCtl$ageCompLik == 0 )
  {
    phases$logitphi1_g  <- -1
    phases$logitpsi_g   <- -1
  }


  maplogitpsi_g   <- maplogitphi1_g + nG

  if( !hypoCtl$hierSel )
  {
    aveSelAlphaMap  <- rep(NA,nG)
    aveSelBetaMap   <- rep(NA,nG)
  }

  # Maps for initial conditions
  lnFinitMap_p <- seq(from = 1e3, by = 1, length.out = nP )
  lnFinitMap_p[initFcode_p == 0] <- NA
  # Maps for initial conditions
  lnRinitMap_p <- seq(from = 1e3+nP+1, by = 1, length.out = nP )
  lnRinitMap_p[initRcode_p == 0] <- NA  

  lnRbarMap_p <- seq(from = 1e3+2*nP+1, by = 1, length.out = nP )

  lnRbarMap_p[avgRcode_p == 0] <- NA    

  # Map for lnSDprobPosIdx_g
  probPosIdxMap_pg <- array(seq( from = 3*1e3, by = 1, length.out = nG*nP ),
                              dim = c(nP,nG))

  mlnqMap_g <- probPosIdxMap_pg[1,] + 3*(nG*nP)
  mlnqMap_g[hypoCtl$qPrior_g[gearNames] != 2] <- NA
  
  probPosIdxMap_pg[,dataCtl$deltaIdx[gearNames] <= 0] <- NA

  if( hypoCtl$mapDeltaModel & nP > 1 )
    for( p in 2:nP )
      probPosIdxMap_pg[p,] <- probPosIdxMap_pg[1,]

  probPosIdxMap_pg[idxWithZeroes_pg == 0] <- NA

  # Map for combo q
  mapComboQ_pg    <- array( seq(from = 4*1e3, by = 1, length.out = nP*nG ), dim =c(nP,nG) )
  mapComboTau_pg  <- mapComboQ_pg + nP * nG + 1
  mapComboQ_pg[,dataCtl$fleetType[gearNames] != 0] <- NA
  mapComboQ_pg[,dataCtl$idxType[gearNames] != 0] <- NA
  mapComboTau_pg[,dataCtl$fleetType[gearNames] != 0] <- NA

  map <- list(  lnSelAlpha_g        = factor(aveSelAlphaMap),
                lnSelBeta_g         = factor(aveSelBetaMap),
                epsSelAlpha_pg      = factor(epsSelAlphaMap_pg),
                epsSelBeta_pg       = factor(epsSelBetaMap_pg),
                lntau2Obs_g         = factor(obstau2map_g),
                lntau2Obs_pg        = factor(obstau2map_pg),
                omegaM_pt           = factor(omegaMmap_pt),
                fDevs_ap            = factor(fDevsMap_ap),
                logitphi1_g         = factor(maplogitphi1_g),
                logitpsi_g          = factor(maplogitpsi_g),
                lnFinit_p           = factor(lnFinitMap_p),
                lnRinit_p           = factor(lnRinitMap_p),
                lnSDProbPosIdx_pg   = factor(probPosIdxMap_pg),
                meanProbPosIdx_pg   = factor(probPosIdxMap_pg + nG*nP),
                lnRbar_p            = factor(lnRbarMap_p),
                mlnq_g              = factor(mlnqMap_g),
                lnqComb_pg          = factor(mapComboQ_pg),
                lntauObsComb_pg     = factor(mapComboTau_pg),
                recDevs_pt          = factor(recDevMap_pt),
                lnqFinit_pg         = factor(mapqF_pg),
                deltaqF_pgt         = factor(mapqFjumps_pgt) )  



  # Run TMBphase
  phzList <- TMBphase(  data = data, 
                        parameters = pars, 
                        random = hypoCtl$RE,
                        phases = phases, 
                        base_map = map,
                        maxPhase = ctrlObj$maxPhase,
                        model_name = "SISCA",
                        optimizer = "nlminb",
                        silent = ctrlObj$quiet,
                        calcSD = ctrlObj$calcSD,
                        maxEval = ctrlObj$maxFunEval,
                        maxIter = ctrlObj$maxIterations, 
                        mcIter  = ctrlObj$mcmcIterations,
                        nChain  = ctrlObj$nChain,
                        mcDisp  = ctrlObj$mcDisp,
                        mcSeed  = ctrlObj$mcSeed,
                        checkPhase = ctrlObj$checkEachPhase,
                        adapt_delta = ctrlObj$adapt_delta,
                        max_treedepth = ctrlObj$max_treedepth,
                        makePostStates = ctrlObj$makePostStates,
                        turnOffRinit = hypoCtl$turnOffRinit,
                        turnOffRbar = hypoCtl$turnOffRbar,
                        forceHMC = ctrlObj$forceHMC ) 
  

  repOpt <- renameReportArrays(phzList$repOpt,data)

  report <- list( repInit   = phzList$repInit,
                  repOpt    = repOpt,
                  sdrepOpt  = phzList$sdrep,
                  posts     = phzList$posts,
                  fYear     = fYear, 
                  lYear     = lYear,
                  gearLabs  = gearNames,
                  stock     = stockNames,
                  data      = data,
                  pars      = pars,
                  initmap   = map,
                  phases    = phases,
                  finalmap  = phzList$map,
                  fitReport = phzList$fitReport,
                  grad      = phzList$grad,
                  stanfit   = phzList$stanfit )

  return(report)
} # END .runSISCA()



# Make a parameter table for fits under the MLEs,
# need to add uncertainties to this
makeParFitTableMLE <- function( repList )
{
  ctlList <- repList$ctlList
  repOpt  <- repList$repOpt
  sdrep   <- repList$sdrepOpt

  if( !is.null(ctlList$ctrl$iscamRep))
  {
    iscamRepPath <- paste("Data/",ctlList$ctrl$iscamRep,sep="")
    ISCAMrep <-read.rep(iscamRepPath)
  }

  

  colNames <- c(  "dataScenario",
                  "modelHyp",
                  "Stock",
                  "$B_0$",   
                  "$R_0$",   
                  "$R_{init}$",
                  "$F_{init}$",
                  "$M_0$",   
                  "$\\overline{M}$",
                  "$q_s$",   
                  "$q_d$",
                  "pdHess" ) 

  pdHess <- sdrep$pdHess

  nP <- repOpt$nP

  parTable <- matrix(NA, nrow = nP, ncol = length(colNames) )
  colnames(parTable) <- colNames
  parTable <- as.data.frame(parTable)

  
  for( p in 1:nP )
  {

    rowIdx <- p
    parTable[rowIdx,"dataScenario"]     <- ctlList$ctrl$dataScenarioName
    parTable[rowIdx,"modelHyp"]         <- ctlList$ctrl$modelHypName
    parTable[rowIdx,"Stock"]            <- repList$stock[p]
    parTable[rowIdx,"$B_0$"]            <- round(repOpt$B0_p[p],3)
    parTable[rowIdx,"$R_0$"]            <- round(repOpt$R0_p[p],3)
    parTable[rowIdx,"$R_{init}$"]       <- round(repOpt$Rinit_p[p],3)
    parTable[rowIdx,"$F_{init}$"]       <- round(repOpt$Finit_p[p],3)
    parTable[rowIdx,"$M_0$"]            <- round(repOpt$M0_p[p],3)
    parTable[rowIdx,"$\\overline{M}$"]  <- round(repOpt$meanM_p[p],3)
    parTable[rowIdx,"$q_s$"]            <- round(repOpt$qhat_pg[p,4],3)
    parTable[rowIdx,"$q_d$"]            <- round(repOpt$qhat_pg[p,5],3)
    parTable[rowIdx,"pdHess"]           <- pdHess
  }

  outList <- list(  parEstTable = parTable)

  if(!is.null(ctlList$ctrl$iscamRep))
  {
    similarityTable <- calcSimilarityISCAMvSISCAH(repList, ISCAMrep)  

    outList$similarityTable <- similarityTable
  }
  outList
}


# Custom TMBphase() function for running hierSCAL in phases. 
# Modified from the version in Kasper Kristensen's TMB_contrib_R github repository 
# https://github.com/kaskr/TMB_contrib_R
# Author:Gavin Fay email: gfay42@gmail.com
# 
# Main modification adds a base map list for array based parameters that
# have (a) identified parameters or (b) only some elements
# activated due to missing data (e.g. selectivity pars for
# fleets in specific areas). Doing it this way reduces number
# of loops in the model, speeding up fitting time.
TMBphase <- function( data, 
                      parameters, 
                      random, 
                      phases, 
                      base_map = list(),
                      maxPhase = NULL,
                      model_name = "assessCA",
                      optimizer = "nlminb",
                      silent = FALSE,
                      calcSD = FALSE,
                      maxEval = 1000,
                      maxIter = 1000,
                      mcIter  = 0,
                      nChain  = 4,
                      mcDisp  = 1,
                      mcSeed  = 123,
                      adapt_delta = 0.8,
                      checkPhase = FALSE,
                      max_treedepth = 12,
                      makePostStates = FALSE,
                      turnOffRinit = FALSE,
                      turnOffRbar = FALSE,
                      forceHMC = FALSE ) 
{

  nP <- dim(data$I_pgt)[1]

  # function to fill list component with a factor
  # of NAs
  fill_vals <- function(x,vals){ factor( rep( vals, length(x) ) ) }

  # compile the model
  DLL_use <- model_name  
  
  # set maximum phase
  if(!is.null(maxPhase))
    maxPhase <- min( maxPhase, max(unlist(phases) + 1 ) )
  else
    maxPhase <- max(unlist(phases)) + 1

  # Make a data.frame that will hold the phase fit info
  fitRepColNames <- c("phase","objFun","maxGrad","nPar","convCode","convMsg", "time", "pdHess")
  fitReport <- matrix(NA, nrow = maxPhase + 1, ncol = length(fitRepColNames))
  colnames(fitReport) <- fitRepColNames
  
  fitReport <- as.data.frame(fitReport)

  fitReport$phase <- c(1:maxPhase,"")

  phaseReports <- vector(mode = "list", length = maxPhase)

  # generate a list of outputs to return
  # to runHierSCAL, initialise 
  # a success flag at TRUE
  outList <- list( success = TRUE )


  for( phase_cur in 1:maxPhase ) 
  {
    # Start timing
    tic("phaseTimer")

    if(phase_cur == maxPhase )
    {
      # Basically want to turn off Rbar estimation
      if(turnOffRbar )
      {
        data$avgRcode_p   <- rep(0,nP)
        phases$lnRbar_p   <- -1
      }
      
      if(turnOffRinit)
      {
        data$initRcode_p  <- rep(0,nP)
        phases$lnRinit_p  <- -1
      }
      # phases$lnm1       <- -1
      # phases$lnM        <- phase_cur


      # I should back-calculate the recdevs?

    }

    # work out the map for this phase
    # if the phase for a parameter is greater than the current phase 
    # or a negative value, then map will contain a factor filled with NAs
    map_use <- base_map
    j <- length(map_use)
    for( i in 1:length(parameters) ) 
    {
      parName <- names(parameters)[i]


      if( parName %in% names(phases) )
      {
        if( (phases[[parName]] > phase_cur) | phases[[parName]] < 0 ) 
        { 
          
          # Check if parName is included in the base_map
          if(parName %in% names(map_use))
            map_use[[parName]] <- fill_vals(parameters[[i]],NA)
          else
          {
            j <- j + 1
            map_use[[j]] <- fill_vals(parameters[[i]],NA)
            names(map_use)[j] <- parName
          }

        }
      } else {
        if( parName %in% names(map_use) )
          map_use[[parName]] <- fill_vals(parameters[[i]],NA)
        else {
          j <- j + 1
          map_use[[j]] <- fill_vals(parameters[[i]],NA)
          names(map_use)[j] <- parName
        }
      }

    }

    #remove the random effects if they are not estimated
    random_use <- random[!random %in% names(map_use)]
  
    # initialize the parameters at values in previous phase
    params_use <- parameters
    if( phase_cur > 1 ) 
      params_use <- obj$env$parList( opt$par )

    if( phase_cur == maxPhase)
    {
      repObj    <- obj$report()
      if(turnOffRbar)
      {
        SRdevs_pt <- repObj$SRdevs_pt[,-1,drop = FALSE]
        newRecDevs_pt <- SRdevs_pt
        # Issue here, as the back-transformation here may
        # differ from the logit transformation inside the cpp. 
        # We should generalise this with ctl file pars
        newRecDevs_pt <- -log( -1 +  10/(SRdevs_pt + 5))
        # last few rec devs are not estimated
        newRecDevs_pt[params_use$recDevs_pt == 0] <- 0
        params_use$recDevs_pt <- newRecDevs_pt
      }

      # Need to re-estimate init rec devs
      if( turnOffRinit )
      {
        N_ap1 <- repObj$N_apt[,,1,drop = FALSE]

        fDevsNew_ap <- array(0,dim = dim(N_ap1)[1:2])
        for( p in 1:nP )
          fDevsNew_ap[,p] <- log(N_ap1[,p,1]/(repObj$R0_p[p] * repObj$initSurv_ap[,p]))/repObj$sigmaR

        params_use$fDevs_ap <- fDevsNew_ap
      }
    }

    map_use <- map_use[names(map_use) %in% names(params_use)]

    # if(phase_cur == 7)
    #   browser()

    # Fit the model
    obj <- TMB::MakeADFun(  data = data,
                            parameters = params_use,
                            random= NULL,
                            DLL= DLL_use,
                            map= map_use,
                            silent = silent )  

    TMB::newtonOption(obj,trace = 10, tol10 = 0.01 )

    if( phase_cur == 1 )
    {
      repInit <- obj$report()

      outList$repInit <- repInit


      checkInit   <- lapply( X = outList$repInit, FUN = is.nan )
      checkFinite <- sapply( X = outList$repInit, FUN = is.finite )

      if(any(unlist(checkInit)) | any(!unlist(checkFinite)))
        browser(beep(expr=cat("NaN items in repInit\n")))
      
    }
  
    
    # Create a control list for the assessment model
    tmbCtrl <- list(  eval.max = maxEval, 
                      iter.max = maxIter  )
    if(!silent)
      cat("\nStarting optimisation for phase ", phase_cur, "\n\n")

    # Try the optimisation
    # options(warn = 1)
    opt <- try( nlminb( start     = obj$par,
                        objective = obj$fn,
                        gradient  = obj$gr,
                        control   = tmbCtrl ) )

    if(!silent)
    {
      message("B0 = ", obj$report()$B0_p,"\n")
      message("Rbar = ", obj$report()$Rbar_p,"\n")
      message("Rinit = ", obj$report()$Rinit_p,"\n")
    }

    # break if there is an issue
    if( class(opt) == "try-error" )
    {

      cat("\nOptimisation halted due to error\n")

      outList$success                   <- FALSE
      phaseReports[[phase_cur]]$opt     <- opt
      phaseReports[[phase_cur]]$success <- FALSE
      outList$maxPhaseComplete          <- phase_cur - 1
      break
    }

    # Update fitReport
    if(class(opt) != "try-error")
    {
      fitReport[phase_cur,]$objFun      <- obj$fn()
      fitReport[phase_cur,]$maxGrad     <- max(abs(obj$gr()))
      fitReport[phase_cur,]$nPar        <- length(opt$par)
      fitReport[phase_cur,]$convCode    <- opt$convergence
      fitReport[phase_cur,]$convMsg     <- opt$message
      fitReport[phase_cur,]$pdHess      <- obj$hessian
    }
    time <- toc(quiet = TRUE)
    fitReport[phase_cur,]$time           <- time$toc - time$tic

    # Save reports and optimisation
    # output
    phaseReports[[phase_cur]]$report  <- obj$report()
    phaseReports[[phase_cur]]$opt     <- opt
    phaseReports[[phase_cur]]$success <- TRUE
    phaseReports[[phase_cur]]$map     <- map_use
    outList$maxPhaseComplete          <- phase_cur

    if(!silent)
    {
      cat(  "\nPhase ", phase_cur, " completed with code ",
            opt$convergence, " and following message:\n", sep = "" )
      cat("\n", opt$message, "\n\n", sep = "" )
    }
    
    # close phase loop
  } # end for( phase_cur in 1:maxPhase ) 

  # Extra Newton steps for to improve MGCs
  if( phase_cur == maxPhase & class(opt) != "try-error" )
  {
    message( "Max phase reached and optimised, trying additional Newton steps to improve gradients.\n")

    params_use <- obj$env$parList(opt$par)

    # # Can we re-weight here?
    # rpt <- obj$report()
    # SRdevs_pt <- rpt$SRdevs_pt
    # effsigR <- sd(SRdevs_pt[SRdevs_pt != 0])

    # Mdevs_pt <- rpt$omegaM_pt
    # effsigM <- sd(Mdevs_pt[Mdevs_pt != 0])

    # params_use$lnsigmaR <- params_use$lnsigmaR + log( round(effsigR, 2) )
    # params_use$lnsigmaM <- params_use$lnsigmaM + log( round(effsigM, 2) )

    # Recreate the object
    obj <- TMB::MakeADFun(  data = data,
                            parameters = params_use,
                            random= NULL,
                            DLL= DLL_use,
                            map= map_use,
                            silent = silent )  

    # Try the optimisation
    opt <- try( nlminb (  start     = obj$par,
                          objective = obj$fn,
                          gradient  = obj$gr,
                          control   = tmbCtrl ) )

    # Update fitReports
    if(class(opt) != "try-error")
    {
      fitReport[maxPhase,]$objFun      <- obj$fn()
      fitReport[maxPhase,]$maxGrad     <- max(abs(obj$gr()))
      fitReport[maxPhase,]$nPar        <- length(opt$par)
      fitReport[maxPhase,]$convCode    <- opt$convergence
      fitReport[maxPhase,]$convMsg     <- opt$message
      fitReport[maxPhase,]$pdHess      <- obj$hessian
    }

    if(!silent)
    {
      cat(  "\nAdditional Newton steps completed with code ",
            opt$convergence, " and following message:\n", sep = "" )
      cat("\n", opt$message, "\n\n", sep = "" )
    }
  }

  # REs
  if( !is.null(random) &  class(opt) != "try-error" )
  { 
    fitReport$phase[maxPhase+1] <- "re"
    tic("re")
    params_use  <- obj$env$parList( opt$par )

    # MakeADFun
    obj <- TMB::MakeADFun(  data = data,
                            parameters = params_use,
                            random = random_use,
                            DLL= DLL_use,
                            map= map_use,
                            silent = silent )  
    TMB::newtonOption(obj,trace = 10, tol10 = 0.0001 )

    # Try the optimisation
    opt <- try( nlminb (  start     = obj$par,
                          objective = obj$fn,
                          gradient  = obj$gr,
                          control   = tmbCtrl ) )

    # Update fitReports
    if(class(opt) != "try-error")
    {
      fitReport[maxPhase + 1,]$objFun      <- obj$fn()
      fitReport[maxPhase + 1,]$maxGrad     <- max(abs(obj$gr()))
      fitReport[maxPhase + 1,]$nPar        <- length(opt$par)
      fitReport[maxPhase + 1,]$convCode    <- opt$convergence
      fitReport[maxPhase + 1,]$convMsg     <- opt$message
      fitReport[maxPhase + 1,]$pdHess      <- obj$hessian
    }

    time <- toc(quiet = TRUE)
    fitReport[maxPhase + 1,]$time          <- time$toc - time$tic

  }

  stanfit <- NULL
  posts <- NULL
  if(outList$success & calcSD )
  {
    tic("posteriors")
    outList$sdrep   <- TMB::sdreport(obj)

    grad <- summary(outList$sdrep)[1:length(opt$par),]
    grad <- cbind( grad, as.vector(obj$gr()), grad[ ,2]/abs(grad[ ,1]) )
    colnames(grad) <- c("est","se","gr","cv")
    grad <- as.data.frame(grad)
    outList$grad <- grad

    MLErpt <- obj$report()
    MLErpt$R2_pt <- array( data=NA, dim=dim(MLErpt$N_apt)[2:3] )
    MLErpt$R2_pt[1, ] <- MLErpt$N_apt[2,1, ]


    if(!mcIter )
    {
      fitReport$phase[maxPhase + 1]         <- "SDreport"
      fitReport[maxPhase + 1,]$pdHess       <- obj$hessian
      time <- toc(quiet = TRUE)
      fitReport$time[maxPhase + 1]          <- time$toc - time$tic
    }

    # Run MCMC
    if( mcIter )
    {

      # Initial conditions for each chain
      npar  <- length(opt$par)
      mcinit <- list()

      if(any(!is.finite(grad$se)) & !forceHMC)
      {
        message("\n\n Non-finite SEs in SD report! Improve convergence and try again.\n\n")

        readline(prompt="\n Press [enter] to skip HMC and create fit report for non-convergent MLEs \n")
      }

      if( all(is.finite(grad$se)) | forceHMC )
      {

        message("\nStarting HMC NUTS algorithm in tmbstan. This may take some time.\n")

        for( i in 1:nChain )
        {
          SDs <- mcDisp*grad$se
          SDs[!is.finite(SDs)] <- 0.0001*abs(grad$est)
          set.seed(mcSeed)
          mcinit[[i]] <- rnorm(n=nrow(grad),mean=grad$est,sd=SDs)

          nTries <- 0
          while(!is.finite(obj$fn(mcinit[[i]])) & nTries <= 40 )
          {
            browser()
            nTries <- nTries+1
            mcinit[[i]] <- rnorm(n=nrow(grad),mean=grad$est,sd=SDs)          
            
          }
        }



        # Draw MCMC samples
        rpt <- obj$report()
        mc_cores <- parallel::detectCores()-1
        stanfit <- tmbstan( obj = obj,
                            chains = nChain,
                            iter = mcIter,
                            init = mcinit,
                            control = list( adapt_delta=adapt_delta,
                                            max_treedepth=max_treedepth ),
                            cores = mc_cores)
        
        posts <- list()

        # Add timing to fit report
        time <- toc(quiet = TRUE)
        fitReport$phase[maxPhase + 1]         <- "HMC"
        fitReport[maxPhase + 1,]$time         <- time$toc - time$tic

        if( makePostStates )
          posts <- makePosts(   stanfit = stanfit,
                                repObj  = MLErpt,
                                data    = data,
                                pars    = parameters,
                                map     = map_use,
                                refPts  = ifelse(data$densityDepM == 1, FALSE, TRUE) )
      }

    }


  }


  outList$phaseReports  <- phaseReports
  outList$repOpt        <- MLErpt
  outList$objfun        <- obj$fn() 
  outList$maxGrad       <- max(abs(obj$gr()))
  outList$fitReport     <- fitReport
  outList$totTime       <- sum(fitReport$time)
  outList$stanfit       <- stanfit
  outList$posts         <- posts
  outList$map           <- map_use
  
  return( outList )  

} # END TMBphase()



# .readParFile   (reads an ASCII file with 1 comment line, header, data frame)
# Purpose:      Reads an ASCII file: 1 comment, 1 header, space-delimited
#               data frame usually containing columns "parameter" and "value".
# Parameters:   parFile is a character string indicating the input file.
# Returns:      result, a data frame.
# Source:       A.R. Kronlund
.readParFile <- function( parFile="inputFile.par" )
{
  # Read the file and store as a dataframe.
  result <- read.table( file=parFile, as.is=TRUE, header=TRUE, skip=1,
                        quote="",sep=" " )
  result
}


.addQuotes <- function( str ) 
{
  # Adds double-quotes to a string.
  return( paste("\"", str, "\"", sep = "") )
}

.convertSlashes <- function( expr, os = .Platform$OS.type ) 
{
  if ( os == "windows" ) 
    expr = gsub("/", "\\\\", expr)
  else expr = gsub("\\\\", "/", expr)
  return(expr)
}

.createList <- function( obj )
{
  # Input  a data frame with columns "parameter" and "value".
  # Output a list with elements named as parameter and values in "value".

  result <- list()

  # Shut off whining, then coerce to numeric to let NA indicate non-numerics.
  options( warn=-1 )
  numericVal <- as.numeric( obj[,"value"] )
  options( warn=0 )

  for ( i in 1:nrow(obj) )
  {
    # Value is numeric, build the parse string.
    if ( !is.na(numericVal[i]) )
      listText <- paste( "result$",obj[i,"parameter"],"=",
                    obj[i,"value"],sep="" )
    # Value is character, build the parse string.
    else
      listText <- paste( "result$",obj[i,"parameter"],"=",
                  obj[i,"value"], sep="" )

    # ARK: At one point I had this code, presumably to handle a different
    #      input format convention, perhaps assuming "value" was all character.
    #                   sQuote(obj[i,"value"]),sep="" )
    
    # Evaluate the parse string.
    eval( parse( text=listText ) )
  }
  result
}

# .calcRampedHCR()
# Calculates ramped HCR based
# on control points, biomass estimate, and reference
# removal rate
.calcRampedHCR <- function( B, LCP, UCP, Fref, lowFmult = 0 )
{
  if( B < LCP )
    F <- lowFmult * Fref
  if( LCP <= B & B < UCP )
    F <- (lowFmult + (1 - lowFmult) * (B - LCP) / (UCP - LCP)) * Fref
  if( UCP <= B )
    F <- Fref

  F
} # END calcRampedHCR

# .calcLegalCatch()
# Legal sized TAC calculation
.calcLegalCatch <- function( B, ...)
{
  targF <- .calcRampedHCR(B = B, ...)

  legalCatch <- targF * B

  return(legalCatch) 
} # END .calcLegalCatch
