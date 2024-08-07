# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
# SISCAtools.R
#
# Functions for loading, saving, and modifying SISCA fit objects
# and associated reports.
#
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

makeEstPhaseTable <- function( repObj = reports )
{
  phases <- repObj$ctlList$phases


  phaseTable <- data.frame(phase = unlist(phases))
  rownames(phaseTable) <- names(phases)

  phaseTable
}


# loadFit()
# Loads the nominated fit reports object into memory, 
# so that plot functions can be called
# inputs:   fit=ordinal indicator of sim in project folder
# ouputs:   NULL
# usage:    Prior to plotting simulation outputs
.loadFit <- function( fit = 1, groupFolder = "", 
                      retObj = FALSE, quiet = FALSE,
                      baseDir = "Outputs" )
{
  groupFolder <- here::here(baseDir,groupFolder)
  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path=groupFolder,full.names = FALSE,
                        recursive=FALSE)
  # Restrict to sim_ folders, pick the nominated simulation
  fitList <- dirList[grep(pattern="fit",x=dirList)]

  # Load fit object
  folder <- fit
  if( is.character(fit))
      if(!grepl(pattern="fit_",x=fit) )
        folder <- paste("fit_",fit,sep = "")

  if(!is.character(fit))
    folder <- fitList[fit]

  # Load the nominated blob
  reportsFileName <- paste(folder,".rds",sep="")
  reportsPath <- file.path(groupFolder,folder,reportsFileName)
  reports <- readRDS( file = reportsPath )

  if(folder != basename(reports$folder))
    reports$folder <- folder

  if(!quiet)
    cat("MSG (loadFit) Reports in ", folder, " loaded from ", groupFolder, "\n", sep="" )

  if(!retObj)
      assign( "reports", reports, pos=1 )

  if(retObj )
    return(reports)
} # END .loadFit

# makePosts()
# Refactored process to make posteriors from 
# tmbstan samples. Reduces copy/paste errors when resampling
# is required
makePosts <- function(  stanfit = reports$stanfit,
                        repObj  = reports$repOpt,
                        data    = reports$data,
                        pars    = reports$pars,
                        map     = reports$map,
                        refPts  = FALSE )
{
  # Pull samples
  samps <- as.data.frame(stanfit)
  samps$lp__ <- NULL

  # Pull model dimensions
  nA <- repObj$nA
  nP <- repObj$nP
  nG <- repObj$nG
  nT <- repObj$nT
  nV <- dim(repObj$catTable_pvk)[2]

  mcIter <- dim(samps)[1]

  mon <- monitor(stanfit)
  samp <- as.array(stanfit)

  npTable <- nuts_params(stanfit)

  divTable <- npTable |>
              filter( Parameter == "divergent__",
                      Value == 1 )




  # Need to create a TMB object
  # Fit the model
  obj <- TMB::MakeADFun(  data = data,
                          parameters = pars,
                          random= NULL,
                          DLL= "SISCA",
                          map= map,
                          silent = FALSE )

  if(refPts)
  {
    tmpRefPts <- calcRefPts(repObj)$refPts

    Fgrid     <- tmpRefPts$refCurves$F
    nFs       <- length(Fgrid)

  }

  colNames <- c("Chain","Iteration",dimnames(samp)[[3]])
  rowNames <- 1:nrow(divTable)
  divergencesTable <- array(NA, dim = c(length(rowNames),length(colNames)),
                                  dimnames = list(rowNames,colNames))
  divGradTable      <- as.data.frame(divergencesTable[,-length(colNames)])
  divergencesTable  <- as.data.frame(divergencesTable)
  if(nrow(divTable) > 0 )
  {
    message("Calculating gradients at each divergent transition...")
    
    for( k in 1:nrow(divTable))
    {
      chain <- divTable[k,1]
      iter  <- divTable[k,2]

      divergencesTable[k,1]       <- divTable[k,1]
      divergencesTable[k,2]       <- divTable[k,2]
      divergencesTable[k,-(1:2)]  <- samp[iter,chain,]

      divGradTable[k,1]       <- divTable[k,1]
      divGradTable[k,2]       <- divTable[k,2]
      divGradTable[k,-(1:2)]  <- obj$gr(samp[iter,chain,-dim(samp)[3]])

    }
  }



  # Posteriors of quantities of interest
  message("Calculating model posteriors from samples...")
  posts <- list(  SB_ipt             = array( data=NA,  dim=c(mcIter,nP,nT+1) ),
                  zComb_ipt          = array( data=NA,  dim=c(mcIter,nP,nT) ),
                  R_ipt              = array( data=NA,  dim=c(mcIter,nP,nT+1) ),
                  R2_ipt             = array( data=NA,  dim=c(mcIter,nP,nT+1) ),
                  bhR_ipt            = array( data=NA,  dim=c(mcIter,nP,nT+1) ),
                  omegaR_ipt         = array( data=NA,  dim=c(mcIter,nP,nT+1) ),
                  sigmaR_i           = array( data=NA,  dim=c(mcIter) ),
                  omegaM_ipt         = array( data=NA,  dim=c(mcIter,nP,nT) ),
                  SRdevs_ipt         = array( data=NA,  dim=c(mcIter,nP,nT) ),
                  M_iapt             = array( data=NA, dim=c(mcIter,nA,nP,nT+1) ),
                  q_ipg              = array( NA, dim = c(mcIter,nP,nG) ),
                  selAlpha_ipg       = array( data=NA, dim=c(mcIter,nP,nG) ),
                  selBeta_ipg        = array( data=NA, dim=c(mcIter,nP,nG) ),
                  B0_ip              = array( NA, dim = c(mcIter,nP)),
                  h_ip               = array( NA, dim = c(mcIter,nP)),
                  totB0_ip           = array( NA, dim = c(mcIter,nP)),
                  M_ip               = array( NA, dim = c(mcIter,nP)),
                  M0_ip              = array( NA, dim = c(mcIter,nP)),
                  m1_ip              = array( NA, dim = c(mcIter,nP)),
                  U_ipgt             = array( NA, dim = c(mcIter,nP,nG,nT)),
                  initNmult_iap      = array( NA, dim = c(mcIter,nA,nP)),
                  Rbar_ip            = array( NA, dim = c(mcIter,nP)),
                  Rinit_ip           = array( NA, dim = c(mcIter,nP)),
                  R0_ip              = array( NA, dim = c(mcIter,nP)),
                  tauComb_ipg        = array( NA, dim = c(mcIter,nP,nG)), 
                  probPosIdx_ipgt    = array( NA, dim = c(mcIter,nP,nG,nT)),
                  meanProbPosIdx_ipg = array( NA, dim = c(mcIter,nP,nG)),
                  SDProbPosIdx_ipg   = array( NA, dim = c(mcIter,nP,nG)),
                  corrM_ipp          = array( NA, dim = c(mcIter,nP,nP)),
                  corrR_ipp          = array( NA, dim = c(mcIter,nP,nP)),
                  catTable_ipvk      = array( NA, dim = c(mcIter,nP,nV,4)),
                  totLike            = rep( NA, mcIter ) )

  if(refPts)
  {
    posts <- c(posts,
                list( Bmsy_ip       = array(NA, dim = c(mcIter,nP)),
                      Umsy_ip       = array(NA, dim = c(mcIter,nP)),
                      MSY_ip        = array(NA, dim = c(mcIter,nP)),
                      SBeq_ipf      = array(NA, dim = c(mcIter,nP,nFs)),
                      Yeq_ipf       = array(NA, dim = c(mcIter,nP,nFs)),
                      Req_ipf       = array(NA, dim = c(mcIter,nP,nFs)),
                      Ueq_ipf       = array(NA, dim = c(mcIter,nP,nFs)) ))

    posts$F <- F
  }

  # Progress bar
  pb <- txtProgressBar( min=0, max=mcIter, style=3 )
  # Loop over each set of parameter estimates
  
  for( i in 1:mcIter )
  {
    r <- obj$report(samps[i, ])

    posts$SB_ipt[i, , ]           <- r$SB_pt
    posts$zComb_ipt[i, , ]        <- r$zComb_pt
    posts$R_ipt[i, , ]            <- r$R_pt
    posts$R2_ipt[i, , ]           <- r$N_apt[2, , ]
    posts$bhR_ipt[i, , ]          <- r$bhR_pt
    posts$omegaR_ipt[i, , ]       <- r$omegaR_pt
    posts$sigmaR_i[i]             <- r$sigmaR
    posts$omegaM_ipt[i, , ]       <- r$omegaM_pt
    posts$SRdevs_ipt[i, , ]       <- r$SRdevs_pt
    posts$M_iapt[i, , , ]         <- r$M_apt
    posts$q_ipg[i, , ]            <- r$qhat_pg
    posts$selAlpha_ipg[i,,]       <- r$SelAlpha_pg
    posts$selBeta_ipg[i,,]        <- r$SelBeta_pg
    posts$B0_ip[i,]               <- r$B0_p
    posts$h_ip[i,]                <- r$rSteepness_p
    posts$totB0_ip[i,]            <- r$totB0_p
    posts$M_ip[i,]                <- r$M_p
    posts$M0_ip[i,]               <- r$M0_p
    posts$m1_ip[i,]               <- r$m1_p
    posts$U_ipgt[i,,,]            <- r$U_pgt 
    posts$initNmult_iap[i,,]      <- r$initN_mult_ap
    posts$Rbar_ip[i,]             <- r$Rbar_p
    posts$Rinit_ip[i,]            <- r$Rinit_p
    posts$R0_ip[i,]               <- r$R0_p
    posts$probPosIdx_ipgt[i,,,]   <- r$probPosIdx_pgt
    posts$meanProbPosIdx_ipg[i,,] <- r$meanProbPosIdx_pg
    posts$SDProbPosIdx_ipg[i,,]   <- r$SDProbPosIdx_pg
    posts$corrM_ipp[i,,]          <- r$corrM_pp
    posts$corrR_ipp[i,,]          <- r$corrR_pp
    posts$catTable_ipvk[i,,,]     <- r$catTable_pvk
    posts$tauComb_ipg[i,,]        <- r$tauComb_pg
    posts$totLike[i]              <- r$totLike


    if(refPts)
    {
      # Might need to beef up r with some info, but try this for now
      rp <- calcRefPts(r)$refPts

      posts$Umsy_ip[i,]   <- rp$FmsyRefPts$Umsy_p
      posts$Bmsy_ip[i,]   <- rp$FmsyRefPts$SBeqFmsy_p
      posts$MSY_ip[i,]    <- rp$FmsyRefPts$YeqFmsy_p

      posts$SBeq_ipf[i,,] <- rp$refCurves$SBeq_pf
      posts$Req_ipf[i,,]  <- rp$refCurves$Req_pf
      posts$Yeq_ipf[i,,]  <- rp$refCurves$Yeq_pf
      posts$Ueq_ipf[i,,]  <- rp$refCurves$Ueq_pf
    }

    setTxtProgressBar(pb, i)



  }
  close(pb)

  posts$divTable      <- divergencesTable
  posts$divGradTable  <- divGradTable

  posts
} # END makePosts()


# redoPosts()
# Will rerun the postierior generating procedure,
# place posteriors into the reports object, and
# save to the fit directory
redoPosts <- function(  fitID = 1, 
                        groupFolder = "",
                        refPts = FALSE,
                        baseDir = "Outputs/fits" )
{

  rpt <- .loadFit(  fitID, groupFolder = groupFolder, 
                    baseDir = baseDir, retObj = TRUE )


  # # TEMP - fit_CC_HMC already has rpt$map
  if( is.null(rpt$map) )
    rpt$map <- rpt$finalmap

  nP <- rpt$repOpt$nP
  # we need to make sure turnOffRbar and turnOffRinit are
  # respected. Same issue as in ms3

  if(rpt$ctlList$hypo$turnOffRbar)
  {
    rpt$data$avgRcode_p   <- rep(0,nP)
    rpt$map$lnRbar_p      <- as.factor(rep(NA,nP))
  }

  if(rpt$ctlList$hypo$turnOffRinit)
  {
    rpt$data$initRcode_p   <- rep(0,nP)
    rpt$map$lnRinit_p      <- as.factor(rep(NA,nP))
  }



  posts <- makePosts( stanfit = rpt$stanfit,
                      repObj  = rpt$repOpt,
                      data    = rpt$data,
                      pars    = rpt$pars,
                      map     = rpt$map,
                      refPts  = refPts )

  rpt$posts <- posts

  # Now we need to save it, but check if the
  # saved path is correct
  groupFolder <- here::here(baseDir,groupFolder)
  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path=groupFolder,full.names = FALSE,
                        recursive=FALSE)
  # Restrict to fit_ folders, pick the nominated simulation
  fitList <- dirList[grep(pattern="fit",x=dirList)]

  # Load fit object
  folder <- fitID
  if( is.character(fitID))
      if(!grepl(pattern="fit_",x=fitID) )
        folder <- paste("fit_",fitID,sep = "")

  if(!is.character(fitID))
    folder <- fitList[fitID]

  saveFile <- paste(folder,".rds",sep = "")
  savePath <- file.path(groupFolder,folder,saveFile)

  divTabPath <- file.path(groupFolder,folder,"divTable.csv")
  write.csv(posts$divTable, file = divTabPath)

  divGradPath <- file.path(groupFolder,folder,"divGradTable.csv")
  write.csv(posts$divGradTable, file = divGradPath)

  saveRDS(rpt, file = savePath )

  message("Posteriors complete, and reports saved to ", savePath, "\n", sep = "")
} # END redoPosts()




# makeCVtable()
# Generates a table of standard errors/CVs
# for each observation series. 
# Path specific to Herring.
makeCVtable <- function( repList )
{
  repObj  <- repList$repOpt
  nP      <- repObj$nP

  # Pull age obs error and survey obs errors
  tauAge_pg   <- sqrt( repObj$tau2Age_pg[,1:3,drop = FALSE] )
  tauSurv_pg  <- repObj$tauComb_pg[,4:5,drop = FALSE]

  # Make the table
  dataCVtab <- matrix(NA, nrow = nP, ncol = 5 )
  colnames( dataCVtab ) <- c( "tauAgeRedFB",
                              "tauAgeSnRoe",
                              "tauAgeGn",
                              "tauSurvSurface",
                              "tauSurvDive"   )

  dataCVtab <- as.data.frame(dataCVtab)
  dataCVtab[,1] <- round(tauAge_pg[,1],3)
  dataCVtab[,2] <- round(tauAge_pg[,2],3)
  dataCVtab[,3] <- round(tauAge_pg[,3],3)
  dataCVtab[,4] <- round(tauSurv_pg[,1],3)
  dataCVtab[,5] <- round(tauSurv_pg[,2],3)

  dataCVtab
} # END makeCVtable

# makeFitReport() 
# Makes an html fit report from the report
# object, based on a fit report template
makeMPReport <- function( fitID = 1, 
                          groupFolder = ".", 
                          clean = TRUE )
{
  # Need to get iscamRep from the control file/list
  # of the fitID
  
  # So, load the fit
  reports <- .loadFit(fit = fitID, groupFolder = groupFolder, retObj = TRUE)
  
  fitFolder <- here::here("Outputs",groupFolder)

  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path=fitFolder,full.names = FALSE,
                        recursive=FALSE)
  # Restrict to fit_ folders, pick the nominated simulation
  fitList <- dirList[grep(pattern="fit",x=dirList)]

  if( is.character(fitID) )
    folder <- fitList[fitList == paste("fit_",fitID,sep = "") ]
  else
    folder <- fitList[fitID]

  # Load the nominated blob
  reportsFileName <- paste(folder,".rds",sep="")
  reportsPath <- file.path(fitFolder,folder,reportsFileName)
  fitFolderPath <- here::here("Outputs",groupFolder,folder)

  # Create parameter list for rendering the document
  params <- list( rootDir= fitFolderPath,
                  RDSfile = reportsFileName)

  # Make an output file name
  outFile <- paste( "mpApplicationReport.html", sep = "")

  # Render
  rmarkdown::render(  input = here::here("Documentation","mpReportTemplate.Rmd"), 
                      output_file = outFile,
                      output_dir = fitFolderPath,
                      params = params,
                      envir = new.env(),
                      clean = clean,
                      output_format = "bookdown::html_document2" )

  # remove temporary files
  dataReportFiles <- "fitReport_files"
  # unlink(file.path(fitFolderPath,dataReportFiles), recursive = TRUE)
} # END makeFitReport()




# renameReportArrays()
# Updates the dimension names of the arrays in the 
# report lists, as a way of making plot code more efficient later.
renameReportArrays <- function( repObj = repInit, datObj = data )
{
  # Just go down the list, but first do the objects with the same names

  repNames <- names(repObj)
  datNames <- names(datObj)

  bothNames <- repNames[ repNames %in% datNames ]

  for( itemName in bothNames )
  {
    dimnames(repObj[[itemName]]) <- dimnames(datObj[[itemName]])
  }

  # Recover names
  yearNames   <- dimnames(datObj$I_pgt)[[3]]
  gearNames   <- dimnames(datObj$I_pgt)[[2]]
  stockNames  <- dimnames(datObj$I_pgt)[[1]]
  ageNames    <- dimnames(datObj$A_apgt)[[1]]

  # Some series project a year into the future
  maxYear <- as.integer(yearNames[length(yearNames)])
  projYear <- maxYear + 1
  yearNamesProj <- c( yearNames, as.character(projYear))


  # Ok, that's the data taken care of. There are still all the
  # new arrays that we created
  # Bio pars
  names(repObj$B0_p)              <- stockNames
  names(repObj$M_p)               <- stockNames
  #dimnames(repObj$initN_mult_ap)  <- list( age = ageNames, stock = stockNames)
  names(repObj$rSteepness_p)      <- stockNames
  # Selectivity parameters
  dimnames(repObj$sel_ag)           <- list( age = ageNames, gear = gearNames)
  dimnames(repObj$sel_apgt)         <- list( age = ageNames, stock = stockNames, gear = gearNames, year = yearNames)
  names(repObj$SelAlpha_g)          <- gearNames
  names(repObj$SelBeta_g)           <- gearNames
  dimnames(repObj$SelAlpha_pgt)    <- list( stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$SelBeta_pgt)     <- list( stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$epsSelAlpha_pgt) <- list( stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$epsSelBeta_pgt)  <- list( stock = stockNames, gear = gearNames, year = yearNames )
  names(repObj$sigmaSelAlpha_g)    <- gearNames
  names(repObj$sigmaSelBeta_g)     <- gearNames

  # State variables
  dimnames(repObj$B_apt)           <- list( age = ageNames, stock = stockNames, year = yearNamesProj)
  dimnames(repObj$N_apt)           <- list( age = ageNames, stock = stockNames, year = yearNamesProj)
  dimnames(repObj$vulnB_pgt)       <- list( stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$vulnB_apgt)      <- list( age = ageNames, stock = stockNames, gear = gearNames,  year = yearNames)
  dimnames(repObj$vulnN_pgt)       <- list( stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$vulnN_apgt)      <- list( age = ageNames, stock = stockNames, gear = gearNames,  year = yearNames)
  dimnames(repObj$SB_pt)           <- list( stock = stockNames, year = yearNamesProj )
  dimnames(repObj$R_pt)            <- list( stock = stockNames, year = yearNamesProj )
  dimnames(repObj$M_apt)           <- list( age = ageNames, stock = stockNames, year = yearNamesProj )
  dimnames(repObj$U_pgt)           <- list( stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$totC_pgt)        <- list( stock = stockNames, gear = gearNames, year = yearNames )
  

  # Pope's approximation
  dimnames(repObj$uAge_apgt)       <- list( age = ageNames, stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$catAge_apgt)     <- list( age = ageNames, stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$predPA_apgt)     <- list( age = ageNames, stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$ageResids_apgt)  <- list( age = ageNames, stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$tcPred_apgt)     <- list( age = ageNames, stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$tcComps_apgt)    <- list( age = ageNames, stock = stockNames, gear = gearNames, year = yearNames )

  # Observation model quantities
  dimnames(repObj$qhat_pg)         <- list( stock = stockNames, gear = gearNames )
  dimnames(repObj$tauObs_pg)       <- list( stock = stockNames, gear = gearNames )
  dimnames(repObj$z_pgt)           <- list( stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$zSum_pg)         <- list( stock = stockNames, gear = gearNames )
  dimnames(repObj$validObs_pg)     <- list( stock = stockNames, gear = gearNames )
  dimnames(repObj$SSR_pg)          <- list( stock = stockNames, gear = gearNames )
  dimnames(repObj$etaSumSq_pg)     <- list( stock = stockNames, gear = gearNames )
  dimnames(repObj$tau2Age_pg)      <- list( stock = stockNames, gear = gearNames )
  dimnames(repObj$nResids_pg)      <- list( stock = stockNames, gear = gearNames )
  dimnames(repObj$nObsAge_pg)      <- list( stock = stockNames, gear = gearNames )

  return(repObj)
}



# diagnostic plots for checking HMC results from stan
#source: https://avehtari.github.io/rhat_ess/rhat_ess.html
mcmcDiagnostics <- function(stanObj, folder)
{

  mon <- monitor(stanObj)

  samp <- as.array(stanObj)
  nPars <- dim(samp)[3]
  nChain = dim(samp)[2]

  npTable <- nuts_params(stanObj)

  divTable <- npTable |>
              filter( Parameter == "divergent__",
                      Value == 1 )



  # Make a table of parameter values at
  # divergent transitions

  
  # Check for ESS less than 100*nChain
  lowBulkESS <- subset(mon, Bulk_ESS< 100*nChain)
  if(nrow(lowBulkESS)>0)
  {
    cat(nrow(lowBulkESS),'pars less than 100*nChain for Bulk ESS \n')
    write.csv(as.data.frame(lowBulkESS),file.path(folder,'lowBulkESSpars.csv') )
  } 
    


  lowTailESS <- subset(mon, Tail_ESS< 100*nChain)
  if(nrow(lowTailESS)>0)
  {
    cat(nrow(lowTailESS),'pars less than 100*nChain for tail ESS \n')
    
    # save csv
    write.csv(as.data.frame(lowTailESS),file.path(folder,'lowTailESSpars.csv') )
  } 

  
    
  # plots for lowest ESS
  which_min_ess <- which.min(mon[, 'Tail_ESS'])

  # local efficiency plot
  pdf(file=file.path(folder,'hmcDiag_localESS.pdf'))
  print(plot_local_ess(fit = stanObj, par = which_min_ess, nalpha = 20))
  dev.off()

  pdf(file=file.path(folder,'hmcDiag_quantileESS.pdf'))
  print(plot_quantile_ess(fit = stanObj, par = which_min_ess, nalpha = 20))
  dev.off()

  # estimated how ESS changes with more iterations
  pdf(file=file.path(folder,'hmcDiag_changeESS.pdf'))
  print(plot_change_ess(fit = stanObj, par = which_min_ess))
  dev.off()

  # rank plots  
  pdf(file=file.path(folder,'hmcDiag_rankPlots.pdf'))
  print(mcmc_hist_r_scale(samp[, , which_min_ess]))
  dev.off()

  # Pairs plot
  # let's break up the pairs plots a bit
  parIdxSeq <- seq(from = 1, to = nPars, by = 15)
  parIdxSeq <- unique(c(parIdxSeq,nPars))

  nPlots <- length(parIdxSeq)-1
  parNames <- dimnames(as.array(stanObj))[[3]]
  graphics.off()
  options(warn = -1)
  for(p in 1:nPlots)
  {
    png( file = file.path(folder,paste0('hmcPairs_',p,'.png')), 
          res = 150, width = 7, height = 7, units = "in")
    # mcmc_pairs(stanObj, pars = parNames[parIdxSeq[p]:(parIdxSeq[p+1]-1)]  )
    pairs(  stanObj, pars = parNames[parIdxSeq[p]:(parIdxSeq[p+1]-1)] )
    dev.off()
  }
  options(warn = 0)
} # End mcmcDiagnostics

