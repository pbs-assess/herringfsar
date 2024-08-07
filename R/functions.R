# functions for making assessment results

# passObj()
# Replace passed objectives with a filled circle dot 
# (use unicode??)
# Inputs: X = numeric value of performance metric
#         target = threshold value
#         comp = "gt" (greater than) or "lt" (less than)
# Outputs: New value that either shows a checkmark if
#           the metric meets the target, or the
#           original value if the target is not met
# Author: SDNJ, LFR
passObj <- function( X, target = .95, comp = "gt" )
{
  out <- character(length = length(X))
  for( l in 1:length(X))
  { 
    if( length(target) > 1 )
      tar <- target[l]
    else tar <- target

    if( comp == "gt")
    {
      if( round(X,2) >= tar )
        out[l] <- paste("$\\checkmark$")
      else out[l] <- paste(X,sep = "")
    }

    if( comp == "lt")
    {
      if( round(X,2) <= tar )
        out[l] <- paste("$\\checkmark$")
      else out[l] <- paste(X,sep = "")
    }    

  }
  return(out)
} # END passObj()

# calcTAC()
# Function to calculate TAC from SISCAH model output.
# Args: repList = SISCAH model reports object from SISCAH-MP
#                 package
# Outputs: data.frame of harvest advice and management pars
# Author: Sean Cox, SDNJ, LFR
calcTAC <- function( repList = reports )
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

  out.df <- array("", dim = c(6,2))
  out.df <- as.data.frame(out.df)
  colnames(out.df) <- c("Variable","Estimate")
  out.df[,1] <- c("$\\hat{B}_{2024}",
                  "$\\hat{B}_{0}",
                  "LCP",
                  "OCP",
                  "THR",
                  "TAC" )

  out.df[,2] <- round(c(Bhat, estB0, 0.3*estB0, 0.6*estB0, targU, Q),3)

  return(out.df)

} # END calcTAC()

# stock_status_text()
# Function to automate stock status reporting from
# a table off ref pts and model output. Currently
# works off MS3 blob but should be updated to be based
# on a weighted posterior (higher sample size)
# Inputs: refPtsTab = table of ensemble model ref pts
#         history   = posterior biomass history
stock_status_text <- function(  refPtsTab = ensRefPtsTable,
                                parTab = ensParTable,
                                history = mpBlobList[[2]],
                                fYear = 1951, lYear = 2023  )
{

  yrs <- fYear:lYear
  lastTdx <- length(yrs)
  goodReps <- history$goodReps
  SB_t <- apply(history$om$SB_ispt[goodReps,1,1,], FUN = median, MARGIN = 2)
  SB_T <- round(SB_t[lastTdx],3)
  SB_Tm1 <- round(SB_t[lastTdx-1],3)

  PBTGtLRP <- parTab$PBTGtLRP

  B0 <- round(refPtsTab$B0,3)

  paste0( "Estimated unfished spawning biomass $SB_0$ is ", B0,
  " kt, and the LRP of $0.3 \\cdot SB_0$ is ", round(0.3 * B0,3) ," kt (posterior medians). Compared to last year, estimated spawning biomass in 2023 $SB_{2023}$ decreased from ",  
  SB_Tm1, " to ", SB_T, " kt (posterior median), and is equivalent to ", 
  100 * round(SB_T/B0,3), 
  " \\% of $SB_0$ (Tables XX & XX). Spawning biomass in 2023 is estimated to be above the LRP with a ", 
  100 * PBTGtLRP, " \\% probability (Table XX).")

} # END stock_status_text


# proj_biomass_text()
# Function to automate biomass forecast reporting from
# model output. Currently works off MLEs but should be 
# updated to be based on posteriors (higher sample size)
# Inputs: MPfit   = biomass history
proj_biomass_text <- function(  mpFit = fit_maxTHR0.14,
                                fYear = 1951, lYear = 2023,
                                B0 = ensRefPtsTable$B0  )
{

  yrs <- fYear:(lYear + 1)
  lastTdx <- length(yrs)-1
  projTdx <- length(yrs)
  
  SB_t <- (mpFit$repOpt$SB_pt[1,])
  SB_T <- round(SB_t[lastTdx],3)
  SB_forecast <- round(SB_t[projTdx],3)

  B0 <- round(B0,3)

  paste0( "In the absence of fishing, spawning biomass in ", lYear+1, 
    "$SB_{", lYear+1,"}$ is estimated to be ", SB_forecast, 
    " kt (maximum likelihood estimate; Table \\@ref(tab:TACtable)). Spawning biomass in ", 
    lYear + 1, " is forecast to be below the LRP of $0.3SB_0$ (",
    round(0.3*B0), " kt) with a X\\% probability, in the absence of fishing (Table XX and Figure XX).")

} # END proj_biomass_text


# makeModelHistTable()
# Generates a model historical biomass table
# from an input MS3 simulation object for plotting
# as a 4-panel figure in the Herring FSAR.
# Args:
#   obj     = MS3 blob object
#   fYear   = first year of model history
#   lYear   = last year of model history
#   B0      = externally estimated unfished biomass
#   USR     = externally estimated USR
#   TAC     = MP determined TAC
#   maxTHR  = maximum target harvest rate
# Value:
#   out.df  = Table of catch, biomass, fishing/natural 
#             mortality, and recruitment for plotting.
makeModelHistTable <- function( obj,
                                fYear = 1951,
                                lYear = 2023,
                                B0 = 90.23,
                                USR = 67.31,
                                TAC = 11.568,
                                Uref = 0.14 )
{

  # Pull good replicates
  goodReps  <- which(obj$goodReps)

  # Years
  yrs       <- fYear:lYear
  tIdx      <- 1:length(yrs)
  
  # Catch
  C_qt      <- apply(obj$om$C_ispt[goodReps,1,1,tIdx], FUN = quantile, probs = c(0.025,0.5,.975),MARGIN = 2, na.rm = TRUE)
  # Biomass
  SB_qt     <- apply(obj$om$SB_ispt[goodReps,1,1,tIdx], FUN = quantile, probs = c(0.025,0.5,.975),MARGIN = 2, na.rm = TRUE)
  # Harvest rates
  U_ispt    <- obj$om$C_ispt/(obj$om$C_ispt + obj$om$SB_ispt)
  U_qt      <- apply(U_ispt[goodReps,1,1,tIdx], FUN = quantile, probs = c(0.025,0.5,.975),MARGIN = 2, na.rm = TRUE)
  # Recruitment
  R_qt      <- apply(obj$om$R_ispt[goodReps,1,1,tIdx], FUN = quantile, probs = c(0.025,0.5,.975),MARGIN = 2, na.rm = TRUE)
  # Natural mortality
  M_qt      <- apply(obj$om$M_iaxspt[goodReps,2,1,1,1,tIdx], FUN = quantile, probs = c(0.025,0.5,.975),MARGIN = 2, na.rm = TRUE)

  # Data.frame of data for plotting.
  out.df    <- data.frame(  Year    = yrs, 
                            Catch   = C_qt[2,],
                            TAC     = mean(C_qt[2,]),
                            SSB_med = SB_qt[2,],
                            SSB_min = SB_qt[1,],
                            SSB_max = SB_qt[3,],
                            LRP     = 0.3*B0,
                            USR     = USR,
                            F_med   = U_qt[2,],
                            F_min   = U_qt[1,],
                            F_max   = U_qt[3,],
                            F_lim   = Uref,
                            M_med   = 1 - exp(-M_qt[2,]),
                            M_max   = 1 - exp(-M_qt[3,]),
                            M_min   = 1 - exp(-M_qt[1,]),
                            R_med   = R_qt[2,],
                            R_max   = R_qt[3,],
                            R_min   = R_qt[1,]
                          )

  out.df
} # END makeModelHistTable()


