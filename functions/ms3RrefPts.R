# --------------------------------------------------------------------------
# refPts.R
# 
# Takes parameters from a report object and calculates
# reference points and curves under an age structured/YPR formulation. 
# Called from within runHierSCAL() after runnings an AM.
# 
# 
# Author: Samuel Johnson
# Date: March 7, 2019
# 
# --------------------------------------------------------------------------

# calcRefPts()
# Calculates biological reference points based on given biological parameters 
# assuming a delay difference biological model.
# inputs:   obj = list of biological parameters
# ouputs:   refPts = list() of reference points
calcRefPts <- function( obj )
{
  # Calculate selectivity
  # obj <- .calcSel_xsp(obj, fleetIdx = 2)
  nS <- obj$om$nS
  nP <- obj$om$nP
  nF <- obj$om$nF

  nSpec   <- nS

  # temporarily use calcPerRecruit to recalc R0
  obj$relF_spf    <- array(0,dim = c(nS,nP,nF))    
  tmp             <- .calcPerRecruit( f = 0, obj = obj )
  yprList         <- tmp$yprList



  # Calculate R0_sp and totB
  h_sp          <- obj$om$h_sp
  B0_sp         <- obj$om$B0_sp

  obj$R0_sp     <- array(0,dim = c(nS,nP))
  obj$totB0_sp  <- array(0,dim = c(nS,nP))

  obj$M0_xsp    <- yprList$Meq_xsp

  obj$R0_sp     <- B0_sp / yprList$ssbpr_sp
  obj$totB0_sp  <- obj$R0_sp * yprList$totbpr_sp

  lastB_sp      <<- obj$totB0_sp

  # Beverton-Holt a/b parameters
  obj$rec.a_sp  <- 4.*h_sp*obj$R0_sp/(B0_sp*(1.-h_sp))
  obj$rec.b_sp  <- (5.*h_sp-1.)/(B0_sp*(1.-h_sp))

  # Need to calculate relF_spf
  # Calculate relative F

  alloc_spf       <- obj$om$alloc_spf
  commIdx_spf     <- alloc_spf > 0

  # In case there's only one fishery...
  commIdx_spf[alloc_spf == 1] <- FALSE
  obj$relF_spf[alloc_spf == 1] <- 1 

  if(any(commIdx_spf))
  {
    nComm           <- sum(commIdx_spf)
    optRelF         <- optim( par = rep(-2,nComm), fn = .getRelFs,
                              method = "BFGS", control=list(maxit = 50),
                              f = mean(obj$om$M_xsp), obj = obj,
                              commIdx_spf = commIdx_spf,
                              allocVar = obj$allocVar )

    # Overwrite relF_pf
    # Overwrite relF_pf
    idx <- 0
    for( s in 1:nS )
    {
      for(p in 1:nP)
      {
        fIdx <- which(commIdx_spf[s,p,])
        nF <- length(fIdx)
        obj$relF_spf[s,p,fIdx] <- exp(optRelF$par[idx + 1:nF])
      }

      idx <- idx + nF 
    }

    # Normalise
    for(s in 1:nS)
      for( p in 1:nP)
          obj$relF_spf[s,p,] <- obj$relF_spf[s,p,]/sum(obj$relF_spf[s,p,])
  }

  # Calculate generation time
  obj$genTime_sp <- .calcGenTime(obj)

  lastB_sp <<- obj$totB0_sp

  # Calculate reference curves
  refCurves <- .calcRefCurves( obj )
  
  # First, let's just do Fmsy reference points
  FmsyRefPts <- .getFmsy_sp(  obj = obj, 
                              refCurves = refCurves )

  

  
  obj$refPts <- list()
  obj$refPts$refCurves    <- refCurves
  obj$refPts$FmsyRefPts   <- FmsyRefPts

  # message("Fmsy = ", FmsyRefPts$Fmsy_sp, "\n")
  # browser()
  # message("Bmsy = ", FmsyRefPts$S_sp, "\n")

  # totB0_sp[1:nS,1:nP] <- refCurves$totBeq_spf[,,1]

  
  # Get survivorship
  obj$refPts$Surv_axsp    <- tmp$Surv_axsp
  obj$refPts$ssbpr_sp     <- yprList$ssbpr_sp
  obj$refPts$R0_sp        <- obj$R0_sp
  obj$refPts$totB0_sp     <- obj$totB0_sp
  obj$refPts$totbpr_sp    <- yprList$totbpr_sp
  obj$refPts$rec.a_sp     <- obj$rec.a_sp
  obj$refPts$rec.b_sp     <- obj$rec.b_sp
  obj$refPts$B0_sp        <- refCurves$Beq_spf[,,1,drop = FALSE]
  obj$refPts$M_xsp        <- obj$om$M_xsp
  obj$refPts$m1_sp        <- obj$om$m1_sp
  obj$refPts$Mb_sp        <- obj$om$Mb_sp
  obj$refPts$h_sp         <- h_sp
  obj$refPts$relF_spf     <- obj$relF_spf

  

  return(obj$refPts)

} # END calcRefPts()

# .getRelFs()
# Function used for determining relative fishing
# mortality rates that meet a nominated allocation
# among gears.
.getRelFs <- function( lnfg, obj, f=0, commIdx_spf, allocVar = "Biomass" )
{
  nS <- obj$om$nS
  nP <- obj$om$nP
  nF <- obj$om$nF

  # Overwrite relF_pf
  idx <- 0
  for( s in 1:nS )
  {
    for(p in 1:nP)
    {
      fIdx <- which(commIdx_spf[s,p,])
      nF <- length(fIdx)
      obj$relF_spf[s,p,fIdx] <- exp(lnfg[idx + 1:nF])
    }

    idx <- idx + nF 
  }


  # Normalise
  for(s in 1:nS)
    for( p in 1:nP)
        obj$relF_spf[s,p,] <- obj$relF_spf[s,p,]/sum(obj$relF_spf[s,p,])

  lastB_sp <<- obj$totB0_sp

  obj$relF_spf[!is.finite(obj$relF_spf)] <- 0
  # Calculate legal YPR
  tmp         <- .calcPerRecruit( f=f, obj )$yprList 

  funcVal <- 0
  for(s in 1:nS)
    for(p in 1:nP)
    {
      commGears <- which(obj$om$alloc_spf[s,p,] > 0)
      # Switch here between egg yield and bio/ponded fish yield
      if(allocVar == "Eggs")
        prop <- tmp$epr_spf[s,p,commGears]/sum(tmp$epr_spf[s,p,commGears])
      
      if( allocVar == "Biomass")
        prop <- tmp$ypr_spf[s,p,commGears]/sum(tmp$ypr_spf[s,p,commGears])

      # cat("s = ",s,"; propYield = ",prop,"\n")
      funcVal <- funcVal + sum((log(prop)-log(obj$om$alloc_spf[s,p,commGears]))^2.)
      # cat("f = ", funcVal,"\n")
    }
        
  
  funcVal
}


# .calcRefCurves()
# Calculates equilibrium curves of equilbrium biomass, numbers,
# yield and recruitment as a function of input fishing mortality rates
# inputs:   obj = list of biological parameters
# ouputs:   refCurves = list() of reference curves (vectors)
.calcRefCurves <- function( obj, nFs = 200 )
{
  # First, compute max F (tolerance of 1e-5)
  nT   <- dim(obj$om$qF_spft)[4]
  maxF <- 10 * max(obj$om$M_xsp,na.rm = T)

  if( obj$condModel == "hierSCAL")
    maxE <- max( maxF / obj$om$qF_spft[,,2,nT])

  # We're going to need to fill each species' ref curves,
  # so labeling and dimensions are needed
  nS          <- obj$om$nS
  nP          <- obj$om$nP
  nA          <- obj$om$nA
  nX          <- obj$om$nX

  specNames   <- c("Herring","Hake")[1:nS]
  stockNames  <- dimnames(obj$M_apt)[[2]]  

  f <- seq( from = 0.0, to = maxF, length = nFs )

  if( obj$condModel == "hierSCAL")
    e <- seq( from = 0.0, to = maxE, length = nFs )



  # Create matrices to hold Recruitment reference curve, name rows and copy
  # for each of Beq, Neq and Yeq
  Req_spf      <- array( NA,  dim = c(nS, nP, nFs),
                              dimnames = list(  species = specNames,
                                                stock = stockNames,
                                                F = f ) )

  surv_axspf  <- array( NA, dim = c(nA,nX,nS,nP,nFs) )
  Z_axspf     <- array( NA, dim = c(nA,nX,nS,nP,nFs) )
  Meq_xspf    <- array( NA, dim = c(nX,nS,nP,nFs) )

  Beq_spf       <- Req_spf
  totBeq_spf    <- Req_spf
  Yeq_spf       <- Req_spf
  ypr_spf       <- Req_spf
  ssbpr_spf     <- Req_spf
  Ueq_spf       <- Req_spf
  Yeq_spf       <- Req_spf
  EYeq_spf      <- Req_spf
  ypr_spf       <- Req_spf
  epr_spf       <- Req_spf
  ssbpr_spf     <- Req_spf
  matbpr_spf    <- Req_spf
  Ueq_spf       <- Req_spf
  Udead_spf     <- Req_spf


  # Loop and fill
  for( i in 1:length(f) )
  {
    tmp               <- .calcEquil( f = f[i], obj = obj )
    Req_spf[,,i]      <- tmp$Req_sp
    Beq_spf[,,i]      <- tmp$Beq_sp
    # expBeq_spf[,,i]   <- tmp$expBeq_sp
    totBeq_spf[,,i]   <- tmp$totBeq_sp
    Yeq_spf[,,i]      <- tmp$Yeq_sp
    EYeq_spf[,,i]     <- tmp$EYeq_sp
    ypr_spf[,,i]      <- tmp$ypr_sp
    epr_spf[,,i]      <- tmp$epr_sp
    ssbpr_spf[,,i]    <- tmp$ssbpr_sp
    Ueq_spf[,,i]      <- tmp$Ueq_sp
    surv_axspf[,,,,i] <- tmp$surv_axsp
    Meq_xspf[,,,i]    <- tmp$Meq_xsp
    Z_axspf[,,,,i]    <- tmp$Z_axsp
  }


  # Save F based ref points
  refCurves <- list()
    refCurves$F           <- f
    refCurves$ypr_spf     <- ypr_spf
    refCurves$Req_spf     <- Req_spf
    refCurves$Beq_spf     <- Beq_spf
    # refCurves$expBeq_spf  <- expBeq_spf
    refCurves$totBeq_spf  <- totBeq_spf
    refCurves$Yeq_spf     <- Yeq_spf
    refCurves$EYeq_spf    <- EYeq_spf
    refCurves$Ueq_spf     <- Yeq_spf
    refCurves$Meq_xspf    <- Meq_xspf
    refCurves$surv_axspf  <- surv_axspf
    refCurves$Z_axspf     <- Z_axspf


  return( refCurves )
} # END .calcRefCurves

# .calcEquil()
# Calculates equilibrium Biomass, Yield and Recruitment 
# with respect to given fishing mortality rates F
# inputs:   f = input fishing mortality rate
#           obj = list of biological parameters
# ouputs:   equil = list() of equilibrium biomass, yield and recruitment
.calcEquil <- function( f = 0, obj, type = "fmort" )
{
  nS    <- obj$om$nS
  nMICE <- obj$om$nMICE
  nP    <- obj$om$nP
  A_s   <- obj$om$A_s


  # Now calculate eqbm recruitment at given f value
  tmp <- .calcPerRecruit( f = f, obj = obj, type = type )
  yprList <- tmp$yprList



  recruits_sp <- ( obj$rec.a_sp * yprList$ssbpr_sp - 1) / (obj$rec.b_sp * yprList$ssbpr_sp)

  
  equil <- list()
    equil$Req_sp      <- recruits_sp
    equil$Beq_sp      <- recruits_sp * yprList$ssbpr_sp
    # equil$expBeq_sp  <- recruits_sp * yprList$expbpr_sp
    equil$totBeq_sp   <- recruits_sp * yprList$totbpr_sp
    equil$Yeq_sp      <- recruits_sp * yprList$ypr_sp
    equil$EYeq_sp     <- recruits_sp * yprList$epr_sp
    equil$ypr_sp      <- yprList$ypr_sp
    equil$epr_sp      <- yprList$epr_sp
    equil$ssbpr_sp    <- yprList$ssbpr_sp
    equil$Ueq_sp      <- equil$Yeq_sp / (equil$Beq_sp + equil$Yeq_sp)
    equil$surv_axsp   <- tmp$Surv_axsp
    equil$Meq_xsp     <- yprList$Meq_xsp
    equil$Z_axsp      <- yprList$Z_axsp

  return(equil)
}


.calcSel_xsp <- function( obj, fleetIdx = 2)
{
  # Model dimensions
  nS    <- obj$nS
  nP    <- obj$nP
  nX    <- obj$nX
  A_s   <- obj$A_s
  nA    <- obj$nA
  nL    <- obj$nL

  # Probability of being a given length-at-age
  probLenAge_laspx <- obj$probLenAge_laspx

  # Selectivity - this is mean over each fleet's 
  # time period, not time-varying
  # Might be useful to take the group mean for the comm fleet...
  xSel50_sp   <- obj$xSel50_spf[ , , fleetIdx, drop = FALSE ]
  xSel95_sp   <- obj$xSel95_spf[ , , fleetIdx, drop = FALSE ]
  xSelStep_sp <- obj$xSelStep_spf[ , , fleetIdx, drop = FALSE ]

  # Harcode for comm.mod for now,
  # and length based selectivity
  selLen_lsp  <- array(NA, dim = c(nL,nS,nP) )
  selAge_axsp <- array(NA, dim = c(nA,nX,nS,nP) )

  # Loop over species and stocks, calculate
  # selLen so we can get selAge
  for( s in 1:nS )
    for(p in 1:nP )
    {
      selLen_lsp[,s,p] <- (1 + exp(-(1:nL - xSel50_sp[s,p,])/xSelStep_sp[s,p,]/log(19)))^(-1)
      for( x in 1:nX)
      {
        for( a in 1:nA )
        {
          selAge_axsp[a,x,s,p] <- sum(probLenAge_laspx[,a,s,p,x] * selLen_lsp[,s,p])
        }
        selAge_axsp[,x,s,p] <- selAge_axsp[,x,s,p] / max(selAge_axsp[,,s,p],na.rm = T)
      }
    }

  obj$selLen_lsp <- selLen_lsp
  obj$selAge_axsp <- selAge_axsp

  obj
}

# .calcPerRecruit
# Purpose:     Calculate all equilibrium per-recruit quantities of interest for an
#              input fishing mortality.
# Parameters:  f=scalar input fishing mortality rate; obj=list of all operating
#              model parameters.
# Returns:     a list with equilibrium quantities - (i) spawning stock biomass-per-recruit
#              and (ii) yield-per-recruit (ypr)
# Source:      S.P. Cox, modified for hierSCAL by SDNJ
.calcPerRecruit <- function( f, obj, type = "fmort" )
{

  # Compute eqbm spawning biomass per recruit for
  # given f and species/stock pars
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nF      <- obj$om$nF
  A_s     <- obj$om$A_s
  nA      <- obj$om$nA
  nL      <- obj$om$nL
  nX      <- obj$om$nX
  M_xsp   <- obj$om$M_xsp
  Mb_sp   <- obj$om$Mb_sp
  m1_sp   <- obj$om$m1_sp
  nT      <- obj$om$nT
  tMP     <- obj$om$tMP
  qF_spf  <- obj$om$qF_spft[,,,nT]

  nSpec   <- nS


  # Life history schedules
  matAge_asp        <- obj$om$matAge_asp
  wtAge_axsp        <- obj$om$meanWtAge_axsp

  # Spawn timing
  spawnTiming_s     <- obj$om$spawnTiming_s

  # Recruitment pars
  if(f > 0)
  {
    rec.a_sp <- obj$rec.a_sp
    rec.b_sp <- obj$rec.b_sp
    totB0_sp <- obj$totB0_sp
  }
  


  # Selectivity 
  selAge_axspf           <- array(NA, dim = c(nA,nX,nS,nP,nF))
  selAge_axspf[1:nA,1:nX,1:nS,1:nP,1:nF]  <- obj$om$sel_axspft[1:nA,1:nX,1:nS,1:nP,1:nF,tMP-1]

  # Fishing mortality
  relF_spf  <- obj$relF_spf
  relF_spf[relF_spf < 1e-5] <- 0
  f_spf <- f * relF_spf

  if(obj$densityDepM == 1)
  {
    if(f == 0 )
    {
      Meq_sp    <- Mb_sp + exp(-m1_sp)
      lastB_sp  <<- obj$totB0_sp
    }
    if( f > 0 )
    {   
      lnB_sp <- log(lastB_sp)
      if( nS > 1 | nP > 1)
      {
        opt <- optim( par = lnB_sp,
                      f = solveForMeq,                       
                      ff = f, obj = obj,
                      fit = TRUE )
        lnB_sp[1:nS,1:nP] <- opt$par 
      }

      if(nS == 1 & nP == 1)
      {
        opt <- optimize(  f = solveForMeq, 
                          interval = c(log(0.5*lastB_sp), lnB_sp),
                          ff = f, obj = obj, fit = TRUE )
        lnB_sp[1:nS,1:nP] <- opt$minimum
      }

      
      lastB_sp[1:nS,1:nP] <<- exp(lnB_sp)
      
      

      Meq_sp <- solveForMeq( lnB_sp = lnB_sp, obj = obj,ff = f, fit = FALSE)
    }

    # cat("lastB_sp = ", lastB_sp, "\n", sep = "")

    M_xsp[1,,] <- Meq_sp
  }

  # message("lastB = ", lastB_sp, "\n")

  # Zero indexing
  juveMage <- obj$juveMage


  # Compute Z_asp
  Z_axsp    <- array( NA, dim = c(nA,nX,nS,nP))
  Surv_axsp <- array( NA, dim = c(nA,nX,nS,nP))
  Surv_axsp[1,,,] <- 1/nX

  
  for( x in 1:nX )
    for( s in 1:nS )
      for( p in 1:nP )
      {
        
        if( juveMage > 0 )
          Z_axsp[1:(juveMage),x,s,p] <- obj$Mjuve_p[p]

        Z_axsp[(juveMage+1):A_s[s],x,s,p] <- M_xsp[x,s,p]
        
        for( a in 1:A_s[s])
        {
          for( f in 1:nF)
            Z_axsp[a,x,s,p] <- Z_axsp[a,x,s,p] + selAge_axspf[a,x,s,p,f] * f_spf[s,p,f]
          
          if( a > 1 )
            Surv_axsp[a,x,s,p] <- Surv_axsp[a-1,x,s,p] * exp( -Z_axsp[a-1,x,s,p])
          if( a == A_s[s])
            Surv_axsp[a,x,s,p] <- Surv_axsp[a,x,s,p] / (1 - exp(-Z_axsp[a,x,s,p]))
        }
      }


  # Calculate yield-per-recruit, ssb per recruit, and exp biomass per recruit
  # by using survival array
  ssbpr_asp     <- array( NA, dim = c(nA,nS,nP) )
  matbpr_asp    <- array( NA, dim = c(nA,nS,nP) )
  expbpr_axspf  <- array( NA, dim = c(nA,nX,nS,nP,nF) )
  totbpr_axsp   <- array( NA, dim = c(nA,nX,nS,nP) )
  C_axspf       <- array(0,   dim = c(nA,nX,nS,nP,nF))
  P_axspf       <- array(0,   dim = c(nA,nX,nS,nP,nF))
  E_axspf       <- array(0,   dim = c(nA,nX,nS,nP,nF))


  for( s in 1:nS )
    for( p in 1:nP)
    {
      ssbpr_asp[,s,p]  <- Surv_axsp[,nX,s,p] * 
                          wtAge_axsp[,nX,s,p] * 
                          matAge_asp[,s,p]

      # Keep track of separate ssbpr and mature biomass. SSB is
      # fish that are used in the SR relationship, while mature biomass
      # includes spawners whose eggs are harvested by SOK (or predators)
      ssbpr_asp[,s,p] <- ssbpr_asp[,s,p] * exp(-spawnTiming_s[s] * Z_axsp[,nX,s,p])        
      matbpr_asp[,s,p] <- ssbpr_asp[,s,p]
      

      for( x in 1:nX )
      {
        for( f in 1:nF)
        {
          # Catch - total ponded fish for SOK fleets
          C_axspf[,x,s,p,f] <-  Surv_axsp[,x,s,p] * wtAge_axsp[,x,s,p] * 
                                selAge_axspf[,x,s,p,f] * f_spf[s,p,f] * 
                                (1 - exp(-Z_axsp[,x,s,p]))/Z_axsp[,x,s,p]

          # Egg yield
          E_axspf[,x,s,p,f]   <- C_axspf[,x,s,p,f] * matAge_asp[,s,p] * 0.5 * 200 * 0.35
          # Convert eggs to SOK and overwite catch, 
          # record "catch" based on F as ponded fish
          if( obj$om$fleetType_f[f] == 2 )
          {
            P_axspf[,x,s,p,f] <-  C_axspf[,x,s,p,f]

            # C_axspf[,x,s,p,f] <-  E_axspf[,x,s,p,f] * 
            #                       obj$om$gamma_f[f]

            # Return ponded fish to mature biomass, but
            # they don't contribute to spawners
            matbpr_asp[,s,p]  <- matbpr_asp[,s,p] + exp(-0.315)* P_axspf[,x,s,p,f]
          }

          expbpr_axspf[,x,s,p,f] <- Surv_axsp[,x,s,p] * 
                                    selAge_axspf[,x,s,p,f] * 
                                    wtAge_axsp[,x,s,p]

          
        }

        
        totbpr_axsp[,x,s,p] <- Surv_axsp[,x,s,p] * wtAge_axsp[,x,s,p]
      }

    }

  # Replace NAs with 0 (unmodeled ages)
  # Z_axsp[is.na(Z_axsp)] <- 0
  # Surv_axsp[is.na(Surv_axsp)] <- 0
  # C_axsp[is.na(C_axsp)] <- 0
  
  # Compute ratios
  ypr_sp        <- apply( X = C_axspf, FUN = sum, MARGIN = c(3,4),na.rm = T)
  ypr_spf       <- apply( X = C_axspf, FUN = sum, MARGIN = c(3:5),na.rm = T)
  epr_sp        <- apply( X = E_axspf, FUN = sum, MARGIN = c(3,4),na.rm = T)
  epr_spf       <- apply( X = E_axspf, FUN = sum, MARGIN = c(3:5),na.rm = T)
  ssbpr_sp      <- apply( X = ssbpr_asp, FUN = sum, MARGIN = c(2,3), na.rm = T )
  matbpr_sp     <- apply( X = matbpr_asp, FUN = sum, MARGIN = c(2,3), na.rm = T )
  expbpr_spf    <- apply( X = expbpr_axspf, FUN = sum, MARGIN = 3:5, na.rm = T )  
  totbpr_sp     <- apply( X = totbpr_axsp[2:nA,,,,drop = FALSE], FUN = sum, MARGIN = c(3,4), na.rm = T )  


  # if(any(expbpr_sp - ypr_sp < 0))

  # browser()  

  # compile output list
  yprList <- list(  ssbpr_sp    = ssbpr_sp,
                    ypr_sp      = ypr_sp,
                    ypr_spf     = ypr_spf,
                    epr_sp      = epr_sp,
                    epr_spf     = epr_spf,
                    expbpr_spf  = expbpr_spf,
                    totbpr_sp   = totbpr_sp,
                    matbpr_sp   = matbpr_sp,
                    Meq_xsp     = M_xsp,
                    Z_axsp      = Z_axsp )

  obj$yprList <- yprList
  obj$Surv_axsp <- Surv_axsp

  return(obj)
}

# solve for density dependent M
solveForMeq <- function( lnB_sp = log(totB0_sp), obj, ff, fit = TRUE )
{
  # Compute eqbm spawning biomass per recruit for
  # given f and species/stock pars
  # Compute eqbm spawning biomass per recruit for
  # given f and species/stock pars
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nF      <- obj$om$nF
  A_s     <- obj$om$A_s
  nA      <- obj$om$nA
  nL      <- obj$om$nL
  nX      <- obj$om$nX
  M_xsp   <- obj$om$M_xsp
  Mb_sp   <- obj$om$Mb_sp
  m1_sp   <- obj$om$m1_sp
  nT      <- obj$om$nT
  tMP     <- obj$om$tMP
  qF_spf  <- obj$om$qF_spft[,,,nT]


  totB0_sp <- obj$totB0_sp

  # Life history schedules
  matAge_asp        <- obj$om$matAge_asp
  wtAge_axsp        <- obj$om$meanWtAge_axsp

  # Spawn timing
  spawnTiming_s     <- obj$om$spawnTiming_s

  # Selectivity 
  selAge_axspf           <- array(NA, dim = c(nA,nX,nS,nP,nF))
  selAge_axspf[1:nA,1:nX,1:nS,1:nP,1:nF]  <- obj$om$sel_axspft[1:nA,1:nX,1:nS,1:nP,1:nF,tMP-1]


  # Fishing mortality
  relF_spf  <- obj$relF_spf
  relF_spf[relF_spf < 1e-3] <- 0
  f_spf <- ff * relF_spf

  # DepM model pars
  Mb_sp         <- obj$om$Mb_sp
  m1_sp         <- obj$om$m1_sp
  
  juveMage      <- obj$juveMage
  Mjuve_p       <- obj$Mjuve_p

  # Recruitment pars
  if(ff > 0)
  {
    rec.a_sp <- obj$rec.a_sp
    rec.b_sp <- obj$rec.b_sp
  }

  
  # message("totB0 = ", totB0_sp, "\n")
  initTotBeq_sp <- exp(lnB_sp)
  Meq_sp        <- Mb_sp + exp(-m1_sp * initTotBeq_sp/totB0_sp)
  
  if(ff == 0)
    return(Meq_sp)

  ssbpr_sp <- array(0,dim = c(nS,nP))
  totbpr_sp <- array(0,dim = c(nS,nP))

  # Compute Z_asp
  Z_axsp    <- array( NA, dim = c(nA,nX,nS,nP))
  Surv_axsp <- array( NA, dim = c(nA,nX,nS,nP))
  Surv_axsp[1,,,] <- 1/nX

  
  for( x in 1:nX )
    for( s in 1:nS )
      for( p in 1:nP )
      {
        
        if( juveMage > 0 )
          Z_axsp[1:(juveMage),x,s,p] <- Mjuve_p[p]

        Z_axsp[(juveMage+1):A_s[s],x,s,p] <- Meq_sp[s,p]
        
        for( a in 1:A_s[s])
        {
          for( fIdx in 1:nF)
            Z_axsp[a,x,s,p] <- Z_axsp[a,x,s,p] + selAge_axspf[a,x,s,p,fIdx] * f_spf[s,p,fIdx]
          
          if( a > 1 )
            Surv_axsp[a,x,s,p] <- Surv_axsp[a-1,x,s,p] * exp( -Z_axsp[a-1,x,s,p])
          if( a == A_s[s])
            Surv_axsp[a,x,s,p] <- Surv_axsp[a,x,s,p] / (1 - exp(-Z_axsp[a,x,s,p]))
        }
      }


  # Calculate ssbpr and totbpr
  for( p in 1:nP )
  {
    for( s in 1:nS )     
    {
      ssbpr_sp[s,p]  <- sum(Surv_axsp[,,s,p] * wtAge_axsp[,,s,p] * matAge_asp[,s,p] * exp(-spawnTiming_s[s] * Z_axsp[,,s,p]))
      totbpr_sp[s,p] <- sum(Surv_axsp[(juveMage+1):nA,,s,p] * wtAge_axsp[(juveMage+1):nA,,s,p])    
    } 
    
  }

  guessR_sp     <- initTotBeq_sp/totbpr_sp
  guessSB_sp    <- guessR_sp * ssbpr_sp
  guessR2_sp    <- (rec.a_sp * ssbpr_sp - 1) / (rec.b_sp * ssbpr_sp)
  
  guessTotB_sp  <- guessR2_sp * totbpr_sp

  # guessR22_p <- rec.a_sp[1,] * guessSB_p / (1 + rec.b_sp[1,] * guessSB_p)

  # message( "totbpr = ", totbpr_sp, "\n")
  # message( "ssbpr = ", ssbpr_sp, "\n")
  # message( "guessB = ", initTotBeq_sp, "\n")
  # message( "guessB2 = ", guessTotB_sp, "\n")
  # message( "guessB/tbpr = ", guessR_sp, "\n")
  # message( "ssbpr Rec = ", guessR2_sp, "\n")
  # message( "ssbpr rec * totbpr = ", guessTotB_sp, "\n")
   
  if( fit )
  {
    resid <- sum((log(guessR2_sp) - log(guessR_sp))^2)
    # message("optB f = ", resid, "\n")
    return(resid)
  }

  return(Meq_sp)

}


# # Calculates recruitment parameters, and equilibrium unfished
# # numbers and recruitment.
# .calcRecPars <- function( obj )
# {
#   # Calculate eqbm parameters
#   # Survival

    
#   # Beverton-Holt a parameters
#   rec.a <- 4.*obj$rSteepness*R0/(obj$B0*(1.-obj$rSteepness))
#   rec.b <- (5.*obj$rSteepness-1.)/(obj$B0*(1.-obj$rSteepness))

#   # Now return everything in a list
#   recList <- list(  S0 = S0,
#                     wbar0 = wbar0,
#                     N0 = N0,
#                     R0 = R0,
#                     rec.a = rec.a,
#                     rec.b = rec.b  )

#   return(recList)
# }

# .getEmsy_p()
# Calculates effort based reference curves
# for each stock and species, then computes
# multispecies optimum yield and effort
# values for each stock area. Returns
# arrays for plotting of yield curves.
.getEmsy_p <- function( obj, refCurves, FmsyRefPts )
{
  nS  <- obj$nS
  nP  <- obj$nP

  Eseq <- refCurves$E

  .getEmsy <- function( yieldCurve, E = Eseq )
  {
    minE <- 0
    maxE <- max(E)

    # Now check if yieldCurve goes negative anywhere
    if( any(yieldCurve < 0) )
    {
      minE <- E[min(which(yieldCurve >= 0))]
      maxE <- E[max(which(yieldCurve >= 0))]
    }

    eySplineFun <- splinefun( x=E, y=yieldCurve )  

    # Find stat point for Fmsy
    Emsy <- try( uniroot( f = eySplineFun, interval = c(minE, maxE),
                          deriv = 1 )$root )
    if(class(Emsy) == "try-error")
    {
      browser(cat("uniroot failed\n\n"))
      Emsy <- 0
    }

    Emsy <- min(Emsy, maxE)

    return(Emsy)
  }

  # calculate Fmsy for each species/stock
  Emsy_sp <-  apply(  X = refCurves$Yeq_spe, FUN = .getEmsy,
                      MARGIN = c(1,2) )

  Yeq_spe <- refCurves$Yeq_spe
  Yeq_spe[Yeq_spe < 0] <- NA

  # Now we want the total stock area yield curve
  Yeq_pe  <- apply( X  = refCurves$Yeq_spe, FUN = sum, MARGIN = c(2,3), 
                    na.rm = T )
  Yeq_pe[Yeq_pe == 0] <- NA
  Yeq_pe[,1] <- 0
  Emsy_p  <- apply( X = Yeq_pe, FUN = .getEmsy, MARGIN = 1)

  # Now create a matrix to hold the species/stock values
  # on each curve
  spMat       <- matrix(0, nrow = nS, ncol = nP)
  dimnames(spMat) <- dimnames(Emsy_sp)

  EmsyRefPts  <- list(  yprEmsy_sp    = spMat,
                        YeqEmsy_sp    = spMat,
                        BeqEmsy_sp    = spMat,
                        ReqEmsy_sp    = spMat,
                        expBeqEmsy_sp = spMat,
                        Umsy_sp       = spMat )

  EmsyMSRefPts <- list( yprEmsy_sp    = spMat,
                        YeqEmsy_sp    = spMat,
                        BeqEmsy_sp    = spMat,
                        ReqEmsy_sp    = spMat,
                        expBeqEmsy_sp = spMat,
                        Umsy_sp       = spMat,
                        YeqEmsy_p     = rep(0,nP) )

  
  # Calculate ref points
  EmsyRefPts$Emsy_sp    <- Emsy_sp
  EmsyMSRefPts$EmsyMS_p   <- Emsy_p

  # Loop and get each stock/species ref pt
  for( p in 1:nP )
  {
    for( s in 1:nS )
    {
      tmp <- .calcEquil( f = Emsy_sp[s,p], obj = obj, type = "effort" )
      EmsyRefPts$yprEmsy_sp[s,p]    <- tmp$ypr_sp[s,p]
      EmsyRefPts$YeqEmsy_sp[s,p]    <- tmp$Yeq_sp[s,p]
      EmsyRefPts$BeqEmsy_sp[s,p]    <- tmp$Beq_sp[s,p]
      EmsyRefPts$expBeqEmsy_sp[s,p] <- tmp$expBeq_sp[s,p]
      EmsyRefPts$NeqEmsy_sp[s,p]    <- tmp$Neq_sp[s,p]
      EmsyRefPts$ReqEmsy_sp[s,p]    <- tmp$Req_sp[s,p]
      EmsyRefPts$Umsy_sp[s,p]       <- tmp$Ueq_sp[s,p]

      tmpMS <- .calcEquil( f = Emsy_p[p], obj = obj, type = "effort" )
      EmsyMSRefPts$yprEmsy_sp[s,p]    <- tmpMS$ypr_sp[s,p]
      EmsyMSRefPts$YeqEmsy_sp[s,p]    <- tmpMS$Yeq_sp[s,p]
      EmsyMSRefPts$BeqEmsy_sp[s,p]    <- tmpMS$Beq_sp[s,p]
      EmsyMSRefPts$expBeqEmsy_sp[s,p] <- tmpMS$expBeq_sp[s,p]
      EmsyMSRefPts$NeqEmsy_sp[s,p]    <- tmpMS$Neq_sp[s,p]
      EmsyMSRefPts$ReqEmsy_sp[s,p]    <- tmpMS$Req_sp[s,p]
      EmsyMSRefPts$Umsy_sp[s,p]       <- tmpMS$Ueq_sp[s,p]

    }
    # Now sum species eq yield for complex opt yield
    EmsyMSRefPts$YeqEmsy_p[p]         <- sum( EmsyMSRefPts$YeqEmsy_sp[,p] )
  }

  outlist <- list(  EmsyRefPts    = EmsyRefPts,
                    EmsyMSRefPts  = EmsyMSRefPts )
}

# .getFmsy     ()
# Purpose:     fit a spline function to f vs yield, then use a root finder to get Fmsy. 
# Parameters:  obj=list of all operating model parameters, schedules, equilibrium functions
# Returns:     a list with all equilibrium quantities for Fmsy
# Source:      S.P. Cox, modified for hierSCAL by SDNJ
.getFmsy_sp <- function( obj, refCurves )
{
  Fseq  <- refCurves$F
  nS    <- obj$om$nS
  nP    <- obj$om$nP
  nMICE <- obj$om$nMICE

  .getFmsy <- function( yieldCurve, F = Fseq )
  {
    minF <- 0
    maxF <- max(F)

    # Now check if yieldCurve goes negative anywhere
    if( any(yieldCurve < 0) )
    {
      minF <- F[min(which(yieldCurve >= 0))]
      maxF <- F[max(which(yieldCurve >= 0))]
    }

    fySplineFun <- splinefun( x=F, y=yieldCurve )  
    # Find stat point for Fmsy
    Fmsy <- try( uniroot( f = fySplineFun, interval = c(minF, maxF),
                          deriv = 1 )$root )
    
    if(class(Fmsy) == "try-error")
    {
      # browser(cat("uniroot failed\n\n"))
      Fmsy <- 0
    }

    Fmsy <- min(Fmsy, maxF)

    return(Fmsy)
  }

  # calculate Fmsy for each species/stock
  Fmsy_sp <-  apply(  X = refCurves$Yeq_spf, FUN = .getFmsy,
                      MARGIN = c(1,2) )
  # Now create a matrix to hold the species/stock values
  # on each curve
  spMat       <- matrix(0, nrow = nS, ncol = nP)
  dimnames(spMat) <- dimnames(Fmsy_sp)

  FmsyRefPts  <- list(  yprFmsy_sp    = spMat,
                        YeqFmsy_sp    = spMat,
                        EYeqFmsy_sp   = spMat,
                        BeqFmsy_sp    = spMat,
                        ReqFmsy_sp    = spMat,
                        totBeqFmsy_sp = spMat,
                        Umsy_sp       = spMat )

  
  # Calculate ref points
  FmsyRefPts$Fmsy_sp    <- Fmsy_sp
  # # Reset lastB_p
  nothing <- .calcEquil( f = 0, obj = obj )


  # Loop and get each stock/species ref pt
  for( s in 1:nS )
    for( p in 1:nP )
    {
      tmp <- .calcEquil( f = Fmsy_sp[s,p], obj = obj )
      FmsyRefPts$yprFmsy_sp[s,p]    <- tmp$ypr_sp[s,p]
      FmsyRefPts$YeqFmsy_sp[s,p]    <- tmp$Yeq_sp[s,p]
      FmsyRefPts$EYeqFmsy_sp[s,p]   <- tmp$EYeq_sp[s,p]
      FmsyRefPts$BeqFmsy_sp[s,p]    <- tmp$Beq_sp[s,p]
      FmsyRefPts$totBeqFmsy_sp[s,p] <- tmp$totBeq_sp[s,p]
      FmsyRefPts$NeqFmsy_sp[s,p]    <- tmp$Neq_sp[s,p]
      FmsyRefPts$ReqFmsy_sp[s,p]    <- tmp$Req_sp[s,p]
      FmsyRefPts$Umsy_sp[s,p]       <- tmp$Ueq_sp[s,p]
    }

  FmsyRefPts
}     # END function .getFmsy


# .calcGenTime
# Purpose:     Calculate generation time as average age of mature stock
# Parameters:  natural mortality and maturity
# Returns:     generation time
# Source:      S.P. Cox
.calcGenTime <- function( obj )
{
  # Compute eqbm spawning biomass per recruit for
  # given f and species/stock pars
  nS    <- obj$om$nS
  nP    <- obj$om$nP
  A_s   <- obj$om$A_s
  nA    <- obj$om$nA
  nL    <- obj$om$nL
  nX    <- obj$om$nX
  M_xsp <- obj$om$M_xsp

  # Life history schedules
  matAge_asp        <- obj$om$matAge_asp

  genTime_sp <- array(0, dim = c(nS,nP))

  for(s in 1:nS )
    for( p in 1:nP )
    {
      surv <- rep(1, length.out = A_s[s])
      a <- c(1:A_s[s])
      surv[a] <- exp( -M_xsp[nX,s,p] * (a - 1) )
      surv[A_s[s]] <- surv[A_s[s]] / (1 - exp(-M_xsp[nX,s,p]))
      genTime_sp[s,p] <- sum( a * matAge_asp[a,s,p] * surv) / sum( surv * matAge_asp[a,s,p])
    }

  genTime_sp

  
}


## Empirical reference point functions

makeStochEmpRefPts <- function( histFolder = 'fit_SOG_2024_Jan23',
                                groupFolder = "SOG_OM1_refPts",
                                maxXspline = .5,
                                plot = TRUE,
                                yrRange = -c(200:1),
                                qProbs = c(0.25,0.5,0.75),
                                area = "SOG")
{
  # First load history
  histRpt <- readRDS(file.path("./history",histFolder,paste0(histFolder,".rds")))

  # then load empirical ref pts
  load(file.path("./Outputs",groupFolder,"empRefCurves.RData"))
  

  refPtsTable <- makeRefPtsTable(   rpt = histRpt,
                                    EmpRefCurves = saveEmpRefCurves,
                                    maxXspline = maxXspline,
                                    plot = plot,
                                    yrRange = yrRange,
                                    qProbs = qProbs, 
                                    groupFolder = groupFolder, 
                                    area = area)

  refPtsTable$Scenario <- histRpt$ctlList$ctrl$modelHypName

  refPtsTable
}


# makeRefPtsTable()
makeRefPtsTable  <- function( rpt = HGrpt,
                              EmpRefCurves = HGempRefCurves,
                              maxXspline = .5,
                              plot = FALSE,
                              yrRange = -(200:1),
                              qProbs = c(0.25, 0.5, 0.75),
                              groupFolder = groupFolder, 
                              area = area)
{
  # First, estimate the USR

  Dstar   <- rpt$ctlList$hypo$Dstar - 1950
  Dscalar <- rpt$ctlList$hypo$Dscalar
  SB_it   <- rpt$posts$SB_ipt[,1,]
  nMCMC   <- dim(SB_it)[1]

  Busr_i  <- apply(X = Dscalar*SB_it[,Dstar], FUN = mean, MARGIN = 1)

  h_i       <- rpt$posts$h_ip[,1]
  B0_i      <- rpt$posts$B0_ip[,1]
  R0_i      <- rpt$posts$R0_ip[,1]

  # ... update solve sim eqbria produce a distribution.
  refPtsLabs <- c("B0","Busr","Uusr","Yusr","Bmsy","Umsy","MSY","Bcrash","Ucrash","DT")
  eqTable <- matrix(NA, nrow = 3, ncol = length(refPtsLabs))
  colnames(eqTable) <- refPtsLabs
  eqTable <- as.data.frame(eqTable)


  
  # Check if we're using PRD
  stockNames <- rpt$ctlList$ctrl$stockNames

  if(! "PRD" %in% stockNames)
  {
    nReps <- dim(EmpRefCurves$C_isptk)[1]
    # Get posterior samples to match the empirical ref curves
    set.seed(1234)
    
    samples <- sample(1:nMCMC, size = nReps, replace = FALSE )
    Busr_i  <- Busr_i[samples[1:nReps]]

    Busr_q  <- quantile(Busr_i, probs = qProbs )
    SB_qT <- quantile(SB_it[,72], probs = qProbs)

    


    # message("SAR is", rpt$ctlList$data$stock)

    eqListQuant <- solveSimEqbriaQuant( empRefCurves = EmpRefCurves,
                                          maxXspline = maxXspline,
                                          USR_q = Busr_q,
                                          yrRange = yrRange,
                                          qProbs = qProbs )  
    if(plot){
      # pdf("Outputs/referencePtsExplainerHG.pdf", width = 7, height = 10)
      # plotSimEqbriaQuant( empRefCurves = EmpRefCurves,
      #                     maxXspline = maxXspline,
      #                     yrRange = yrRange,
      #                     USR_q = Busr_q, qProbs = qProbs )
      # dev.off()
      pdf(paste0("Outputs/",area,"_A.pdf"), width = 9, height = 7)
      plotRefPtsExpA(empRefCurves = EmpRefCurves,
                    maxXspline = maxXspline,
                    yrRange = yrRange,
                    USR_q = Busr_q, qProbs = qProbs,
                    area = area)
      dev.off()
      pdf(paste0("Outputs/",area,"_B.pdf"), width = 7, height = 9)
      plotRefPtsExpB(empRefCurves = EmpRefCurves,
                     maxXspline = maxXspline,
                     yrRange = yrRange,
                     USR_q = Busr_q, qProbs = qProbs,
                     area = area)
      dev.off()
      pdf(paste0("Outputs/",area,"_C.pdf"), width = 9, height = 7)
      plotRefPtsExpC(empRefCurves = EmpRefCurves,
                     maxXspline = maxXspline,
                     yrRange = yrRange,
                     USR_q = Busr_q, qProbs = qProbs,
                     area = area,
                     history = T, folder = groupFolder)
      dev.off()
    }
    
  }

  if("PRD" %in% stockNames)
  {
    Busr_q  <- quantile(Busr_i, probs = qProbs )
    
    # Build new eqListQuant with ref pts

    eqListQuant <- list()
    eqListQuant$B0_q      <- quantile(rpt$posts$B0_ip[,1], probs = qProbs)
    
    Yeq_if  <- rpt$posts$Yeq_ipf[,1,]
    Ueq_if  <- rpt$posts$Ueq_ipf[,1,]
    Beq_if  <- rpt$posts$SBeq_ipf[,1,]

    Uusr_i  <- Yusr_i <- rep(0,nMCMC)
    Ucra_i  <- Bcra_i <- rep(0,nMCMC)

    # Need to estimate Uusr, Busr, Yusr
    for( i in 1:nMCMC )
    {
      YUspline <- splinefun(  x = Ueq_if[i,], y = Yeq_if[i,] )
      Useq <- seq(from = 0, to = max(Ueq_if[i,]), length.out = 100)
     
      if(Beq_if[i,1] > Busr_i[i])
      {
        BUspline <- splinefun(  x = Ueq_if[i,], y = Beq_if[i,] - Busr_i[i] ) 
        Uusr <- try(uniroot(f = BUspline, interval = c(0,max(Ueq_if[i,])))$root)
        if(class(Uusr) != "try-error")
          Uusr_i[i] <- Uusr
        else browser()

      }
      else Uusr_i[i] <- 0
      
      Yusr_i[i] <- YUspline(Uusr_i[i])

      crashIdx  <- max(which(Beq_if[i,] > 0) )
      Ucra_i[i] <- Ueq_if[i,crashIdx]
      Bcra_i[i] <- Beq_if[i,crashIdx]
    }


    eqListQuant$Uusr_q    <- quantile(Uusr_i, probs = qProbs)
    eqListQuant$YeqUSR_q  <- quantile(Yusr_i, probs = qProbs)

    eqListQuant$Bmsy_q    <- quantile(rpt$posts$Bmsy_ip[,1], probs = qProbs)
    eqListQuant$Umsy_q    <- quantile(rpt$posts$Umsy_ip[,1], probs = qProbs)
    eqListQuant$MSY_q     <- quantile(rpt$posts$MSY_ip[,1], probs = qProbs)

    # Need to estimate Bcrash, Ucrash
    eqListQuant$Ucrash_q  <- quantile(Ucra_i, probs = qProbs)
    eqListQuant$Bcrash_q  <- quantile(Bcra_i, probs = qProbs)

    SB_qT <- quantile(SB_it[,72], probs = qProbs)
  }


 


  

  for( i in 1:3)
  {
   
    eqTable[i,"B0"]     <- eqListQuant$B0_q[i]
    eqTable[i,"Busr"]   <- Busr_q[i]
    eqTable[i,"Uusr"]   <- eqListQuant$Uusr_q[i]
    eqTable[i,"Yusr"]   <- eqListQuant$YeqUSR_q[i]
    eqTable[i,"Bmsy"]   <- eqListQuant$Bmsy_q[i]
    eqTable[i,"Umsy"]   <- eqListQuant$Umsy_q[i]
    eqTable[i,"MSY"]    <- eqListQuant$MSY_q[i]
    eqTable[i,"Bcrash"] <- eqListQuant$Bcrash_q[i]
    eqTable[i,"Ucrash"] <- eqListQuant$Ucrash_q[i]
    eqTable[i,"DT"]     <- SB_qT[i]/eqListQuant$B0_q[i]

  }

  halfIQR <- function(x)
  {
    halfIQR <- 0.5*abs(range(x)[2] - range(x)[1])
    halfIQR
  }

  refPtsTable <- eqTable |>
                  summarise(  mB0      = format(round(median(B0),2), nsmall=2),
                              sdB0     = format(round(halfIQR(B0),2), nsmall=2), 
                              mBusr    = format(round(median(Busr),2), nsmall=2),
                              sdBusr   = format(round(halfIQR(Busr),2), nsmall=2), 
                              mUusr    = format(round(median(Uusr),3), nsmall=3),
                              sdUusr   = format(round(halfIQR(Uusr),3), nsmall=3), 
                              mYusr    = format(round(median(Yusr),3), nsmall=3),
                              sdYusr   = format(round(halfIQR(Yusr),3), nsmall=3),
                              mBmsy    = format(round(median(Bmsy),2), nsmall=2),
                              sdBmsy   = format(round(halfIQR(Bmsy),2), nsmall=2), 
                              mUmsy    = format(round(median(Umsy),3), nsmall=3),
                              sdUmsy   = format(round(halfIQR(Umsy),3), nsmall=3), 
                              mMSY     = format(round(median(MSY),3), nsmall=3),
                              sdMSY    = format(round(halfIQR(MSY),3), nsmall=3),
                              mBcrash  = format(round(median(Bcrash),2), nsmall=2),
                              sdBcrash = format(round(halfIQR(Bcrash),2), nsmall=2),  
                              mUcrash  = format(round(median(Ucrash),3), nsmall=3),
                              sdUcrash = format(round(halfIQR(Ucrash),3), nsmall=3),
                              mDT      = format(round(median(DT),2), nsmall=2),
                              sdDT     = format(round(halfIQR(DT),2), nsmall=2)  ) |>
                  mutate( B0      = paste(mB0, " (", sdB0 ,")",sep = ""),
                          Busr    = paste(mBusr, " (", sdBusr ,")",sep = ""),
                          Uusr    = paste(mUusr, " (", sdUusr ,")",sep = ""),
                          Yusr    = paste(mYusr, " (", sdYusr ,")",sep = ""),
                          Bmsy    = paste(mBmsy, " (", sdBmsy ,")",sep = ""),
                          Umsy    = paste(mUmsy, " (", sdUmsy ,")",sep = ""),
                          MSY     = paste(mMSY, " (", sdMSY ,")",sep = ""),
                          Bcrash  = paste(mBcrash, " (", sdBcrash ,")",sep = ""),
                          Ucrash  = paste(mUcrash, " (", sdUcrash ,")",sep = ""),
                          DT      = paste(mDT, " (", sdDT ,")",sep = "") ) |>
                dplyr::select(  B0,
                                Busr,
                                Uusr,
                                Yusr,
                                Bmsy,
                                Umsy,
                                MSY,
                                Bcrash,
                                Ucrash,
                                DT )


  refPtsTable
}

compareISCAMrefPts <- function( rpt           = HGrpt,
                                empRefCurves  = HGempRefCurves,
                                iscamRpt      = iscamHG,
                                yrRange       = -(150:50) )
{

  # First, estimate the USR

  Dstar   <- rpt$ctlList$hypo$Dstar - 1950
  Dscalar <- rpt$ctlList$hypo$Dscalar
  SB_it   <- rpt$posts$SB_ipt[,1,]

  qProbs <- c(0.25,0.5,0.75)
  SB_qT <- quantile(SB_it[,72], probs = qProbs)

  # Get posterior samples to match the empirical ref curves
  set.seed(1234)
  nMCMC <- dim(SB_it)[1]
  nReps <- dim(empRefCurves$C_isptk)[1]

  samples <- sample(1:nMCMC, size = nReps, replace = FALSE )

  Busr_i  <- apply(X = Dscalar*SB_it[,Dstar], FUN = mean, MARGIN = 1)

  h_i       <- rpt$posts$h_ip[,1]
  B0_i      <- rpt$posts$B0_ip[,1]
  R0_i      <- rpt$posts$R0_ip[,1]



  Busr_i  <- Busr_i[samples[1:nReps]]
  Busr_q  <- quantile(Busr_i, probs = qProbs )

  eqListQuant <- solveSimEqbriaQuant( empRefCurves = empRefCurves,
                                      maxXspline = .5,
                                      USR_q = Busr_q,
                                      yrRange = yrRange,
                                      qProbs = qProbs )  


  iscamMSY  <- round(iscamRpt$msy[1],2)
  iscamBmsy <- round(iscamRpt$bmsy,2)

  iscamUmsy <- round(iscamMSY/(iscamMSY + iscamBmsy),3)

  outTable <- array(0, dim = c(2,4))
  colnames(outTable) <- c("Model","MSY","Bmsy","Umsy")
  outTable <- as.data.frame(outTable)

  outTable[,"Model"] <- c("SISCAH","ISCAM")

  outTable[1,"MSY"]   <- round(eqListQuant$MSY_q[2],2)
  outTable[1,"Bmsy"]  <- round(eqListQuant$Bmsy_q[2],2)
  outTable[1,"Umsy"]  <- round(eqListQuant$Umsy_q[2],3)

  outTable[2,2:4] <- c(iscamMSY,iscamBmsy,iscamUmsy)

  outTable
}


solveSimEqbriaQuant <- function(  empRefCurves,
                                  maxXspline,
                                  sIdx = 1, pIdx = 1,
                                  USR_q = 38, qProbs,
                                  reps = NULL,
                                  yrRange = -(150:50))
{



  # Pull empirical ref curves
  if(is.null(reps))
    reps <- 1:dim(empRefCurves$C_isptk)[1]

  pT <- dim(empRefCurves$C_isptk)[4]
  yrRange <- pT+yrRange
  yrRange <- round(yrRange)

  C_qtk <- apply(X = empRefCurves$C_isptk[reps,sIdx,pIdx,yrRange,,drop = FALSE], FUN = quantile, MARGIN = 4:5, probs = qProbs)
  Y_qtk <- apply(X = empRefCurves$Y_isptk[reps,sIdx,pIdx,yrRange,,drop = FALSE], FUN = quantile, MARGIN = 4:5, probs = qProbs)
  B_qtk <- apply(X = empRefCurves$B_isptk[reps,sIdx,pIdx,yrRange,,drop = FALSE], FUN = quantile, MARGIN = 4:5, probs = qProbs)
  F_qtk <- apply(X = empRefCurves$F_isptk[reps,sIdx,pIdx,yrRange,,drop = FALSE], FUN = quantile, MARGIN = 4:5, probs = qProbs)
  U_qtk <- apply(X = empRefCurves$U_isptk[reps,sIdx,pIdx,yrRange,,drop = FALSE], FUN = quantile, MARGIN = 4:5, probs = qProbs)
  M_qtk <- apply(X = empRefCurves$M_isptk[reps,sIdx,pIdx,yrRange,,drop = FALSE], FUN = quantile, MARGIN = 4:5, probs = qProbs)
  R_qtk <- apply(X = empRefCurves$R_isptk[reps,sIdx,pIdx,yrRange,,drop = FALSE], FUN = quantile, MARGIN = 4:5, probs = qProbs)

  C_qk <- apply(X = C_qtk, FUN = median, MARGIN = c(1,3))
  Y_qk <- apply(X = Y_qtk, FUN = median, MARGIN = c(1,3))
  B_qk <- apply(X = B_qtk, FUN = median, MARGIN = c(1,3))
  F_qk <- apply(X = F_qtk, FUN = median, MARGIN = c(1,3))
  U_qk <- apply(X = U_qtk, FUN = median, MARGIN = c(1,3))
  M_qk <- apply(X = M_qtk, FUN = median, MARGIN = c(1,3))
  R_qk <- apply(X = R_qtk, FUN = median, MARGIN = c(1,3))

  # Order
  Forder <- order(F_qk[which(qProbs == 0.5),])
  C_qk <- C_qk[,Forder]
  Y_qk <- Y_qk[,Forder]
  B_qk <- B_qk[,Forder]
  F_qk <- F_qk[,Forder]
  U_qk <- U_qk[,Forder]
  M_qk <- M_qk[,Forder]
  R_qk <- R_qk[,Forder]
  
  B0_q      <- c(0)
  Bmsy_q    <- c(0)
  Umsy_q    <- c(0)
  Bmsy_q    <- c(0)
  MSY_q     <- c(0)
  Ucrash_q  <- c(0)
  Ccrash_q  <- c(0)
  Bcrash_q  <- c(0)
  Mcrash_q  <- c(0)
  Rcrash_q  <- c(0)
  Uusr_q    <- c(0)
  YeqUSR_q  <- c(0)


  # Make spline
  for( qq in 1:length(qProbs))
  {


    CUspline <- splinefun(x = U_qk[qq,], y = C_qk[qq,])
    YUspline <- splinefun(x = U_qk[qq,], y = Y_qk[qq,])
    BUspline <- splinefun(x = U_qk[qq,], y = B_qk[qq,])
    MUspline <- splinefun(x = U_qk[qq,], y = M_qk[qq,])
    RUspline <- splinefun(x = U_qk[qq,], y = R_qk[qq,])
    
    BUSRspline <- splinefun(x = U_qk[qq,], y = B_qk[qq,] - USR_q[qq])

    
    maxUidx <- length(B_qk[qq,])

    # max(which(B_qk[qq,] > 1e-4))
    # cat("maxU = ", U_k[maxUidx],"\n")
    # if(U_k[maxUidx] < 1e-4)
    #   browser()
    # Note, the curve isn't smooth, so we need
    # to account for local optima when finding
    # these crash points - optimise and uniroot are
    # not super clever.
    minUcrash <- max(0,max(Ucrash_q))
    Useq <- seq(from = 0, to = U_qk[qq,maxUidx], length.out = 1000)
    Yseq <- YUspline(Useq,deriv = 1)

    Ucrash <- try(optimize( interval = c(minUcrash,U_qk[qq,maxUidx]), f=YUspline, deriv=1))

    
    if( class(Ucrash) == "try-error"){
      Ucrash <- list()
      Ucrash$minimum <- 0
    }

    
    k <- 0


    while(any(Yseq < YUspline(Ucrash$minimum,deriv = 1)))
    {
      minUidx <- which.min(Yseq)
      if(minUidx == length(Yseq))
        break
      
      newInt <- c(Useq[minUidx - 3], Useq[minUidx + 3])

      Ucrash <- try(optimize( interval = newInt, f=YUspline, deriv=1))
      k <- k+1

      if(k == 4)
        break
    }
    
    Ucrash <- Ucrash$minimum
    

    if(Ucrash > 0  )
    {
      minUmsy <- max(0.01,max(Umsy_q))

      if(sign(YUspline(minUmsy,deriv = 1)) != sign(YUspline(Ucrash,deriv = 1)))
      {
        Umsy  <- try(uniroot(   interval = c(minUmsy,Ucrash),
                                f = YUspline,
                                deriv = 1 )$root)

      }
      else Umsy <- Ucrash

      # Need to add a check that the estimate of MSY from this is great

      Yseq <- YUspline(Useq)
      k <- 0

      while(any(Yseq > YUspline(Umsy)))
      {
      
        maxUidx <- which.max(Yseq)
        newInt <- c(Useq[maxUidx-5],Useq[maxUidx+5])

        Umsy  <- try(uniroot( interval = newInt,
                              f = YUspline,
                              deriv = 1 )$root)
        k <- k+1

        if(k == 4)
          break
      }

    }
    else Umsy <- 0

    Uusr <- 0
    if(USR_q[qq] < B_qk[qq,1])
      Uusr <- try(uniroot(  interval = c(0.0,Umsy),
                            f = BUSRspline,
                            deriv = 0 )$root)
    
    if( class(Umsy) == "try-error")
      Umsy <- 0
    
    if( class(Uusr) == "try-error")
    {
      browser()
      Uusr <- 0
    }
    
    
    
    B0_q[qq]      <- BUspline(0)
    Bmsy_q[qq]    <- BUspline(Umsy) 
    Umsy_q[qq]    <- Umsy
    MSY_q[qq]     <- YUspline(Umsy) 
    Ucrash_q[qq]  <- Ucrash
    Ccrash_q[qq]  <- YUspline(Ucrash) 
    Bcrash_q[qq]  <- BUspline(Ucrash) 
    Mcrash_q[qq]  <- MUspline(Ucrash) 
    Rcrash_q[qq]  <- RUspline(Ucrash) 
    Uusr_q[qq]    <- Uusr
    YeqUSR_q[qq]  <- YUspline(Uusr)
    
  }
  

  outList <- list(  C_qk = C_qk,
                    Y_qk = Y_qk,
                    U_qk = U_qk,
                    F_qk = F_qk,
                    B_qk = B_qk,
                    M_qk = M_qk,
                    R_qk = R_qk,
                    B0_q      = B0_q,
                    Bmsy_q    = Bmsy_q,
                    Umsy_q    = Umsy_q,
                    MSY_q     = MSY_q,
                    Ucrash_q  = Ucrash_q,
                    Ccrash_q  = Ccrash_q,
                    Bcrash_q  = Bcrash_q,
                    Mcrash_q  = Mcrash_q,
                    Rcrash_q  = Rcrash_q,
                    Uusr_q    = Uusr_q,
                    YeqUSR_q  = YeqUSR_q )


  outList
}





solveSimEqbriaRep <- function(  empRefCurves,
                                maxXspline,
                                sIdx = 1, pIdx = 1,
                                USR = 38,
                                iRep = 1)
{



  # Pull empirical ref curves
  pT <- dim(empRefCurves$C_isptk)[4]
  yrRange <- c(pT/4,3*pT/4)
  yrRange <- round(yrRange)
  C_k <- apply(X = empRefCurves$C_isptk[iRep,sIdx,pIdx,yrRange,], FUN = median, MARGIN = 2)
  Y_k <- apply(X = empRefCurves$Y_isptk[iRep,sIdx,pIdx,yrRange,], FUN = median, MARGIN = 2)
  B_k <- apply(X = empRefCurves$B_isptk[iRep,sIdx,pIdx,yrRange,], FUN = median, MARGIN = 2)
  F_k <- apply(X = empRefCurves$F_isptk[iRep,sIdx,pIdx,yrRange,], FUN = median, MARGIN = 2)
  U_k <- apply(X = empRefCurves$U_isptk[iRep,sIdx,pIdx,yrRange,], FUN = median, MARGIN = 2)
  M_k <- apply(X = empRefCurves$M_isptk[iRep,sIdx,pIdx,yrRange,], FUN = median, MARGIN = 2)
  R_k <- apply(X = empRefCurves$R_isptk[iRep,sIdx,pIdx,yrRange,], FUN = median, MARGIN = 2)
  


  # Order
  Forder <- order(F_k)
  C_k <- C_k[Forder]
  Y_k <- Y_k[Forder]
  B_k <- B_k[Forder]
  F_k <- F_k[Forder]
  U_k <- U_k[Forder]
  M_k <- M_k[Forder]
  R_k <- R_k[Forder]

  # Make spline
  CUspline <- splinefun(x = U_k, y = C_k)
  BUspline <- splinefun(x = U_k, y = B_k)
  MUspline <- splinefun(x = U_k, y = M_k)
  RUspline <- splinefun(x = U_k, y = R_k)
  
  BUSRspline <- splinefun(x = U_k, y = B_k - USR)

  maxU <- max(U_k)
  maxC <- max(C_k)
  maxB <- max(B_k)
  maxM <- max(M_k)
  maxR <- max(R_k)

  maxUidx <- max(which(B_k > 1e-6))
  # cat("maxU = ", U_k[maxUidx],"\n")
  # if(U_k[maxUidx] < 1e-4)
  #   browser()

  Ucrash <- try(optimize( interval = c(0,U_k[maxUidx]), f=CUspline, deriv=1))
  
  if( class(Ucrash) == "try-error"){
    Ucrash <- 0
  }else{
    Ucrash <- Ucrash$minimum
  }

  if(Ucrash > 0  )
  {
    if(sign(CUspline(0.01,deriv = 1)) != sign(CUspline(Ucrash,deriv = 1)))
      Umsy  <- try(uniroot(   interval = c(0.01,Ucrash),
                              f = CUspline,
                              deriv = 1 )$root)
    else Umsy <- Ucrash
  }
  else Umsy <- 0

  Uusr <- 0
  if(USR < B_k[1])
    Uusr <- try(uniroot(  interval = c(0.0,Umsy),
                          f = BUSRspline,
                          deriv = 0 )$root)
  
  if( class(Umsy) == "try-error")
    Umsy <- 0
  
  if( class(Uusr) == "try-error")
  {
    browser()
    Uusr <- 0
  }
  
  
  # Ucrash  <- (U_k[maxUidx])
  
  # U_critical <- try(uniroot(  interval = c(Umsy,U_k[maxUidx]),
  #                       f = CUspline,
  #                       deriv = 2, extendInt = "yes")$root)

  
  
  

  MSY  <- CUspline(Umsy)
  Bmsy <- BUspline(Umsy)

  Ccrash <- CUspline(Ucrash)
  Bcrash <- BUspline(Ucrash)
  Mcrash <- MUspline(Ucrash)
  Rcrash <- RUspline(Ucrash)
  
  YeqUSR  <- CUspline(Uusr)
  

  outList <- list(  C_k = C_k,
                    Y_k = Y_k,
                    U_k = U_k,
                    F_k = F_k,
                    B_k = B_k,
                    M_k = M_k,
                    R_k = R_k,

                    Umsy  = Umsy,
                    Bmsy  = Bmsy,
                    MSY   = MSY ,
                    Ucrash = Ucrash,
                    Ccrash = Ccrash,
                    Bcrash = Bcrash,
                    Mcrash = Mcrash,
                    Rcrash = Rcrash,
                    Uusr    = Uusr,
                    YeqUSR  = YeqUSR  )


  outList
}

plotSimEqbriaQuant <- function( empRefCurves = SOGempRefCurves,
                                maxXspline = 0.18,
                                sIdx = 1, pIdx = 1,
                                USR_q = 38, qProbs,
                                plotRange = FALSE,
                                yrRange = -(450:250)  )
{
  
  eqListQuant <- solveSimEqbriaQuant( empRefCurves = empRefCurves, yrRange = yrRange,
                                      maxXspline = maxXspline, USR_q = USR_q, qProbs = qProbs )



  # ref curves
  C_qk <- eqListQuant$C_qk
  Y_qk <- eqListQuant$Y_qk
  B_qk <- eqListQuant$B_qk
  F_qk <- eqListQuant$F_qk
  U_qk <- eqListQuant$U_qk
  M_qk <- eqListQuant$M_qk
  R_qk <- eqListQuant$R_qk

  maxU <- max(U_qk[2,])
  maxC <- max(C_qk[2,])
  maxB <- max(B_qk[2,])
  maxR <- max(R_qk[2,])
  maxM <- max(M_qk[2,])

  # Ref pts
  Umsy_q  <- eqListQuant$Umsy_q
  Bmsy_q  <- eqListQuant$Bmsy_q
  MSY_q   <- eqListQuant$MSY_q


  Ucrash_q <- eqListQuant$Ucrash_q
  Ccrash_q <- eqListQuant$Ccrash_q
  Bcrash_q <- eqListQuant$Bcrash_q
  Mcrash_q <- eqListQuant$Mcrash_q
  Rcrash_q <- eqListQuant$Rcrash_q
  
##
  Uusr_q    <- eqListQuant$Uusr_q
  YeqUSR_q  <- eqListQuant$YeqUSR_q

##
  empUmsy_q <- round(Umsy_q,2)
  empUusr_q <- round(Uusr_q,3)
  empBmsy_q <- round(Bmsy_q,2)
  empB0_q   <- round(B_qk[,1],2)
  empMSY_q  <- round(MSY_q,2)
  YeqUSR_q  <- round(YeqUSR_q,2)

  maxUidx <- length(U_qk[2,])
  # Ucrash  <- (U_k[maxUidx])
  # Bcrash  <- (B_k[maxUidx])
  # Ccrash  <- (C_k[maxUidx])
  # Mcrash  <- (M_k[maxUidx])
  # Rcrash  <- (R_k[maxUidx])

  # Solve for U_USR


  maxU <- 1.3*max(Ucrash_q)

  par(mfrow = c(4,1),
      mar = c(.1,2,.1,2),
      oma = c(4,4,1,1) )

  qIdx <- which(qProbs %in% range(qProbs))
  medIdx <- which(qProbs == 0.5)

  # First, plot yield
  plot( x = c(0,maxU), y = c(0,1.2*maxC),
        axes = FALSE, type = "n" )
    axis(side = 2, las =1)
    grid()
    box()
    if(plotRange)
      for( qq in qIdx )
      {
        lines( x = U_qk[qq,1:maxUidx], y = Y_qk[qq,1:maxUidx], col = "black", lwd = 3, lty = 2)
        points( x = empUmsy_q[qq], y = empMSY_q[qq], pch = 16, col = "grey40", cex = 2)
        points( x = Ucrash_q[qq], y = Ccrash_q[qq], pch = 21, bg = "red", cex = 2)
      }
    lines( x = U_qk[medIdx,1:maxUidx], y = Y_qk[medIdx,1:maxUidx], col = "black", lwd = 3, lty = 1)
    points( x = empUmsy_q[medIdx], y = empMSY_q[medIdx], pch = 16, col = "grey40", cex = 2)
    points( x = Ucrash_q[medIdx], y = Ccrash_q[medIdx], pch = 21, bg = "red", cex = 2)
    
    mtext(side = 2, text = "Yield (kt)", line = 3)
    segments(x0 = Uusr_q[medIdx], y0 = 0, y1 = YeqUSR_q[medIdx], lty = 2, col = "darkgreen", lwd = 2)
    segments(x0 = 0, x1 = Uusr_q[medIdx], y0 = YeqUSR_q[medIdx], lty = 2, col = "darkgreen", lwd = 2)
    
    legend( "topright", bty = "n", cex = 2,
            legend = c( bquote(U[MSY]*"   = " ~ .(empUmsy_q[medIdx])),
                        paste(" MSY   = ", empMSY_q[medIdx], sep = ""),
                        bquote(U[crash]*" = " ~ .(round(Ucrash_q[medIdx],2))),
                        paste("Yield at USR = ", YeqUSR_q[medIdx], sep = "") ),
            pch = c(NA,16,21,NA),
            lty = c(NA,NA,NA,2),
            lwd = c(NA,NA,NA,2),
            col = c(NA,"grey40","black","darkgreen"),
            pt.bg = c(NA,NA,"red",NA) )

  # Biomass
  plot( x = c(0,maxU), y = c(0,maxB),
        axes = FALSE, type = "n" )
    axis(side = 2, las =1)
    grid()
    box()
    if(plotRange)
      for(qq in qIdx)
      {
        lines( x = U_qk[qq,1:maxUidx], y = B_qk[qq,1:maxUidx], col = "black", lwd = 3, lty = 2)  
        points( x = empUmsy_q[qq], y = empBmsy_q[qq], pch = 16, col = "grey40", cex = 2)
        points( x = Ucrash_q[qq], y = Bcrash_q[qq], pch = 21, bg = "red", cex = 2)  
      }
    lines( x = U_qk[medIdx,1:maxUidx], y = B_qk[medIdx,1:maxUidx], col = "black", lwd = 3, lty = 1)  
    points( x = empUmsy_q[medIdx], y = empBmsy_q[medIdx], pch = 16, col = "grey40", cex = 2)
    points( x = Ucrash_q[medIdx], y = Bcrash_q[medIdx], pch = 21, bg = "red", cex = 2)
    
    segments(x0 = Uusr_q[medIdx], y0 = 0, y1 = USR_q[medIdx], lty = 2, col = "darkgreen", lwd = 2)
    segments(x0 = 0, x1 = Uusr_q[medIdx], y0 = USR_q[medIdx], lty = 2, col = "darkgreen", lwd = 2)
    

    # points( x = U_critical, y = B_critical, pch = 24, bg = "green", cex = 2)
    legend( "topright", bty = "n", cex = 2,
            legend = c( bquote(B[0]*" = " ~ .(empB0_q[2])),
                        bquote(B[MSY]*" = " ~ .(empBmsy_q[2])),
                        bquote(B[USR]*" = " ~ .(round(USR_q[2],2))),
                        bquote(B[crash]*" = " ~ .(round(Bcrash_q[2],2))),
                        bquote(U[USR]*" = " ~ .(empUusr_q[2])) ),
            pch = c(16,NA,21,NA),
            lty = c(NA,2,NA,2),
            lwd = c(NA,2,NA,2),
            col = c("grey40","darkgreen","black","darkgreen"),
            pt.bg = c(NA,NA,"red",NA)

             )
    mtext(side = 2, text = "Spawning\nBiomass (kt)", line = 2.5)

  # Mortality
  plot( x = c(0,maxU), y = c(0,maxM),
        axes = FALSE, type = "n" )
    axis(side = 2, las =1)
    grid()
    box()
    if(plotRange)
      for( qq in qIdx)
      {
        lines( x = U_qk[qq,1:maxUidx], y = M_qk[qq,1:maxUidx], col = "black", lwd = 3, lty = 2)
        points( x = Ucrash_q[qq], y = Mcrash_q[qq], pch = 21, bg = "red", cex = 2)
      }
    lines( x = U_qk[medIdx,1:maxUidx], y = M_qk[medIdx,1:maxUidx], col = "black", lwd = 3, lty = 1)
    points( x = Ucrash_q[medIdx], y = Mcrash_q[medIdx], pch = 21, bg = "red", cex = 2)
    mtext(side = 2, text = "Natural\nMortality (/yr)", line = 2.5)

  # Recruitment
  plot( x = c(0,maxU), y = c(0,maxR),
        axes = FALSE, type = "n" )
    axis(side = 2, las =1)
    axis(side = 1)
    grid()
    box()
    if(plotRange)
      for( qq in qIdx)
      {
        lines( x = U_qk[qq,1:maxUidx], y = R_qk[qq,1:maxUidx], col = "black", lwd = 3, lty = 2)
        points( x = Ucrash_q[qq], y = Rcrash_q[qq], pch = 21, bg = "red", cex = 2)
      }
    lines( x = U_qk[medIdx,1:maxUidx], y = R_qk[medIdx,1:maxUidx], col = "black", lwd = 3, lty = 1)
    points( x = Ucrash_q[medIdx], y = Rcrash_q[medIdx], pch = 21, bg = "red", cex = 2)
    mtext(side = 2, text = "Recruitment (1e6)", line = 3)

    mtext(side = 1, text = "Harvest rate C / (C + B) (/yr)", line = 2.5)

}

plotSimEqbriaRep <- function( saveEmpRefCurves = SOGempRefCurves,
                              maxXspline = 0.18,
                              sIdx = 1, pIdx = 1,
                              USR = 38, iRep = 1 )
{
  
  pT <- dim(saveEmpRefCurves$M_sptk)[3]
  saveEmpRefCurves$M_spk  <- apply(X = saveEmpRefCurves$M_sptk[,,pT-0:49,,drop = FALSE], FUN = median, MARGIN = c(1,2,4), na.rm = T )
  saveEmpRefCurves$R_spk  <- apply(X = saveEmpRefCurves$R_sptk[,,pT-0:49,,drop = FALSE], FUN = median, MARGIN = c(1,2,4), na.rm = T )
  saveEmpRefCurves$Y_spk  <- apply(X = saveEmpRefCurves$Y_sptk[,,pT-0:49,,drop = FALSE], FUN = median, MARGIN = c(1,2,4), na.rm = T )
  



  simEqList <- solveSimEqbriaRep( empRefCurves = saveEmpRefCurves,
                                  maxXspline = maxXspline, USR = USR,
                                  iRep = iRep )



  # ref curves
  C_k <- simEqList$C_k
  Y_k <- simEqList$Y_k
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


  Ucrash <- simEqList$Ucrash
  Ccrash <- simEqList$Ccrash
  Bcrash <- simEqList$Bcrash
  Mcrash <- simEqList$Mcrash
  Rcrash <- simEqList$Rcrash
  
##
  Uusr    <- simEqList$Uusr
  YeqUSR  <- simEqList$YeqUSR

##
  empUmsy <- round(Umsy,2)
  empUusr <- round(Uusr,3)
  empBmsy <- round(Bmsy,2)
  empMSY  <- round(MSY,2)

  maxUidx <- length(U_k)
  # Ucrash  <- (U_k[maxUidx])
  # Bcrash  <- (B_k[maxUidx])
  # Ccrash  <- (C_k[maxUidx])
  # Mcrash  <- (M_k[maxUidx])
  # Rcrash  <- (R_k[maxUidx])

  # Solve for U_USR


  maxU <- 1.3*Ucrash

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
    # points( x = U_critical, y = C_critical, pch = 24, bg = "green", cex = 2)
    # points( x = U_k[1:maxUidx], y = C_k[1:maxUidx], col = "black")
    mtext(side = 2, text = "Yield (kt)", line = 3)
    segments(x0 = Uusr, y0 = 0, y1 = YeqUSR, lty = 2, col = "darkgreen", lwd = 2)
    segments(x0 = 0, x1 = Uusr, y0 = YeqUSR, lty = 2, col = "darkgreen", lwd = 2)
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
    segments(x0 = Uusr, y0 = 0, y1 = USR, lty = 2, col = "darkgreen", lwd = 2)
    segments(x0 = 0, x1 = Uusr, y0 = USR, lty = 2, col = "darkgreen", lwd = 2)
    points( x = Ucrash, y = Bcrash, pch = 21, bg = "red", cex = 2)
    # points( x = U_critical, y = B_critical, pch = 24, bg = "green", cex = 2)
    legend( "topright", bty = "n", cex = 2,
            legend = c( paste("Bmsy = ", empBmsy, sep = ""),
                        paste("Busr = ", round(USR,2), sep = ""),
                        paste("Bcrash = ", round(Bcrash,2), sep = ""),
                        paste("Uusr = ", empUusr, sep = "") ) )
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

# plot2dimEmpYieldCurves()
plotAllocRefCurves <- function( filepath = "data/CC_allocRefCurves.RData" )
{
  
  # Load ref curve list
  load(filepath)


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

  pT      <- dim(C_sptk)[3]
  nS      <- dim(C_sptk)[1]
  nP      <- dim(C_sptk)[2]
  
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
        mar = c(.25,1.5,.25,1.5),
        oma = c(5,5,3,3) )


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
    

      X <- U
      xLab  <- "Whole fish harvest rate"
      
      CXspline <- splinefun(x = X, y = Y)
      BXspline <- splinefun(x = X, y = B)

      maxXidx <- max(which(Y > 1e-2))
      maxYidx <- which.max(Y)

      empXmsy  <- try(uniroot(  interval = c(0,X[min(maxYidx + 1, maxXidx)]),
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

      legend( "topleft", bty = "n",cex = 2,
              legend = c( paste( "     ",100*(1 - thisAlloc), "% SOK")) )      

      legend( "topright", bty = "n",cex = 2,
              legend = c( bquote("      "*U[MSY]*" =   " ~ .(format(empXmsy,nsmall=2))),
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

} # END plotEmpYieldCurves()


plotSimEqbriaQuant <- function( empRefCurves = SOGempRefCurves,
                                maxXspline = 0.18,
                                sIdx = 1, pIdx = 1,
                                USR_q = 38, qProbs,
                                plotRange = FALSE,
                                yrRange = -(450:250)  )
{
  
  eqListQuant <- solveSimEqbriaQuant( empRefCurves = empRefCurves, yrRange = yrRange,
                                      maxXspline = maxXspline, USR_q = USR_q, qProbs = qProbs )
  
  
  
  # ref curves
  C_qk <- eqListQuant$C_qk
  Y_qk <- eqListQuant$Y_qk
  B_qk <- eqListQuant$B_qk
  F_qk <- eqListQuant$F_qk
  U_qk <- eqListQuant$U_qk
  M_qk <- eqListQuant$M_qk
  R_qk <- eqListQuant$R_qk
  
  maxU <- max(U_qk[2,])
  maxC <- max(C_qk[2,])
  maxB <- max(B_qk[2,])
  maxR <- max(R_qk[2,])
  maxM <- max(M_qk[2,])
  
  # Ref pts
  Umsy_q  <- eqListQuant$Umsy_q
  Bmsy_q  <- eqListQuant$Bmsy_q
  MSY_q   <- eqListQuant$MSY_q
  
  
  Ucrash_q <- eqListQuant$Ucrash_q
  Ccrash_q <- eqListQuant$Ccrash_q
  Bcrash_q <- eqListQuant$Bcrash_q
  Mcrash_q <- eqListQuant$Mcrash_q
  Rcrash_q <- eqListQuant$Rcrash_q
  
  ##
  Uusr_q    <- eqListQuant$Uusr_q
  YeqUSR_q  <- eqListQuant$YeqUSR_q
  
  ##
  empUmsy_q <- round(Umsy_q,2)
  empUusr_q <- round(Uusr_q,3)
  empBmsy_q <- round(Bmsy_q,2)
  empB0_q   <- round(B_qk[,1],2)
  empMSY_q  <- round(MSY_q,2)
  YeqUSR_q  <- round(YeqUSR_q,2)
  
  maxUidx <- length(U_qk[2,])
  # Ucrash  <- (U_k[maxUidx])
  # Bcrash  <- (B_k[maxUidx])
  # Ccrash  <- (C_k[maxUidx])
  # Mcrash  <- (M_k[maxUidx])
  # Rcrash  <- (R_k[maxUidx])
  
  # Solve for U_USR
  
  
  maxU <- 1.3*max(Ucrash_q)
  
  par(mfrow = c(4,1),
      mar = c(.1,2,.1,2),
      oma = c(4,4,1,1) )
  
  qIdx <- which(qProbs %in% range(qProbs))
  medIdx <- which(qProbs == 0.5)
  
  # First, plot yield
  plot( x = c(0,maxU), y = c(0,1.2*maxC),
        axes = FALSE, type = "n" )
  axis(side = 2, las =1)
  grid()
  box()
  if(plotRange)
    for( qq in qIdx )
    {
      lines( x = U_qk[qq,1:maxUidx], y = Y_qk[qq,1:maxUidx], col = "black", lwd = 3, lty = 2)
      points( x = empUmsy_q[qq], y = empMSY_q[qq], pch = 16, col = "grey40", cex = 2)
      points( x = Ucrash_q[qq], y = Ccrash_q[qq], pch = 21, bg = "red", cex = 2)
    }
  lines( x = U_qk[medIdx,1:maxUidx], y = Y_qk[medIdx,1:maxUidx], col = "black", lwd = 3, lty = 1)
  points( x = empUmsy_q[medIdx], y = empMSY_q[medIdx], pch = 16, col = "grey40", cex = 2)
  points( x = Ucrash_q[medIdx], y = Ccrash_q[medIdx], pch = 21, bg = "red", cex = 2)
  
  mtext(side = 2, text = "Yield (kt)", line = 3)
  segments(x0 = Uusr_q[medIdx], y0 = 0, y1 = YeqUSR_q[medIdx], lty = 2, col = "darkgreen", lwd = 2)
  segments(x0 = 0, x1 = Uusr_q[medIdx], y0 = YeqUSR_q[medIdx], lty = 2, col = "darkgreen", lwd = 2)
  
  legend( "topright", bty = "n", cex = 2,
          legend = c( bquote(U[MSY]*"   = " ~ .(empUmsy_q[medIdx])),
                      paste(" MSY   = ", empMSY_q[medIdx], sep = ""),
                      bquote(U[crash]*" = " ~ .(round(Ucrash_q[medIdx],2))),
                      paste("Yield at USR = ", YeqUSR_q[medIdx], sep = "") ),
          pch = c(NA,16,21,NA),
          lty = c(NA,NA,NA,2),
          lwd = c(NA,NA,NA,2),
          col = c(NA,"grey40","black","darkgreen"),
          pt.bg = c(NA,NA,"red",NA) )
  
  # Biomass
  plot( x = c(0,maxU), y = c(0,maxB),
        axes = FALSE, type = "n" )
  axis(side = 2, las =1)
  grid()
  box()
  if(plotRange)
    for(qq in qIdx)
    {
      lines( x = U_qk[qq,1:maxUidx], y = B_qk[qq,1:maxUidx], col = "black", lwd = 3, lty = 2)  
      points( x = empUmsy_q[qq], y = empBmsy_q[qq], pch = 16, col = "grey40", cex = 2)
      points( x = Ucrash_q[qq], y = Bcrash_q[qq], pch = 21, bg = "red", cex = 2)  
    }
  lines( x = U_qk[medIdx,1:maxUidx], y = B_qk[medIdx,1:maxUidx], col = "black", lwd = 3, lty = 1)  
  points( x = empUmsy_q[medIdx], y = empBmsy_q[medIdx], pch = 16, col = "grey40", cex = 2)
  points( x = Ucrash_q[medIdx], y = Bcrash_q[medIdx], pch = 21, bg = "red", cex = 2)
  
  segments(x0 = Uusr_q[medIdx], y0 = 0, y1 = USR_q[medIdx], lty = 2, col = "darkgreen", lwd = 2)
  segments(x0 = 0, x1 = Uusr_q[medIdx], y0 = USR_q[medIdx], lty = 2, col = "darkgreen", lwd = 2)
  
  
  # points( x = U_critical, y = B_critical, pch = 24, bg = "green", cex = 2)
  legend( "topright", bty = "n", cex = 2,
          legend = c( bquote(B[0]*" = " ~ .(empB0_q[2])),
                      bquote(B[MSY]*" = " ~ .(empBmsy_q[2])),
                      bquote(B[USR]*" = " ~ .(round(USR_q[2],2))),
                      bquote(B[crash]*" = " ~ .(round(Bcrash_q[2],2))),
                      bquote(U[USR]*" = " ~ .(empUusr_q[2])) ),
          pch = c(16,NA,21,NA),
          lty = c(NA,2,NA,2),
          lwd = c(NA,2,NA,2),
          col = c("grey40","darkgreen","black","darkgreen"),
          pt.bg = c(NA,NA,"red",NA)
          
  )
  mtext(side = 2, text = "Spawning\nBiomass (kt)", line = 2.5)
  
  # Mortality
  plot( x = c(0,maxU), y = c(0,maxM),
        axes = FALSE, type = "n" )
  axis(side = 2, las =1)
  grid()
  box()
  if(plotRange)
    for( qq in qIdx)
    {
      lines( x = U_qk[qq,1:maxUidx], y = M_qk[qq,1:maxUidx], col = "black", lwd = 3, lty = 2)
      points( x = Ucrash_q[qq], y = Mcrash_q[qq], pch = 21, bg = "red", cex = 2)
    }
  lines( x = U_qk[medIdx,1:maxUidx], y = M_qk[medIdx,1:maxUidx], col = "black", lwd = 3, lty = 1)
  points( x = Ucrash_q[medIdx], y = Mcrash_q[medIdx], pch = 21, bg = "red", cex = 2)
  mtext(side = 2, text = "Natural\nMortality (/yr)", line = 2.5)
  
  # Recruitment
  plot( x = c(0,maxU), y = c(0,maxR),
        axes = FALSE, type = "n" )
  axis(side = 2, las =1)
  axis(side = 1)
  grid()
  box()
  if(plotRange)
    for( qq in qIdx)
    {
      lines( x = U_qk[qq,1:maxUidx], y = R_qk[qq,1:maxUidx], col = "black", lwd = 3, lty = 2)
      points( x = Ucrash_q[qq], y = Rcrash_q[qq], pch = 21, bg = "red", cex = 2)
    }
  lines( x = U_qk[medIdx,1:maxUidx], y = R_qk[medIdx,1:maxUidx], col = "black", lwd = 3, lty = 1)
  points( x = Ucrash_q[medIdx], y = Rcrash_q[medIdx], pch = 21, bg = "red", cex = 2)
  mtext(side = 2, text = "Recruitment (1e6)", line = 3)
  
  mtext(side = 1, text = "Harvest rate C / (C + B) (/yr)", line = 2.5)
  
}

# plotRefPtsExpA()
# Plots a multi-panel plot
# Reference points explainer for SOG herring
plotRefPtsExpA <- function( empRefCurves = SOGempRefCurves,
                           maxXspline = 0.18,
                           yrRange = yrRange,
                           sIdx = 1, pIdx = 1,
                           USR_q, qProbs,
                           projY = 200,
                           area = 'SOG')
{
  
  eqListQuant <- solveSimEqbriaQuant( empRefCurves = empRefCurves, yrRange = yrRange,
                                      maxXspline = maxXspline, USR_q = USR_q, qProbs = qProbs )
  
  # ref curves
  C_qk <- eqListQuant$C_qk
  Y_qk <- eqListQuant$Y_qk
  B_qk <- eqListQuant$B_qk
  F_qk <- eqListQuant$F_qk
  U_qk <- eqListQuant$U_qk
  
  maxU <- max(U_qk[2,])
  maxC <- max(C_qk[2,])
  maxB <- max(B_qk[2,])
  
  # Ref pts
  Umsy_q  <- eqListQuant$Umsy_q
  Bmsy_q  <- eqListQuant$Bmsy_q
  MSY_q   <- eqListQuant$MSY_q
  
  
  Ucrash_q <- eqListQuant$Ucrash_q
  Ccrash_q <- eqListQuant$Ccrash_q
  Bcrash_q <- eqListQuant$Bcrash_q
  
  maxUidx <- length(U_qk[2,])
  
  qIdx <- which(qProbs %in% range(qProbs))
  medIdx <- which(qProbs == 0.5)
  
  B_itk <- empRefCurves$B_isptk[,1,1,,]
  
  if (area == 'SOG'){
    sim <- c(1,23,47,6)
    orderedsim <- c(1,3,6,14)
    U <- c(0,0.7,0.14,0.26)
    YLIMs <- c(200,160,110,80)
    labels <- c(" HR = 0",expression(' HR = 0.5U'['MSY | OM1']*' = 7%'),
                expression(' HR = U'['MSY | OM1']*' = 14%'), bquote(atop('HR = '*U['crash | OM1'],
                                                                         "= 26%")))
  }
  else if (area == 'HG'){
    sim <- c(1,24,57,9)
    orderedsim <- c(1,3,6,16)
    U <- c(0,0.07,0.14,0.31)
    YLIMs <- c(20,14,12,9)
    labels <- c(" HR = 0",expression(' HR = 0.5U'['MSY | OM1']*' = 7%'),
                expression(' HR = U'['MSY | OM1']*' = 14%'), bquote(atop('HR = '*U['crash | OM1'],
                                                                         "= 31%")))
  }
  
  
  label <- c("HR = 0",expression('HR = 0.5U'['MSY | OM1']),
             expression('HR = U'['MSY | OM1']),expression('HR = U'['crash | OM1']))
  tPlotIdx <- 1:projY
  
  plotLayout <- matrix( c(1,1,5,5,5,5,
                          2,2,5,5,5,5,
                          3,3,5,5,5,5,
                          4,4,5,5,5,5),
                        ncol = 6, byrow = T)
  layout(plotLayout)
  
  par(mar = c(.25,.25,.25,.25),
      oma = c(4,5,1,4), xaxs = 'r', yaxs = 'r', las = 1)
  
  yrs <- 2024:2523
  clrs <- paletteer_d("RColorBrewer::Dark2",4)
  
  for( i in 1:4 )
  {
    # get B for 4 scenarios and plot them
    B_it <- B_itk[,,sim[i]]
    
    B_qt <- apply( X = B_it, FUN = quantile,
                   MARGIN = c(2), probs = c(0.025, 0.5, 0.975),
                   na.rm = T )
    if(i == 4){
      plot( x = range(yrs[tPlotIdx]),
            y = c(-1,YLIMs[i]),
            type = "n", axes = F)
      box()
      grid()
      abline(h = 0,lty = 2, lwd = 2)
    }
    else{
      plot( x = range(yrs[tPlotIdx]),
            y = c(0,YLIMs[i]),
            type = "n", axes = F)
      box()
      grid()
      abline(h = 0,lty = 2, lwd = 2)
    }
    
    
    mfg <- par("mfg")
    if( mfg[1] == mfg[3] )
      axis( side = 1, cex.axis = 1.3 )
    axis( side = 2, las = 1, cex.axis = 1.3)
    # polygon(  x = c(yrs[tPlotIdx], rev(yrs[tPlotIdx])),
    #           y = c(B_qt[1,tPlotIdx], rev(B_qt[3,tPlotIdx])),
    #           col = "grey65", border = NA )
    lines( x = yrs[tPlotIdx], y = B_qt[2,tPlotIdx], lwd = 3, col = clrs[i])
    legend("topright", bty = 'n', legend = label[i],cex = 1.5)
    
    
    if( i==4 ){
      # legend( x = "topright", bty = "n",
      #         legend = c( "Median Spawning Biomass", 
      #                     "Central 95%",
      #                     "Replicate Traces"),
      #         # "Unfished Biomass",
      #         # expression(paste("LRP = 0.3",B[0]))),
      #         # "Median Projection SB" ),
      #         # expression(B[MSY,MS])),
      #         col = c(  "black", "grey65", "black",
      #                   "grey50", "red","salmon" ),
      #         pch = c(NA,15, NA ),# NA, NA,NA,NA),
      #         lty = c(1, NA, 1), # 3, 2, 2),
      #         lwd = c(3, NA, .8) )#, 2, 3, 3 ) )
      mtext( side = 1, outer = F, text = "Year",
             line = 2.5, font = 2)
    }
  }
  mtext( side = 2, outer = TRUE, text = "Spawning Biomass (kt)",
         line = 3, font = 2, las = 0)
  
  # plot yield
  plot( x = c(0,0.4), y = c(0,maxB),
        axes = FALSE, type = "n" )
  axis(side = 4, las =1, cex.axis = 1.3)
  axis(side = 1, cex.axis = 1.3)
  grid()
  box()
  lines( x = U_qk[medIdx,1:maxUidx], y = Y_qk[medIdx,1:maxUidx], col = "black", lwd = 3, lty = 1)
  points(x = U_qk[medIdx,orderedsim[1]], y = Y_qk[medIdx,orderedsim[1]], pch = 16, col = clrs[1], cex = 2)
  points(x = U_qk[medIdx,orderedsim[2]], y = Y_qk[medIdx,orderedsim[2]], pch = 16, col = clrs[2], cex = 2)
  points(x = U_qk[medIdx,orderedsim[3]], y = Y_qk[medIdx,orderedsim[3]], pch = 16, col = clrs[3], cex = 2)
  points(x = Ucrash_q[medIdx], y = Ccrash_q[medIdx], pch = 16, col = clrs[4], cex = 2)
  # points(x = U_qk[medIdx,13], y = Y_qk[medIdx,13], pch = 16, col = clrs[4], cex = 2)
  
  # legend( "topright", bty = "n", cex = 1.5,
  #         legend = c("U = 0","1/2 Umsy = 0.1","Umsy = 0.2","U = 0.25"),
  #         pch = c(16,16,16,16),
  #         col = clrs)
  
  # Biomass
  lines( x = U_qk[medIdx,1:maxUidx], y = B_qk[medIdx,1:maxUidx], col = "black", lwd = 3, lty = 1) 
  
  # with U0
  segments( x0 = U_qk[medIdx,orderedsim[1]], y0 = 0, y1 = B_qk[medIdx,orderedsim[1]], lty = 2, col = clrs[1], lwd = 2)
  points(x = U_qk[medIdx,orderedsim[1]], y = B_qk[medIdx,orderedsim[1]], pch = 16, col = clrs[1], cex = 2)
  text(x = U_qk[medIdx,orderedsim[1]], y = B_qk[medIdx,orderedsim[1]], labels[1], pos = 4, col = clrs[1], cex = 1.5)
  
  # with U = 0.09
  segments(x0 = U_qk[medIdx,orderedsim[2]], y0 = 0, y1 = B_qk[medIdx,orderedsim[2]], lty = 2, col = clrs[2], lwd = 2)
  points(x = U_qk[medIdx,orderedsim[2]], y = B_qk[medIdx,orderedsim[2]], pch = 16, col = clrs[2], cex = 2)
  text(x = U_qk[medIdx,orderedsim[2]], y = B_qk[medIdx,orderedsim[2]], labels[2], pos = 4, col = clrs[2], cex = 1.5)
  
  # with U = 0.2, Umsy
  segments( x0 = U_qk[medIdx,orderedsim[3]], y0 = 0, y1 = B_qk[medIdx,orderedsim[3]], lty = 2, col = clrs[3], lwd = 2)
  points(x = U_qk[medIdx,orderedsim[3]], y = B_qk[medIdx,orderedsim[3]], pch = 16, col = clrs[3], cex = 2)
  text(x = U_qk[medIdx,orderedsim[3]], y = B_qk[medIdx,orderedsim[3]], labels[3], pos = 4, col = clrs[3], cex = 1.5)
  
  # with U = Ucrash
  segments( x0 = Ucrash_q[medIdx], y0 = 0, y1 = Bcrash_q[medIdx], lty = 2, col = clrs[4], lwd = 2)
  points(x = Ucrash_q[medIdx], y = Bcrash_q[medIdx], pch = 16, col = clrs[4], cex = 2)
  text(x = Ucrash_q[medIdx], y = Bcrash_q[medIdx], labels[4], pos = 4, col = clrs[4], cex = 1.5, xpd = NA)
  
  # with U = 0.25
  # segments( x0 = U_qk[medIdx,13], y0 = 0, y1 = B_qk[medIdx,13], lty = 2, col = clrs[4], lwd = 2)
  # points(x = U_qk[medIdx,13], y = B_qk[medIdx,13], pch = 16, col = clrs[4], cex = 2)
  # text(x = U_qk[medIdx,13], y = B_qk[medIdx,13], "25.3%", pos = 4, col = clrs[4], cex = 2)
  
  # label the lines
  if (area == 'SOG'){
    text(x = 0.052, y = 80.5, "Spawning Biomass", pos = 3, srt = -56, cex = 1.5)
    text(x = 0.05, y = 4, "Yield", pos = 3, srt = 19, cex = 1.5)
  }
  else if (area == 'HG'){
    text(x = 0.04, y = 5.5, "Spawning Biomass", pos = 3, srt = -74, cex = 1.5)
    text(x = 0.1, y = 0.2, "Yield", pos = 3, srt = 0, cex = 1.5)
  }
  
  mtext(side = 4, text = "Thousands of tonnes", line = 2.5, font = 2, las = 3)
  mtext(side = 1, text = "Harvest rate", line = 2.5, font = 2)
  
  
}# plotRefPtsExpA()

draw_segments <- function(x, y, threshold, color) {
  # Find indices where y crosses the threshold
  below_threshold <- y < threshold
  
  # Less of the line
  for (i in 1:(length(x) - 1)) {
    if (below_threshold[i] && below_threshold[i + 1]) {
      lines(x[i:(i + 1)], y[i:(i + 1)], col = "red")
    } else {
      lines(x[i:(i + 1)], y[i:(i + 1)], col = color)
    }
  }
  
  # more of the line
  # for (i in 1:(length(x) - 1)) {
  #   # Determine if the current segment crosses the threshold
  #   if ((below_threshold[i] && !below_threshold[i + 1]) || (!below_threshold[i] && below_threshold[i + 1])) {
  #     # Draw a segment with the threshold color
  #     lines(x[i:(i + 1)], y[i:(i + 1)], col = "red")
  #   } else if (below_threshold[i] && below_threshold[i + 1]) {
  #     # Both points are below the threshold
  #     lines(x[i:(i + 1)], y[i:(i + 1)], col = "red")
  #   } else {
  #     # Both points are above the threshold
  #     lines(x[i:(i + 1)], y[i:(i + 1)], col = color)
  #   }
  
}

# plotRefPtsExpB()
plotRefPtsExpB <- function( empRefCurves = saveEmpRefCurves,
                           maxXspline = 0.18,
                           yrRange = yrRange,
                           sIdx = 1, pIdx = 1,
                           USR_q, qProbs,
                           trace = 45,
                           projY = 15,
                           area = 'SOG')
{
  
  eqListQuant <- solveSimEqbriaQuant( empRefCurves = empRefCurves, yrRange = yrRange,
                                      maxXspline = maxXspline, USR_q = USR_q, qProbs = qProbs )
  
  if (area == 'SOG'){
    yr <- 73
    sim <- c(1,23,47) # as ordered in empRefCurves
    orderedsim <- c(1,3,6) # as ordered in sims and eqListQuant
    U <- c(0,0.07,0.14)
    YLIM <- c(0,220)
    label <- c("HR = 0",expression('HR = 0.5U'['MSY | OM1']*' = 7%'),
               expression('HR = U'['MSY | OM1']*' = 14%'))
    # trace <- 45
  }
  else if (area == 'HG'){
    yr <- 72
    sim <- c(1,24,57)
    orderedsim <- c(1,3,6)
    U <- c(0,0.07,0.14)
    YLIM <- c(-10,120)
    label <- c("HR = 0",expression('HR = 0.5U'['MSY | OM1']*' = 7%'),
               expression('HR = U'['MSY | OM1']*' = 14%'))
    # trace <- 105
  }
  
  tPlotIdx <- 1:(yr+projY)
  
  B_qk <- eqListQuant$B_qk
  Bmsy_q  <- eqListQuant$Bmsy_q
  Uusr_q    <- eqListQuant$Uusr_q
  YeqUSR_q  <- eqListQuant$YeqUSR_q
  
  B_itk <- array(0,dim = c(200,323,length(sim)))
  
  for (i in 1:length(sim)){
    if (area == 'SOG'){
      load(paste0("Outputs/SOG_OM1_refPts/sim_parBatSOG_OM1_refPts",orderedsim[i],"/sim_parBatSOG_OM1_refPts",orderedsim[i],".RData"))
    }
    else if (area == 'HG'){
      load(paste0("Outputs/stochRefPtsLong_HG/sim_parBatHG_StochRefPts_",orderedsim[i],"/sim_parBatHG_StochRefPts_",orderedsim[i],".RData"))
    }
    B_itk[,,i] <- blob$om$SB_ispt[,1,1,]
  }
  B_qtk <- apply( X = B_itk, FUN = quantile,
                  MARGIN = c(2,3), probs = c(0.25, 0.5, 0.75),
                  na.rm = T )[,1:88,]
  B_qtk95 <- apply( X = B_itk, FUN = quantile,
                  MARGIN = c(2,3), probs = c(0.025, 0.5, 0.975),
                  na.rm = T )[,1:88,]
  
  # prepare box and whickers
  # C <- list(B_qtk[1,(yr+1):88,1],B_qtk[1,(yr+1):88,2],B_qtk[1,(yr+1):88,3])
  C <- list(B_itk[,(yr+1):88,1],B_itk[,(yr+1):88,2],B_itk[,(yr+1):88,3])
  names(C) <- c("HR = 0","0.5 Umsy","Umsy")
  # Calculate percentiles
  calculate_percentiles <- function(x) {
    quants <- quantile(x, probs = c(0.025, 0.25, 0.75, 0.975))
    return(c(quants[2], quants[3], median(x), quants[1], quants[4]))
  }
  
  # Create a boxplot object
  boxplot_data <- lapply(C, calculate_percentiles)
  names(boxplot_data) <- names(C)
  

  clrs <- paletteer_d("RColorBrewer::Dark2",length(sim)+5)
  clrs <- c(clrs[1],clrs[3],clrs[6])
  
  yrs <- 1951:2523
  
  plotLayout <- matrix( c(1,1,1,1,1,1,1,1,2,
                          1,1,1,1,1,1,1,1,2,
                          3,3,3,3,3,3,3,3,4,
                          3,3,3,3,3,3,3,3,4,
                          5,5,5,5,5,5,5,5,6,
                          5,5,5,5,5,5,5,5,6),
                        ncol = 9, byrow = T)
  layout(plotLayout)
  
  par(mar = c(.25,.25,.25,.25), oma = c(3,4,1,1), 
      xaxs = 'i', yaxs = 'i', las = 1, tck=-.015,
      mgp=c(2, 0.4,0))
  
  for( i in 1:length(sim)){
    plot( x = range(yrs[tPlotIdx]),
          y = c(0,190),
          xlim = c(1951,2038),
          ylim = YLIM,
          type = "n", axes = F)
    grid(ny = NA)
    axis(side = 2)
    if(i == 3){
      axis(side = 1)
    }
    
    # plot history
    # polygon(x = c(yrs[1:yr], rev(yrs[1:yr])),
    #         y = c(B_qtk95[1,1:yr,1], rev(B_qtk95[3,1:yr,1])),
    #         col = "grey80", border = NA )
    polygon(x = c(yrs[1:yr], rev(yrs[1:yr])),
            y = c(B_qtk[1,1:yr,1], rev(B_qtk[3,1:yr,1])),
            col = "grey50", border = NA )
    # lines(x = yrs[1:yr], y = B_qtk[1,1:yr,1], col = "black", lwd = 2.5)
    
    # plot projection
    # polygon(x = c(yrs[yr:length(tPlotIdx)], rev(yrs[yr:length(tPlotIdx)])),
    #         y = c(B_qtk95[1,yr:length(tPlotIdx),i], rev(B_qtk95[3,yr:length(tPlotIdx),i])),
    #         col = adjustcolor(clrs[i],alpha.f = 0.2), border = NA )
    polygon(x = c(yrs[yr:length(tPlotIdx)], rev(yrs[yr:length(tPlotIdx)])),
            y = c(B_qtk[1,yr:length(tPlotIdx),i], rev(B_qtk[3,yr:length(tPlotIdx),i])),
            col = adjustcolor(clrs[i], alpha.f = 0.5), border = NA )
    
    
    points(x = 1951, y = Bmsy_q[2], col = 'black', pch = 16, xpd = NA, cex = 1.5)
    text(x = 1951, y = Bmsy_q[2],expression('B'['MSY']), pos = 4, cex = 1.5)
    abline(h = 0.3*B_qk[2,1], lty = 2)
    
    # historical median
    # lines(x = yrs[1:yr], y = B_qtk95[2,1:yr,1], col = "black", lwd = 3.5, lty = 1)
    # projection median
    # lines(x = yrs[yr:length(tPlotIdx)], y = B_qtk[2,yr:length(tPlotIdx),i], col = clrs[i], lwd = 3.5)
    
    # 1 trace
    # SOG: 45, HG: 
    # lines(x = yrs[1:length(tPlotIdx)], y = B_itk[trace,1:length(tPlotIdx),i], col = clrs[i], lwd = 2)
    draw_segments(x = yrs[1:length(tPlotIdx)], y = B_itk[trace,1:length(tPlotIdx),i], 
                  threshold = 0.3*B_qk[2,1], color = clrs[i])
    
    legend('topleft',label[i], bty = 'n', cex = 1.5)
    box()
    
    if (area == 'SOG'){
      text(x = 1995, y = 0.3*B_qk[2,1], "30% B0", pos = 3, cex = 1.5)
      text(x = 1995, y = 0.3*B_qk[2,1], "Limit Reference Point", pos = 1, cex = 1.5)
    }
    else if (area == 'HG'){
      text(x = 1980, y = 0.3*B_qk[2,1], "30% B0", pos = 3, cex = 1.5)
      text(x = 1980, y = 0.3*B_qk[2,1], "Limit Reference Point", pos = 1, cex = 1.5)
    }
    
    # points and horizontal abline
    # avrg <- apply(X = B_qtk[1,(yr+1):88,], FUN = mean, MARGIN = 2)
    # points(x = rep(2038,length(sim)), y = avrg, col = clrs, pch = 16, xpd = NA, cex = 1.5)
    # for (i in 1:length(U)) {
      # segments(x0 = 2024, y0 = avrg[i], x1 = 2038, y1 = avrg[i], lty = 2, col = clrs[i], lwd = 1)
      # mtext(U[i], side=4, at=avrg[i], col=clrs[i], cex=1, las=1, line = 1)
    # }
    
    # box and whisker
    tmp <- boxplot_data[i]
  
    # Convert to a proper boxplot object
    boxplot_stats <- boxplot(C[i], plot = FALSE)
    boxplot_stats$stats[1,1] <- tmp[[1]][4] # 2.5th percentile
    boxplot_stats$stats[2,1] <- tmp[[1]][1] # 25th percentile
    boxplot_stats$stats[3,1] <- tmp[[1]][3] # median
    boxplot_stats$stats[4,1] <- tmp[[1]][2] # 75th percentile
    boxplot_stats$stats[5,1] <- tmp[[1]][5] # 97.5th percentile
    boxplot_stats$out <- numeric(0) # Remove outliers
    
    # Plot the boxplot
    bxp(boxplot_stats, border = clrs[i], col = clrs[i], axes = FALSE, ylim = YLIM)
    
  }
  mtext("Spawning Biomass (kt)", side = 2, line = 2, las = 0, outer = T)
  
}

# yield curve with SB on X axis and yield on Y, with 20 years of historic data
plotRefPtsExpC <- function( empRefCurves = SOGempRefCurves,
                            maxXspline = 0.18,
                            yrRange = yrRange,
                            sIdx = 1, pIdx = 1,
                            USR_q, qProbs,
                            area = 'SOG',
                            history = T,
                            histYrs = 20,
                            folder = 'SOG_OM1_refPts')
{
  
  eqListQuant <- solveSimEqbriaQuant( empRefCurves = empRefCurves, yrRange = yrRange,
                                      maxXspline = maxXspline, USR_q = USR_q, qProbs = qProbs )
  
  # ref curves
  C_qk <- eqListQuant$C_qk
  Y_qk <- eqListQuant$Y_qk
  B_qk <- eqListQuant$B_qk
  F_qk <- eqListQuant$F_qk
  U_qk <- eqListQuant$U_qk
  
  maxU <- max(U_qk[2,])
  maxC <- max(C_qk[2,])
  maxB <- max(B_qk[2,])
  
  # Ref pts
  Umsy_q  <- eqListQuant$Umsy_q
  Bmsy_q  <- eqListQuant$Bmsy_q
  MSY_q   <- eqListQuant$MSY_q
  
  
  Ucrash_q <- eqListQuant$Ucrash_q
  Ccrash_q <- eqListQuant$Ccrash_q
  Bcrash_q <- eqListQuant$Bcrash_q
  
  ##
  empUmsy_q <- round(Umsy_q,2)
  empBmsy_q <- round(Bmsy_q,2)
  empB0_q   <- round(B_qk[,1],2)
  empMSY_q  <- round(MSY_q,2)
  
  Uusr_q    <- eqListQuant$Uusr_q
  YeqUSR_q  <- eqListQuant$YeqUSR_q
  
  maxUidx <- length(U_qk[2,])
  
  maxU <- 1.3*max(Ucrash_q)
  
  qIdx <- which(qProbs %in% range(qProbs))
  medIdx <- which(qProbs == 0.5)
  
  B_itk <- empRefCurves$B_isptk[,1,1,,]
  
  if (area == 'SOG'){
    orderedsim <- c(1,6)
    U <- c(0,0.2)
    if(history == T){
      YLIM <- c(0,26)
      XLIM <- c(0,110)
    }
    else{
      YLIM <- c(0,8)
      XLIM <- c(0,90)
    }
    adjustpos <- 2.4
    msy_x = 5
  }
  else if (area == 'HG'){
    orderedsim <- c(1,6)
    U <- c(0,0.14)
    YLIM <- c(0, 0.25)
    adjustpos <- 0.2
    msy_x = .5
  }
  
  
  label <- c(expression(B[0]),expression("0.3"*B[0]),expression(B['MSY | OM 1']),
             expression("0.4"*B['MSY | OM 1']),expression("0.8"*B['MSY | OM 1']),
             "USR | OM 1")
  
  plotLayout <- matrix( c(1),
                        ncol = 1, byrow = T)
  layout(plotLayout)
  
  par(mar = c(.25,.25,.25,.25),
      oma = c(4,4,1,1), xaxs = 'r', yaxs = 'r', las = 1)
  
  
  yrs <- 2024:2523
  clrs <- paletteer_d("RColorBrewer::Dark2",length(label))
  
  # plot yield
  plot( x = c(-0.3,maxB), y = c(0,0.25),ylim = YLIM,xlim = XLIM,
        axes = FALSE, type = "n" )
  axis(side = 2, las =1, cex.axis = 1)
  axis(side = 1, cex.axis = 1)
  grid()
  box()
  
  lines( x = B_qk[medIdx,1:maxUidx], y = Y_qk[medIdx,1:maxUidx], col = "black", lwd = 3, lty = 1)
  
  fit <- smooth.spline(x = B_qk[medIdx,maxUidx:1], y = Y_qk[medIdx,maxUidx:1], nknots = 10)
  # lines(fit, col = "black", lwd = 3)
  
  # B0
  points(x = empB0_q[2], y = predict(fit, empB0_q[2])$y, pch = 16, col = clrs[1], cex = 1.5)
  text(x = empB0_q[2], y = predict(fit, empB0_q[2])$y, label[1], pos = 2)
  
  # 0.3B0
  segments(x0 = predict(fit, 0.3*empB0_q[2])$x, y0 = -1, 
           x1 = predict(fit, 0.3*empB0_q[2])$x, y1 = predict(fit, 0.3*empB0_q[2])$y, lty = 2)
  points(predict(fit, 0.3*empB0_q[2]), pch = 16, col = clrs[2], cex = 1.5)
  text(x = 0.3*empB0_q[2]-adjustpos, y = 0, label[2], pos = 4, srt = 90)
  
  # MSY and label
  segments(x0 = -5, y0 = empMSY_q[2], x1 = empBmsy_q[2], y1 = empMSY_q[2], lty = 2)
  text(x = msy_x, y = empMSY_q[2], "MSY | OM 1", pos = 3)
  
  
  # Bmsy
  segments(x0 = empBmsy_q[2], y0 = -1, x1 = empBmsy_q[2], y1 = empMSY_q[2], lty = 2)
  points(x = empBmsy_q[2], y = empMSY_q[2], pch = 16, col = clrs[3], cex = 1.5)
  text(x = empBmsy_q[2]-adjustpos, y = 0, label[3], pos = 4, srt = 90)
  
  # 0.4Bmsy
  segments(x0 = 0.4*empBmsy_q[2], y0 = -1, 
           x1 = 0.4*empBmsy_q[2], y1 = predict(fit, 0.4*empBmsy_q[2])$y, lty = 2)
  points(predict(fit, 0.4*empBmsy_q[2]), pch = 16, col = clrs[4], cex = 1.5)
  text(x = 0.4*empBmsy_q[2]-adjustpos, y = 0, label[4], pos = 4, srt = 90)
  
  # 0.8Bmsy
  segments(x0 = 0.8*empBmsy_q[2], y0 = -1, 
           x1 = 0.8*empBmsy_q[2], y1 = predict(fit, 0.8*empBmsy_q[2])$y, lty = 2)
  points(predict(fit, 0.8*empBmsy_q[2]), pch = 16, col = clrs[5], cex = 1.5)
  text(x = 0.8*empBmsy_q[2]-adjustpos, y = 0, label[5], pos = 4, srt = 90)
  
  # USR
  segments(x0 = USR_q[2], y0 = -1, 
           x1 = USR_q[2], y1 = predict(fit, USR_q[2])$y, lty = 2)
  points(x = USR_q[2], y = predict(fit, USR_q[2])$y, col = clrs[6], pch = 16, cex = 1.5)
  text(x = USR_q[2]-adjustpos, y = 0, label[6], pos = 4, srt = 90, cex = 0.8)
  
  if(history == T){
    # historic data: 2004 - 2023
    endYr <- 2023
    startYr <- endYr - histYrs +1
    B_it <- array(0,dim = c(200,histYrs))
    Y_it <- array(0,dim = c(200,histYrs))
    load(paste0("Outputs/",folder,"/sim_parBat",folder,1,"/sim_parBat",folder,1,".RData"))
    startIdx <- 73 - histYrs +1
    B_it <- blob$om$SB_ispt[,1,1,startIdx:73]
    Y_it <- blob$om$C_ispt[,1,1,startIdx:73]
    
    B_t <- apply(B_it, MARGIN = 2, FUN = median)
    catch <- read.csv('history/catchData.csv')
    group <- group_by(catch, Year) %>% summarise(C = sum(Value)) %>% 
      filter(Year >= startYr)
    divClrs <- paletteer_c('grDevices::Grays', n = length(B_t))
    
    lines(B_t, group$C, col = 'gray', lwd = 2)
    points(B_t, group$C, col = 'black', bg = divClrs, pch = 21, )
    text(B_t[1],group$C[1], startYr, pos = 1)
    text(B_t[histYrs],group$C[histYrs], endYr, pos = 1)
  }
  
  
  mtext( side = 1, outer = TRUE, text = "Spawning Biomass (kt)",
         line = 2.5, font = 2, las = 0)
  mtext(side = 2, text = "Yield (kt)", line = 2.5, font = 2, las = 0)
  
  
}# plotRefPtsExpC()
