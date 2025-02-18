// ********************************************************************
// Module: SISCA.cpp
// Authors: S. D. N. Johnson, S.P. Cox (LFR)
// Procedure: Spatially Integrated Statistical Catch at Age model
//
// References: 
//
// Revised from single fleet to multifleet by S. D. N. Johnson
// Date Last Revised: March 3, 2018
// 
// Spatially Integrated Statistical Catch at Age model. Developed
// by Landmark Fisheries Research for assessment of spatially
// diverse finfish stocks with discrete fisheries. Can
// support movement for each age class, or model stocks
// as independent sub-populations. Shrinkage priors are applied
// to improve estimates of biological and observation model
// parameters (steepness, maybe catchability).
// 
// Features:
//  1. Discrete fisheries (Pope's approx)
//  2. Hierarchical shrinkage priors for bio and obs model pars
//  3. Uncentred hierarchical distributions for easy identification
//      of pars among stock replicates and better tmbStan 
//      performance.
//  4. Logistic normal likelihood for compositional data, with
//      optional correlated residuals at lag-1.
//  5. Conditional maximum likelihood estimates for catchability
//      are used when appropriate. 
//  6. Single sex
//  7. Average weight-at-age observations by fleet.
// 
// Style/usage Notes for future authors:
//  1. All arrays are given as a variable name followed by
//      subscripts after an underscore, e.g. biomass for 
//      stock p and time t is B_pt
//  2. Observation error variances are tau, process error variances
//      are sigma
//  3. Order of subscripts: apgt (age, pop, gear, time)
// 
// LFR provides this code with ABSOLUTELY NO WARRANTY. This is source
// code for custom software, i.e. this is always the beta version.
//      
// 
// NOTES for SDNJ:
// 1. What can we refactor into functions? Mixed and specific data?
// *******************************************************************

#include <TMB.hpp>       // Links in the TMB libraries
#include <iostream>

// posfun
template<class Type>
Type posfun(Type x, Type eps, Type &pen)
{
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2)-x/eps));
}

// // invLogit
// template<class Type>
// Type invLogit(Type x, Type scale, Type trans)
// {
//   return scale/(Type(1.0) + exp(-Type(1.0)*x)) - trans;
// }

// square
template<class Type>
Type square(Type x)
{
  return x*x;
}


// calcSpawn()
// When called depletes numbers at age to
// the appropriate level for the spawn timing, calculates
// spawning biomass and produces (age-1) recruitment for
// the following time step.
// inputs:    SB_pt       = array of spawning biomas by pop/time
//            spawnN_ap   = array of numbers at age by pop at spawning time
//            M_apt       = array of natural mortality by age
//            Mjuve       = Scalar of juvenile mortality rate
//            wtAge_ap    = array of weight at age by pop
//            mat         = vector proportion mature at age
//            R_pt        = array of recruitment
//            omegaR_pt   = array of recruitment process errors
//            prevTime    = The fractional time step of last event (fishing/movement)
//            t           = current time step
//            tInitModel_p= Initial time step of model for pop p
// outputs:   tmpN_ap     = updated array of post-spawning numbers at age.
// Usage:     For applying a movement hypothesis in a multi-stock system
// Source:    S. D. N. Johnson
// Reference: FIND ONE
// Side-effects:  - prevTime is updated to be moveTime,
//                - SB_pt is updated with the spawning biomass for time t
//                - R_pt is updated with recruitment at time t+1
template<class Type>
void calcSpawn( array<Type>   spawnN_ap,
                array<Type>&  SB_pt,
                array<Type>   M_ap,
                array<Type>   wtAge_ap,
                vector<Type>  mat,
                array<Type>&  R_pt,
                array<Type>&  Eggs_pt,
                Type          fec,
                vector<Type>  reca_p,
                vector<Type>  recb_p,
                int           t,
                vector<int>   tInitModel_p,
                int           SRindVar )
{
  // Get model dims
  int nP = SB_pt.dim(0);

  // Deplet tmpN_ap by natural mortality
  for( int p = 0; p < nP; p++ )
  {
    // Caclulate spawning biomass at spawn time
    SB_pt(p,t)    = (spawnN_ap.col(p) * wtAge_ap.col(p) * mat ).sum();
    // Calculate eggs in millions (fec is eggs per g, so scale
    // by 1e3 to get eggs per kg)
    Eggs_pt(p,t)  = SB_pt(p,t) * fec * 0.5 * 1e3;

    // Update expected BH recruitment vector
    if( t >= tInitModel_p(p) )
    {
      if( SRindVar == 2)
        R_pt(p,t+1)  = reca_p(p)*Eggs_pt(p,t)/( 1. + recb_p(p)*Eggs_pt(p,t) );
      if( SRindVar == 1)
        R_pt(p,t+1)  = reca_p(p)*SB_pt(p,t)/( 1. + recb_p(p)*SB_pt(p,t) );
    }
  }
} // END calcSpawn()


// projSpawnBio()
// When called, depletes numbers at age by natural mortality 
// and a given level of catch in each fleet, producing
// a projected spawning biomass ay spawn timing. Assumes 
// deterministic recruitment, and that all fleets are
// discrete timing (to simplify matters)
// inputs:    N_ap          = array of start of year numbers-at-age
//            M_ap          = array of natural mortality by age
//            wtAge_ap      = array of weight at age by pop
//            mat           = vector proportion mature at age
//            fleetTiming_g = fleet timing
//            chronIdx_g    = Chronological order of fleets
//            sel_apg       = Selectivity of each fleet
//            alloc_pg      = Allocation of catch among fleets
//            C_v           = vector of total catch levels
//            postPondM_g   = Post-ponding mortality rate
// outputs:   catTable_pvk  = array that is p copies of a v x k mtx:
//                            Catch | SSB | Dep | HR
// Usage:     For projecting biomass to the end of the following 
//            year under a vector of catches C_v (basically a catch table)
// Source:    S. D. N. Johnson
// Side-effects:  - prevTime is updated to be moveTime,
//                - SB_pt is updated with the spawning biomass for time t
//                - R_pt is updated with recruitment at time t+1
template<class Type>
void projCatchBio(  array<Type>   initN_ap,
                    array<Type>   M_ap,
                    array<Type>   wtAge_ap,
                    array<Type>   wtAge_agp,
                    vector<Type>  mat_a,
                    vector<Type>  fleetTiming_g,
                    vector<int>   fleetType_g,
                    vector<int>   chronIdx_g,
                    Type          spawnTiming,
                    array<Type>   sel_apg,
                    array<Type>   alloc_pg,
                    vector<Type>  C_v,
                    vector<Type>  postPondM_g,
                    vector<Type>  B0_p,
                    array<Type>&  catTable_pvk,
                    int           t)
{
  // Get model dims
  int nA = initN_ap.dim(0);
  int nP = initN_ap.dim(1);
  int nG = sel_apg.dim(2);
  int nV = C_v.size();

  // Loop over the catch vector
  for( int v = 0; v < nV; v++ )
  {
    // First, convert catch to fishery specific
    // catches
    array<Type> C_pg(nP,nG);
    C_pg = C_v(v) * alloc_pg;

    // Need to loop over fisheries, deplete by M and
    // then remove catch. This might be a useful function
    // for further refactoring...
    // need an array of HRs and vuln bio
    array<Type>   U_pg(nP,nG);
    array<Type>   vulnB_pg(nP,nG);
    vector<Type>  spawnB_p(nP);
    U_pg.setZero();
    vulnB_pg.setZero();
    spawnB_p.setZero();

    // Deplete initN_ap by natural mortality and
    // fishing mortality
    for( int p = 0; p < nP; p++ )
    {
      for(int gg = 0; gg < nG; gg++ )
      {
        int g = chronIdx_g(gg);

        // tmp vector of num-at-age
        vector<Type> tmpN_a(nA);
        tmpN_a = initN_ap.col(p);

        // remove prev catch
        if(gg > 0)
          for( int ggg = 0; ggg < gg; ggg++)
            tmpN_a *= (1 - U_pg(p,ggg) * sel_apg.col(ggg).col(p));



        // Now deplete by M
        tmpN_a *= exp(-1 * fleetTiming_g(g) * M_ap.col(p));


        // Calculate vuln bio and 
        vulnB_pg(p,g) = (tmpN_a * wtAge_agp.col(p).col(g) * sel_apg.col(g).col(p)).sum();
        U_pg(p,g) = C_pg(p,g) / vulnB_pg(p,g); 


      }

      // Now do spawn time numbers
      vector<Type> spawnN_a = initN_ap.col(p);
      // And ponded fish
      vector<Type> pond_a(nA);
      pond_a.setZero();

      // deplete by fishing
      for( int g = 0; g < nG; g++ )
      {
        
        if( fleetType_g(g) == 2 | fleetType_g(g) == 3 )
          pond_a += U_pg(p,g) * sel_apg.col(g).col(p) * spawnN_a * exp(-1 * (fleetTiming_g(g) * M_ap.col(p) + postPondM_g(g)));

        if( fleetTiming_g(g) < spawnTiming)
          spawnN_a *= (1 - U_pg(p,g) * sel_apg.col(g).col(p));


      }

      // Deplete by M, add surviving ponded fish back in
      spawnN_a *= exp( -1 * spawnTiming *M_ap.col(p));
      spawnN_a += pond_a;

      // Add up biomass
      spawnB_p(p) = (spawnN_a * mat_a * wtAge_ap.col(p)).sum();

      // Now fill catTable_pvk
      catTable_pvk(p,v,0) = C_v(v);
      catTable_pvk(p,v,1) = spawnB_p(p);
      catTable_pvk(p,v,2) = spawnB_p(p)/B0_p(p);
      catTable_pvk(p,v,3) = C_v(v)/(spawnB_p(p) + C_v(v));
    }
      
  }
} // END projCatchBio()


// applyMovement()
// When called applies a family of Markov movement 
// matrices to the temporary state variable tmpN_ap 
// (numbers at age by pop) within the primary population
// dynamics process.
// inputs:    M_apt     = Type array of natural mortality by age
//            tmpN_ap   = Type array of numbers at age by pop
//            mov_ppa   = Type array of movement matrices (each a slice is a p x p matrix)
//            prevTime  = The most recent fractional time step of an event (fishing, spawning)
//            moveTime  = The fractional time step where movement is applied
// outputs:   tmpN_ap   = updated array of post-movement numbers at age.
// Usage:     For applying a movement hypothesis in a multi-stock system
// Source:    S. D. N. Johnson
// Reference: FIND ONE
// Side-effects: prevTime is updated to be moveTime.
template<class Type>
array<Type> applyMovement(  array<Type>   M_ap,
                            array<Type>   tmpN_ap,
                            array<Type>   mov_ppa,
                            array<Type>&  initN_ap,
                            array<Type>&  termN_ap,
                            Type&         prevTime,
                            Type          moveTiming )
{
  // Pull dimensions
  int nA = tmpN_ap.dim(0);
  int nP = tmpN_ap.dim(1);

  // Apply mortality
  // Get fraction of time step since prevTime
  Type fracM = 0;
  fracM = moveTiming - prevTime;

  // Now pull movement matrix for each age
  // class and apply to the P-vector
  // of individual numbers at age
  for( int a = 0; a < nA; a++ )
  {
    matrix<Type> mov_pp(nP,nP);
    mov_pp = mov_ppa.col(a).matrix();

    vector<Type> initNum_p(nP);
    vector<Type> termNum_p(nP);
    initNum_p.setZero();
    termNum_p.setZero();

    for( int p = 0; p < nP; p++ )
    {
      initNum_p(p) = tmpN_ap(a,p) * exp( -fracM * M_ap(a,p) );
      initN_ap(a,p) = initNum_p(p);
    }


    termNum_p = mov_pp * initNum_p;

    // Now place terminal numbers back into
    // the tmpN_ap array
    for( int p = 0; p < nP; p ++)
    {
      termN_ap(a,p) = termNum_p(p);
      tmpN_ap(a,p) = termNum_p(p);
    }

    prevTime = moveTiming;
  }

  // Return new numbers at age
  return(termN_ap);

} // END applyMovement()

// convertSOK()
// Converts Spawn-on-Kelp catch to ponded fish (in numbers)
// given values of various conversion parameters.
// inputs:    mat_a   = maturity at age vector
//            sel_a   = selectivity-at-age vector
//            C       = catch of SOK
//            pEff    = proportion effective
//            pFem    = proportion female
//            vN_a    = vector of numbers at age
//            propMat = proportion mature
//            fec     = fecundity
//            wtAge_a = weight-at-age
//            gamma   = SOK conversion factor (eggs to SOK)
//            psi     = ratio of biomass of ponded fish to SOK
template<class Type>
Type convertSOK(  Type          C,
                  vector<Type>  mat_a, 
                  vector<Type>  sel_a,
                  Type          pEff,
                  Type          pFem,
                  vector<Type>  vN_a,
                  Type&         propMat,
                  Type          fec,
                  array<Type>   wtAge_a,
                  Type          gamma,
                  Type&         psi,
                  Type          initF )
{
  // We are generating ponded fish from SOK, so
  // we need to work backwards from the ISCAM example...
  // Need to use proportion mature of population
  // rather than ponded fish, slight bias as function
  // of f, but within 5 basis points under a single
  // fleet fishery, not variable wrt age structure
  int nA = vN_a.size();
  vector<Type> pondC_a(nA);
  vector<Type> Z_a(nA);

  // Temp variable of mature biomass
  Type tmpBmat = 0;
  Type tmpB = 0;

  pondC_a.setZero();

  for( int a = 0; a < nA; a++ )
  {
    if( sel_a(a) > 0 & mat_a(a) > 0 )
    {
      Z_a(a) = initF * sel_a(a);
      pondC_a(a) = (1 - exp(-Z_a(a))) * vN_a(a) * initF/Z_a(a);
    }

    tmpBmat  += pondC_a(a) * mat_a(a) * wtAge_a(a);
    tmpB     += pondC_a(a) * wtAge_a(a);
  }

  // Update proportion mature
  propMat = tmpBmat / tmpB;

  // Calculate psi
  psi = pEff * pFem * gamma * fec * propMat;

  // Now convert SOK to ponded fish
  Type pondedBio = C / psi;

  return(pondedBio);
} // END convertSOK()

// calcLogistNormLikelihood()
// Calculates the logistic normal likelihood for compositional data.
// Automatically accumulates proportions below a given threshold. Takes
// cumulative sum of squared resids and number of classes for which
// residuals are calculated as pointers, so that these can be accumulated
// across years for conditional MLE of variance. 
// Will extend to correlated residuals and time-varying variance later.
// inputs:    yObs    = Type vector of observed compositions (samples or proportions)
//            pPred   = Type vector of predicted parameters (true class proportions)
//            minProp = minimum proportion threshold to accumulate classes above
//            etaSumSq= cumulative sum of squared 0-mean resids
//            nResids = cumulative sum of composition classes (post-accumulation)
// outputs:   resids = vector of resids (accumulated to match bins >= minProp)
// Usage:     For computing the likelihood of observed compositional data
// Source:    S. D. N. Johnson
// Reference: Schnute and Haigh, 2007; Francis, 2014
template<class Type>
vector<Type> calcLogistNormLikelihood(  vector<Type>& yObs, 
                                        vector<Type>& pPred,
                                        Type minProp,
                                        Type& etaSumSq,
                                        Type& nResids )
{
  // Get size of vector 
  int nX = yObs.size();

  // Normalise the samples and predictions in case they are numbers
  // and not proportions
  yObs /= yObs.sum();
  pPred /= pPred.sum();

  // Create vector of residuals to return
  vector<Type> resids(nX);
  vector<Type> aboveInd(nX);
  resids.setZero();
  aboveInd.setZero();

  // Count number of observations less
  // than minProp
  int nAbove = 0;
  for( int x = 0; x < nX; x++)
    if(yObs(x) >= minProp)
    {
      nAbove++;
      aboveInd(x) = 1;
    }

  // Now loop and fill
  vector<Type> newObs(nAbove);
  vector<Type> newPred(nAbove);
  newObs.setZero();
  newPred.setZero();
  int k = 0;
  for( int x = 0; x < nX; x++)
  {
    // accumulate observations
    newObs(k) += yObs(x);
    newPred(k) += pPred(x);

    // Increment k if we reach a bin
    // with higher than minProp observations
    // and we aren't yet in the last bin
    if(yObs(x) >= minProp & k < nAbove - 1)
      k++;    
  }

  // Create a residual vector
  vector<Type> res(nAbove);
  res.setZero(); 
  
  for( int k = 0; k < nAbove; k++)
    if( newObs(k) > 0  &  newPred(k) > 0)
      res(k) = log(newObs(k)) - log(newPred(k));

  // Calculate mean residual
  Type meanRes = res.sum()/nAbove;
  // centre residuals
  res -= meanRes;


  // Now loop again, and fill in
  // the residuals vector
  // for plotting later
  k = 0;
  for( int x =0; x < nX; x++)
    if( aboveInd(x) == 1)
    {
      resids(x) = res(k);
      k++;
    }

  
  // Now add squared resids to etaSumSq and nRes to nResids
  etaSumSq  += (res*res).sum();
  nResids   += nAbove;

  return(resids);
} // end calcLogistNormLikelihood()


// calcCorrLogistNormLikelihood()
// Calculates the logistic normal likelihood for compositional data.
// Automatically accumulates proportions below a given threshold. Takes
// cumulative sum of squared resids and number of classes for which
// residuals are calculated as pointers, so that these can be accumulated
// across years for conditional MLE of variance. 
// Will extend to correlated residuals and time-varying variance later.
// inputs:    yObs    = Type vector of observed compositions (samples or proportions)
//            pPred   = Type vector of predicted parameters (true class proportions)
//            minProp = minimum proportion threshold to accumulate classes above
//            etaSumSq= cumulative sum of squared resids (possibly correlated)
//            nResids = cumulative sum of composition classes (post-accumulation)
//            meanSampSize = mean size of annual composition samples
//            compLikeFun = switch for compositional likelihood from Francis 2014
//                            0 => no correlation (AR0, LN1)
//                            1 => AR1 model (LN2)
//                            2 => AR2 model (LN3)
//                            3 => ARMA model (LN3m)
//            phi1    = Correlation parameter, depends on compLikeFun
//            psi     = Correlation parameter, depends on compLikeFun
//            compObsNLL = compositional observation NLL contribution (external var)
// outputs:   resids = vector of standardised resids (accumulated to match bins >= minProp)
// Usage:     For computing the likelihood of observed compositional data
// Source:    S. D. N. Johnson
// Reference: Schnute and Haigh, 2007; Francis, 2014
template<class Type>
vector<Type> calcCorrLogistNormLikelihood(  vector<Type>   yObs, 
                                            vector<Type>   pPred,
                                            Type           minProp,
                                            Type&          etaSumSq,
                                            Type&          nResids,
                                            Type           meanSampSize,
                                            int            compLikeFun,
                                            matrix<Type>   Corr_aa,
                                            Type&          compObsNLL,
                                            Type&          compWt,
                                            vector<Type>&  tcComps,
                                            vector<Type>&  tcPred,
                                            Type&          gmObs,
                                            Type&          gmPred,
                                            int&           nBins,
                                            Type&          logdetV,
                                            matrix<Type>&  saveVchk,
                                            matrix<Type>&  saveCorr,
                                            matrix<Type>&  saveKmat,
                                            matrix<Type>&  saveHmat,
                                            matrix<Type>&  saveFmat,
                                            matrix<Type>&  saveGamma )
{

  // NOTE: some of the notation departs from
  // that defined in the model header. This is
  // to keep it close to Francis 2014 while I work
  // out the model

  // Get size of vector 
  int nX = yObs.size();

  // Get size of this sample
  Type thisSampSize = yObs.sum();

  // Calculate this year's weighting
  Type Wy = sqrt(meanSampSize / thisSampSize);
  compWt = Wy;

  // Normalise the samples and predictions in case they are numbers
  // and not proportions
  yObs /= yObs.sum();
  pPred /= pPred.sum();

  // Create vector of residuals to return
  vector<Type> resids(nX);
  vector<Type> aboveInd(nX);
  resids.setZero();
  aboveInd.setZero();

  // Count number of observations less
  // than minProp
  int nAbove = 0;
  for( int x = 0; x < nX; x++)
    if(yObs(x) >= minProp)
    {
      nAbove++;
      aboveInd(x) = 1;
    }

  nBins = nAbove - 1;

  // If there is only 1 bin after tail compression
  // then there is no info in this year. Move on.
  if(nAbove == 1)
  {
    vector<Type> resids(nX);
    resids.setZero();
    etaSumSq = 0;
    nResids = 0;
    compObsNLL = 0;
    compWt = 1;

    return(resids);
  }

  

  vector<int> idxAbove(nAbove);
  int k = 0;
  for( int x =0; x < nX; x++ )
    if( aboveInd(x) == 1)
    {
      idxAbove(k) = x;
      k++;
    }

  // Now loop and fill
  vector<Type> newObs(nAbove);
  vector<Type> newPred(nAbove);
  newObs.setZero();
  newPred.setZero();

  // Need to add a switch for external or
  // internal tail compression, then
  // need to compress the predictions based on 
  // how we're doing it.

  // Run tail compression below.
  // now go up from the left
  k = 0;
  for( int x = 0; x < nX; x++)
  {
    // accumulate observations
    newObs(k) += yObs(x);
    newPred(k) += pPred(x);

    // Increment k if we reach a bin
    // with higher than minProp observations
    // and we aren't yet in the last bin
    if( (yObs(x) >= minProp) & (k < nAbove - 1) )
      k++;    
  }

  
  // Both should produce 2 vectors of
  // newObs and newPred with small observations
  // compressed

  // OK, now we create a nAbove x nAbove correlation matrix Corr
  // from the above correlation vector, the same
  // size matrix V, and the multiplier matrix K
  matrix<Type> Corr(nAbove,nAbove);
  matrix<Type> V(nAbove,nAbove);
  matrix<Type> Vinv(nAbove,nAbove);
  matrix<Type> Gamma(nAbove,nAbove);
  matrix<Type> K(nAbove-1,nAbove);
  matrix<Type> F(nAbove-1,nAbove);
  matrix<Type> H(nAbove-1,nAbove-1);
  matrix<Type> Hinv(nAbove-1,nAbove-1);


  // Fill with zeroes
  Corr.setZero();
  V.setZero();
  Vinv.setZero();
  K.setZero();
  F.setZero();
  H.fill(1);

  // Get submatrix of Corr_aa
  for( int k = 0; k < nAbove; k++ )
  {
    Corr(k,k) += Corr_aa(idxAbove(k),idxAbove(k));
    for( int kk = k + 1; kk < nAbove; kk++ )
    {
      Corr(k,kk) += Corr_aa(idxAbove(k),idxAbove(kk));
      Corr(kk,k) += Corr_aa(idxAbove(kk),idxAbove(k));
    }
  }
  

  // Now we fill
  for( int rIdx = 0; rIdx < nAbove; rIdx ++ )
  {

    // Now fill in K
    if( rIdx < nAbove-1 )
    {
      K( rIdx, rIdx ) = 1;
      F( rIdx, rIdx ) = 1;
      H( rIdx, rIdx) += 1;
    }

    if( rIdx == nAbove-1 )
    {
      K.col(rIdx).fill(-1);
      F.col(rIdx).fill(1);
    }

    
  }

  // Generate V and its inverse
  V = K * Corr * K.transpose();
  Vinv = atomic::matinvpd( V, logdetV );
  Hinv = atomic::matinv( H );

  Gamma = F.transpose() * Hinv * V * Hinv * F;

  // Fill saveVchk
  saveVchk.topLeftCorner(nAbove-1,nAbove-1) = V;
  saveCorr.topLeftCorner(nAbove,nAbove)     = Corr;
  saveKmat.topLeftCorner(nAbove-1,nAbove)   = K;
  saveHmat.topLeftCorner(nAbove-1,nAbove-1) = H;
  saveFmat.topLeftCorner(nAbove-1,nAbove)   = F;
  saveGamma.topLeftCorner(nAbove,nAbove)    = Gamma;

  // Calculate w_b, independent resids
  matrix<Type> w_b(nAbove-1,1);
  w_b.setZero();
  for( int vIdx = 0; vIdx < nAbove-1; vIdx ++)
    if( newObs(vIdx) > 0 & newPred(vIdx) > 0)
      w_b(vIdx,0) = log( newObs(vIdx) / newObs(nAbove-1) ) - 
                  log( newPred(vIdx) / newPred(nAbove-1) );

  
  // Now add correlated squared resids to etaSumSq and nRes to nResids
  matrix<Type> tmpeta(1,1);
  tmpeta = (w_b.transpose() * Vinv * w_b) / square(Wy); 
  etaSumSq  += tmpeta(0,0);
  nResids   += nAbove - 1;

  compObsNLL += 0.5 * logdetV + (nAbove-1) * log(Wy); // + log( newObs.segment(0,nAbove -1) ).sum();

  // Now expand residuals back out
  // Now loop again, and fill in
  // the residuals vector
  // for plotting later
  vector<Type> res(nAbove);
  res.setZero();
  gmObs  =  exp( log(newObs).sum()  / nAbove  );
  gmPred =  exp( log(newPred).sum() / nAbove  );
  // gmPred =  pow(newPred.prod(),1/nAbove);

  for( int k = 0; k < nAbove; k++ )
  {
    res    =  log( newObs / gmObs ) - 
              log( newPred / gmPred );  
  }
  
  k = 0;
  for( int x =0; x < nX; x++)
    if( aboveInd(x) == 1)
    {
      resids(x)   = res(k) / Wy / sqrt(Gamma(k,k)) ;
      tcComps(x)  = newObs(k) * thisSampSize;
      tcPred(x)   = newPred(k);
      // For AR2 model (not working yet)
      // if( compLikeFun > 1)
      //   resids(x) /=  Wy * sqrt(Gamma(k,k)) ;
      k++;
    }

  return(resids);
}  // END calcCorrLogistNormLikelihood()


// calcFleetOrder()
// Refactored out of main objective function
// so we can apply every year and adjust ordering.
template<class Type>
void calcFleetOrder(  vector<Type> fleetTiming_g,
                      vector<int>& chronIdx_g,
                      vector<int>  mortType_g )
{
  int nG = fleetTiming_g.size();

  // Loop through all fleets in 2 nested loops, and build a 
  // new vector of indices in chronological order
  Type          minTime = 0;
  Type          prevFleetTime = 0;
  vector<int>   usedFleet(nG);

  usedFleet.fill(int(0));
  chronIdx_g.fill(int(0));

  // Initialise fleet timing M fracs
  Type fracM = 0.;
  Type lastFrac = 0.;

  for( int gg = 0; gg < nG; gg++ )
  {
    // Only calculate timing/order for discrete fisheries
    
    // Set max time to end of year (hiT), and
    // lowG (idx with next smallest time) as the first fleet
    Type        hiT = 1;
    int         lowG = 0;
    // Loop over fleet timing vector, pull
    // get gear idx of minimum year fraction
    for( int gIdx = 0; gIdx < nG; gIdx++ )
    {
      // if( usedFleet(gIdx) == 1 )
      //   next();

      if( ( fleetTiming_g(gIdx) >= minTime) & 
          ( fleetTiming_g(gIdx) <= hiT ) &
          ( usedFleet(gIdx) == 0) &
          ( mortType_g(gIdx) == 0) )
      {
        // Record index of new lowest time that hasn't been used
        lowG      = gIdx;
        hiT       = fleetTiming_g(gIdx);
      }

      if(mortType_g(gIdx) == 1)
        usedFleet(gIdx) = int(1);
    }


    chronIdx_g(gg)  = lowG;
    usedFleet(lowG) = int(1);
    prevFleetTime   = minTime;
    minTime         = fleetTiming_g(lowG);



  }

} // END calcFleetOrder()



// <><><><><> Objective Function <><><><><> //
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Call namespaces //
  using namespace density;

  /*\/\/\/\/\ DATA SECTION /\/\/\/\*/
  // Data Structures
  DATA_ARRAY(I_pgt);                    // CPUE data by pop, gear, time
  DATA_ARRAY(C_pgt);                    // Catch data by pop, gear, time
  DATA_ARRAY(A_apgt);                   // Age observations in table format (age, pop, gear, time)
  DATA_ARRAY(W_apgt);                   // Observed weight-at-age by pop, gear, time (gear specific to estimate better numerical removals)
  DATA_ARRAY(W_apt);                    // Observed weight-at-age by pop, time (not gear specific for pop dynamics)
  DATA_ARRAY(mI_gt);                    // Observation data for mixed stocks by gear and time
  DATA_ARRAY(mC_gt);                    // Catch data for mixed stocks by gear and time (likely SOK)
  DATA_ARRAY(mA_agt);                   // Age observations in table format for mixed stocks, by age, gear, and time
  DATA_ARRAY(rI_pgt);                   // Proportion of total idx by pop, gear, time
  DATA_ARRAY(combI_pt);                 // A single index that combines the surface and dive surveys 
  DATA_ARRAY(diluteN_apt);              // An array of input numbers-at-age to dilute the catch from particular fleets/predators
  DATA_ARRAY(E_pgt);                    // Effort/abundance for gear/predator g at time t in area p (when mixed catch, use the same effort for each p)

  // To cut down to tables, we require
  // nP, nT, nG and nA included as data
  // Model dimensions
  int nA = W_apgt.dim(0);
  int nP = W_apgt.dim(1);
  int nG = W_apgt.dim(2);
  int nT = W_apgt.dim(3);

  
  // Model switches
  DATA_IVECTOR(survType_g);             // Type of index (0 = vuln bio, 1 = vuln numbers, 2 = spawn bio)
  DATA_IVECTOR(indexType_g);            // Type of survey (0 = relative, 1 = absolute)
  DATA_IARRAY(deltaIdx_pg);             // Apply a delta model for 0s in indices (0 == NO, 1 == YES)
  DATA_INTEGER(deltaIndVar);            // Independent variable for delta logistic regression (1 == Bio, 2 == Depletion)
  DATA_IVECTOR(qPrior_g);               // Prior on catchability? (0 == No, 1 == Yes)
  DATA_IVECTOR(calcIndex_g);            // Calculate fleet index (0 = no, yes = 1)
  DATA_INTEGER(condMLEtauObs);          // Estimate index observation error as a free par (0) or as the conditional MLE (1)
  DATA_IVECTOR(selType_g);              // Type of selectivity (0 = asymptotic, 1 = domed (normal), 2 = decreasing asymptotic, 3 = domed (dbl asymptotic))
  DATA_ARRAY(scaleSel_gt);              // Selectivity curve scalars (used for certain predator selectivity)
  DATA_INTEGER(hierSel);                // Switch for hierarchical shrinkage prior (1) or single-level prior (0) in selectivity models
  DATA_IVECTOR(selX_g);                 // Len or Age base selectivity (0 = age, 1 = length)
  DATA_IVECTOR(tvSelFleets);            // Switch for time-varying selectivity or not
  DATA_VECTOR(ageCompWeight_g);         // Age composition likelihood weight scalar
  DATA_VECTOR(idxLikeWeight_g);         // observation idx likelihood weighting scalar
  DATA_VECTOR(catLikeWeight_g);         // catch observation idx likelihood weighting scalar
  DATA_VECTOR(fleetTiming_g);           // Fraction of year before catch/survey observation
  DATA_IVECTOR(fleetType_g);            // Survey (0), Fishery/Predator (1), CP spawn (2), spawn OP or pred (3)
  DATA_IVECTOR(mortType_g);             // Type of fishing mortality (0 = discrete, 1 = continuous)
  DATA_IVECTOR(catSeriesType_g);        // Removals (0) or effort (1)
  DATA_INTEGER(minAgeDilute);           // Min age (0-indexed) for dilution
  DATA_SCALAR(alphaDilute);             // Proportion of diluteN_apt to add to modeled stock for removals  
  DATA_VECTOR(betaDilute_g)             // Proportion of year that fleet overlaps with diluted portion of the stock (i.e., proportion of catch/consumption taken from mixed stock)
  DATA_IVECTOR(initCode_p);             // initialise at 0 => unfished, 1=> fished
  DATA_IVECTOR(initFcode_p);            // Use Finit or not (0 ==> No, g > 0 ==> Use fleet g-1 selectivity)
  DATA_IVECTOR(initRcode_p);            // Use Rinit or not (0 ==> No, 1 ==> Use Rinit)
  DATA_STRING(initMethod);              // Fished initialisation method (see header, options are surv and nums)
  DATA_IVECTOR(tInitModel_p);           // Time step for initialising population dynamices model for pop p
  DATA_SCALAR(posPenFactor);            // Positive-penalty multiplication factor
  DATA_IVECTOR(firstRecDev_p);          // First recruitment deviation by pop
  DATA_IVECTOR(lastRecDev_p);           // Last rec dev by pop
  DATA_SCALAR(omegaRradius);            // Radius of bounded omegaR values
  DATA_SCALAR(minPropAge);              // Minimum proportion in age comps (< minPropage are accumulated)
  DATA_IVECTOR(minAge_g);               // minimum age to be considered in age comp data for gear g
  DATA_INTEGER(nYearsProjAve);          // number of years to average M and weight-at-age over for projected biomass
  DATA_SCALAR(spawnTiming);             // fraction of year that passes before spawning
  DATA_SCALAR(moveTiming);              // Fraction of year that passes before movement
  DATA_SCALAR(useMovement);             // Switch for applying movement matrix (1 == yes, 0 == no)
  DATA_INTEGER(juveMage);               // Integer max age for juvenile M
  DATA_INTEGER(juveMsource);            // Integer switch for juve M source (0 == Mbar, 1 == input/est)
  DATA_IVECTOR(avgRcode_p);             // Integer switch for using average R recruitment (1) or BH recruitment (0)
  DATA_INTEGER(SRindVar);               // Stock-recruit model based on 1 == Biomass, 2 == Eggs
  DATA_IVECTOR(whichCombIdx_g);         // Indicator of which surveys are combined in combI_pt
  DATA_INTEGER(densityDepM);            // Indicator of whether tvM is RW and DD
  DATA_INTEGER(corrMdevs);              // Indicator of correlated M devs (1 == corr, 0 == uncorr)
  DATA_INTEGER(corrRdevs);              // Indicator of correlated R devs (1 == corr, 0 == uncorr)
  DATA_SCALAR(corrParWeight);           // Correlation matrix parameter penalty weight - might need a prior here
  DATA_IVECTOR(mixComps_g);             // Mix compositional data by averaging over years with data, using sample size as weighting for expected values
  DATA_IVECTOR(catAllocYrs);            // start and end time step for catch allocation calculations
  DATA_SCALAR(dataLikeWt);              // Data likelihood weighting
  DATA_SCALAR(priorDensWt);             // Prior density weighting

  // 1-year ahead catch tables
  DATA_VECTOR(C_v);                     // Vector of future catches to estimate the following year's depletion level
  // DATA_VECTOR(allocYrRange);            // Vector of beginning and end years for allocation

  // SOK values
  DATA_SCALAR(fec);                     // Eggs per fish
  DATA_VECTOR(gamma_g);                 // Conversion factor (# of eggs to weight of landings/consumption)
  DATA_SCALAR(pFem);                    // Proportion of ponded fish that are female
  DATA_VECTOR(postPondM_g);             // Post ponding mortality rate for each gear
  DATA_SCALAR(sokInitF);                // Initial F used for estimating propMature
  DATA_IVECTOR(calcPropEff_t);          // switch vector for calculating proportion effective (0 => no, 1 => yes)

  // Prior weights
  DATA_SCALAR(jeffWtB0);                // Weighting on B0 Jeffries prior
  DATA_SCALAR(lnm1PriorWt);             // weighting on lnm1_p jeffries prior
  DATA_SCALAR(meanRdevWt);              // weighting on lnm1_p jeffries prior

  // Compositional likelihood inputs
  DATA_INTEGER(compLikeFun);            // Compositional data likelihood function (correlation)
  DATA_ARRAY(meanSampSize_pg);          // Mean sample sizes for stock-specific compositional data
  DATA_VECTOR(meanSampSize_g);          // Mean sample sizes for mixed stock compositional data


  /*\/\/\/\/\ PARAMETER SECTION /\/\/\/\*/
  // Leading biological parameters //
  PARAMETER_VECTOR(lnB0_p);             // Unfished spawning biomass for each stock
  PARAMETER_VECTOR(lnRinit_p);          // Fished initialisation recruitment for each stock
  PARAMETER_VECTOR(lnRbar_p);           // Average recruitment for each stock
  PARAMETER(logit_ySteepness);          // SR steepness on (0,1)
  PARAMETER(lnM);                       // Natural mortality rates
  PARAMETER(lnMjuve);                   // Natural mortality rate for juveniles
  PARAMETER(lnm1);                      // mean DDM rate
  PARAMETER_VECTOR(epslnm1_p);          // DDM rate by area - may fix equal among areas.
  PARAMETER_ARRAY(fDevs_ap);            // Non-eq initial numbers multiplier (nAxnP array)
  PARAMETER_VECTOR(lnFinit_p);          // Non-eq initial F.

  // Random stock effects on biological pars
  PARAMETER_VECTOR(epsM_p);             // Stock specific M deviation
  PARAMETER(lnsigmaStockM);             // Stock M effect log-SD
  PARAMETER_VECTOR(epsSteep_p);         // Stock specific steepness deviation on logit scale
  PARAMETER(lnsigmaStockSteep);         // Stock steepness effect log-SD

  // Observation models //
  // Selectivity
  PARAMETER_VECTOR(lnSelAlpha_g);       // log scale selectivity alpha par (ageSel50, lenSel50, or mu in N model)
  PARAMETER_VECTOR(lnSelBeta_g);        // log scale selectivity beta par (age/len sel step, or sigma in N model)
  PARAMETER_ARRAY(epsSelAlpha_pg);      // log-normal deviation in x-at-50%/mode sel for each stock
  PARAMETER_ARRAY(epsSelBeta_pg);       // log-normal deviation in x-at-95%/std. dev sel for each stock
  PARAMETER_VECTOR(epsSelAlpha_vec);    // log-normal error in x-at-50%/mode sel
  PARAMETER_VECTOR(epsSelBeta_vec);     // log-normal error in x-at-95%/std dev sel
  PARAMETER_VECTOR(lndSelAlpha_g);      // log scale selectivity alpha par (ageSel50, lenSel50, or mu in N model)
  PARAMETER_VECTOR(lndSelBeta_g);       // log scale selectivity beta par (age/len sel step, or sigma in N model)
  PARAMETER_ARRAY(epsdSelAlpha_pg);     // log-normal deviation in x-at-50%/mode sel for each stock
  PARAMETER_ARRAY(epsdSelBeta_pg);      // log-normal deviation in x-at-95%/std. dev sel for each stock
  PARAMETER_VECTOR(epsdSelAlpha_vec);   // log-normal error in x-at-50%/mode sel
  PARAMETER_VECTOR(epsdSelBeta_vec);    // log-normal error in x-at-95%/std dev sel
  PARAMETER_VECTOR(lnsigmaSelAlpha_g);  // log-SD of SelAlpha deviations
  PARAMETER_VECTOR(lnsigmaSelBeta_g);   // log-SD of SelBeta deviations
  PARAMETER(lnsigmaTVsel);              // SD on deviations in time-varying selectivity

  // Fishery removals
  PARAMETER_ARRAY(lnqFinit_pg);         // Array of initial catchabilities to scale effort to removals
  PARAMETER_ARRAY(deltaqF_pgt);         // Array of log-normal random jumps for qF
  PARAMETER_VECTOR(resCatCV_g);         // Vector of catch observation CVs


  // Survey obs error
  PARAMETER_ARRAY(lntau2Obs_pg);        // Explicitly model observation error variance (stock spec. indices)
  PARAMETER_VECTOR(lntau2Obs_g);        // Explicitly model observation error variance (mixed stock indices)
  PARAMETER_ARRAY(lntauObsComb_pg);     // survey-specific obs error SD for blended index


  // Process errors //
  PARAMETER_ARRAY(recDevs_pt);          // Recruitment log-deviations by area/time
  PARAMETER(lnsigmaR);                  // log scale recruitment SD          
  PARAMETER(priorSigR);                 // Prior mean on sigmaR (sigmaR/sqrt2*sqrtpi)
  PARAMETER_ARRAY(omegaM_pt);           // Natural mortality log-deviations
  PARAMETER(lnsigmaM);                  // log scale natural mortality SD
  PARAMETER_VECTOR(off_diag_R);         // Entries for cholesky factor of recruitment corr matrix
  PARAMETER_VECTOR(off_diag_M);         // Entries for cholesky factor of mortality corr matrix

  // Priors //
  PARAMETER_VECTOR(obstau2IGa);         // Inverse Gamma Prior alpha for stock index obs error var prior
  PARAMETER_VECTOR(obstau2IGb);         // Inverse Gamma Prior beta for stock index obs error var prior
  PARAMETER_VECTOR(sig2RPrior);         // Hyperparameters for sig2R prior - (IGa,IGb) or (mean,var)
  PARAMETER_VECTOR(sig2MPrior);         // Hyperparameters for sig2M prior - (IGa,IGb) or (mean,var)
  PARAMETER_VECTOR(rSteepBetaPrior);    // pars for Beta dist steepness prior
  PARAMETER_VECTOR(Mprior);             // pars for M prior (mean,sd), initial M for RWM, and M0 for DDM
  PARAMETER_VECTOR(m1Prior);            // pars for density dependent m1 prior (mean,sd)
  PARAMETER_VECTOR(mlnq_g);             // Hierarchical prior mean catchability
  PARAMETER_VECTOR(sdlnq_g);            // Hierarchical prior sd catchability
  PARAMETER_VECTOR(mq);                 // prior mean catchability 
  PARAMETER_VECTOR(sdq);                // prior mean catchability 
  PARAMETER_ARRAY(lnqComb_pg);          // survey specific log-catchability for use in combined idx
  PARAMETER_VECTOR(mlnqF_g);            // Fleet specific prior mean on removal catchability
  PARAMETER_VECTOR(sdlnqF_g);           // Fleet specific prior sd on initial removal catchability
  PARAMETER_VECTOR(sddeltaqF_g);        // Fleet specific prior SD on removal catchability jumps


  // Fixed LH parameters - if provided
  PARAMETER_VECTOR( mat_a );            // Maturity at age
  PARAMETER_VECTOR( fec_a );            // fecundity-at-age
  PARAMETER( inputL1 );                 // Input length-at-age 1 to overwrite the vonB model
  PARAMETER( Linf );                    // Asymptotic length
  PARAMETER( L1 );                      // Length at age 1
  PARAMETER( vonK );                    // vonB growth rate
  PARAMETER_VECTOR( lenWt );            // Allometric length/weight c1, c2 parameters

  // Max selectivity age
  PARAMETER_VECTOR(mlnSelAlpha_g);      // prior mean selectivity Alpha
  PARAMETER_VECTOR(mlnSelBeta_g);       // prior mean selectivity Alpha
  PARAMETER_VECTOR(sdSel_g);            // selectivity prior SD
  PARAMETER_VECTOR(selPriorWt_g);       // selectivity prior weighting

  // Markov movement model - a square matrix for every
  // age group read in as a nP x nP x nA dim array 
  PARAMETER_ARRAY(mov_ppa);             // Markov movement matrix

  // Closed ponding Spawn-on-Kelp (SOK) parameters and variables
  
  PARAMETER_VECTOR(propEffBounds);      // 2-vector of upper and lower bounds on propEff
  PARAMETER_VECTOR(logitPropEff_vec);   // Vector of logit transformed proportion effective - will be distributed to the right years/stocks
  PARAMETER(mPsi);                      // Average scalar converting between ponded fish and SOK
  PARAMETER(sdPsi);                     // SD on scalar converting between ponded fish and SOK

  // Compositional data likelihood correlation matrix parameters
  PARAMETER_VECTOR(logitphi1_g);        // AR1 correlation coefficient, or AR2 parameter
  PARAMETER_VECTOR(logitpsi_g);         // AR2 parameter

  // DLN bernoulli model parameters
  PARAMETER_ARRAY(lnSDProbPosIdx_pg);    // SD of probability of detecting spawn indices
  PARAMETER_ARRAY(meanProbPosIdx_pg);    // Mean probability of detecting spawn indices
  PARAMETER_VECTOR(muSDProbPosIdx_g);    // Prior mean for the delta LN SD variable
  PARAMETER_VECTOR(muMeanProbPosIdx_g);  // Prior mean for the delta LN Mean variable
  PARAMETER_VECTOR(sigSDProbPosIdx_g);   // Prior mean for the delta LN SD variable
  PARAMETER_VECTOR(sigMeanProbPosIdx_g); // Prior mean for the delta LN Mean variable

  // Convex combo DeltaLN pars
  // PARAMETER_VECTOR(lnSDProbPosCombIdx_p);// SD of probability of detecting spawn indices
  // PARAMETER_VECTOR(meanProbPosCombIdx_p);// Mean probability of detecting spawn indices
  // PARAMETER(muSDProbPosCombIdx);         // Prior mean for the delta LN SD variable
  // PARAMETER(muMeanProbPosCombIdx);       // Prior mean for the delta LN Mean variable
  // PARAMETER(sigSDProbPosCombIdx);        // Prior mean for the delta LN SD variable
  // PARAMETER(sigMeanProbPosCombIdx);      // Prior mean for the delta LN Mean variable

  // Derived Variables //
  // Transformed scalar parameters
  vector<Type>  B0_p            = exp(lnB0_p);
  Type          M               = exp(lnM);
  Type          sigmaR          = exp(lnsigmaR);
  Type          sigmaM          = exp(lnsigmaM);
  Type          sigmaStockM     = exp(lnsigmaStockM);
  Type          sigmaStockSteep = exp(lnsigmaStockSteep);
  vector<Type>  M_p             = exp( lnM + sigmaStockM * epsM_p );
  Type          ySteepness      = invlogit(logit_ySteepness);
  Type          rSteepness      = 0.2 + 0.78 * ySteepness;

  // "unfished" mortality
  vector<Type>  M0_p(nP);
  M0_p.setZero();
  
  // Initial Fs
  vector<Type>  Finit_p(nP);
  Finit_p.setZero();

  vector<Type>  Rinit_p(nP);
  Rinit_p.setZero();

  vector<Type>  Rbar_p(nP);
  Rbar_p = exp(lnRbar_p);


  vector<Type>  Mjuve_p(nP);
  Mjuve_p.setZero();
  if( juveMsource == 1)
    Mjuve_p    += exp(lnMjuve);
  

  // Correlation matrices for Rec and Mort
  // Make unstructured corr
  vector<Type> MoffDiag(nP * (nP-1)/2);
  MoffDiag.fill(0);
  MoffDiag += off_diag_M;

  vector<Type> RoffDiag(nP*(nP-1)/2);
  RoffDiag.fill(0);
  RoffDiag += off_diag_R;
  UNSTRUCTURED_CORR_t<Type> tvM_mvnStd(MoffDiag);
  UNSTRUCTURED_CORR_t<Type> SR_mvnStd(RoffDiag);

  // Pull correlation matrix
  matrix<Type> corrR_pp(nP,nP);
  matrix<Type> corrM_pp(nP,nP);
  corrM_pp.setIdentity();
  corrR_pp.setIdentity();

  if( corrMdevs == 1)
    corrM_pp = tvM_mvnStd.cov();
  if( corrRdevs == 1)
    corrR_pp = SR_mvnStd.cov();


  matrix<Type> cholR_pp(nP,nP);
  matrix<Type> cholM_pp(nP,nP);
  cholR_pp.setIdentity();
  cholM_pp.setIdentity();

  matrix<Type> sigmaR_pp(nP,nP);
  matrix<Type> sigmaM_pp(nP,nP);
  sigmaR_pp.setIdentity();
  sigmaM_pp.setIdentity();

  sigmaR_pp *= sigmaR;
  sigmaM_pp *= sigmaM;

  // Now fill correlation matrices
  int k = 0;
  for( int i=1; i<nP; i++)
  {
    Type Norm2_R=cholR_pp(i,i);
    Type Norm2_M=cholM_pp(i,i);
    for(int j=0; j<=i-1;j++)
    {
      cholR_pp(i,j) = off_diag_R(k);
      Norm2_R += cholR_pp(i,j)*cholR_pp(i,j); 

      cholM_pp(i,j) = off_diag_M(k++);
      Norm2_M += cholM_pp(i,j)*cholM_pp(i,j); 
    }
    for(int j=0; j<=i; j++)
    {
      cholR_pp(i,j) /= sqrt(Norm2_R);
      cholM_pp(i,j) /= sqrt(Norm2_M);
    }
  }
  // // Make correlation matrices
  // corrR_pp = cholR_pp * cholR_pp.transpose();
  // corrM_pp = cholM_pp * cholM_pp.transpose();

  // Make covariance matrices for output
  matrix<Type> SigmaM_pp = sigmaM_pp * corrM_pp * sigmaM_pp;
  matrix<Type> SigmaR_pp = sigmaR_pp * corrR_pp * sigmaR_pp;



  // Compositional likelihood correlation matrix
  // parameters
  vector<Type> phi1_g(nG);
  phi1_g.fill(0);
  vector<Type> psi_g(nG);
  psi_g.fill(.5);
  vector<Type> phi2_g(nG);
  phi2_g.fill(0);

  // Compute correlations out to nX in case we need
  // them
  // Start by defining the LN1 correlation vector
  array<Type> corr_ga(nG,nA);
  corr_ga.fill(0);
  corr_ga.col(0) += 1;
  // Now fill in the rest: 
  if( compLikeFun > 0)
  {
    // The LN2 function
    if( compLikeFun == 1 )
    {
      phi1_g = 2 * invlogit(logitphi1_g) - 1;
      for( int a = 1; a < nA; a++ )
      {
        corr_ga.col(a) = corr_ga.col(a-1) * phi1_g;
      }
    }  

    // The LN3 function
    if( compLikeFun == 2 )
    { 
      phi1_g = 4 * invlogit(logitphi1_g) - 2;
      psi_g  = invlogit(logitpsi_g);
      phi2_g = -1 + (2 - sqrt(square(phi1_g)) ) * psi_g;
      
      corr_ga.col(1) = phi1_g / phi2_g;
      for( int a = 2; a < nA; a++ )
        corr_ga.col(a) = phi1_g * corr_ga.col(a-1) + phi2_g * corr_ga.col(a-2);
    }

    // LN3m, ARMA model
    if( compLikeFun == 3 )
    {
      phi1_g    = 2 * invlogit(logitphi1_g) - 1;
      psi_g    += logitpsi_g;

      vector<Type> sumPhiPsi = phi1_g + psi_g;
      
      corr_ga.col(1) = phi1_g + psi_g * (1 + (sumPhiPsi * sumPhiPsi )/(1 - phi1_g * phi1_g ) );
      
      for( int a = 2; a < nA; a++ )
        corr_ga.col(a) = phi1_g * corr_ga.col(a-2);

    }
  }

  // Now make correlation matrices
  array<Type> Corr_gaa(nG,nA,nA);
  Corr_gaa.setZero();

  for( int g = 0; g < nG; g++ )
  {
    for( int rIdx = 0; rIdx < nA; rIdx ++ )
    {
      Corr_gaa(g,rIdx,rIdx) = 1;
      for( int cIdx = rIdx + 1; cIdx < nA; cIdx ++)
      {
        int idxDiff = cIdx - rIdx;

        Corr_gaa(g,rIdx,cIdx) = corr_ga(g,idxDiff);
        Corr_gaa(g,cIdx,rIdx) = corr_ga(g,idxDiff);

      }
    }
  }
  
  // Make stock-specific steepness pars
  vector<Type>  ySteepness_p(nP);
  for( int p = 0; p < nP; p++)
    ySteepness_p(p) = invlogit(logit_ySteepness + sigmaStockSteep*epsSteep_p(p) );
  vector<Type>  rSteepness_p  = 0.2 + 0.78 * ySteepness_p;

  // Transformed Selectivity model parameter vectors
  vector<Type>  SelAlpha_g       = exp(lnSelAlpha_g);
  vector<Type>  SelBeta_g        = exp(lnSelBeta_g);
  vector<Type>  dSelAlpha_g      = exp(lndSelAlpha_g);
  vector<Type>  dSelBeta_g       = exp(lndSelBeta_g);
  vector<Type>  sigmaSelAlpha_g  = exp(lnsigmaSelAlpha_g);
  vector<Type>  sigmaSelBeta_g   = exp(lnsigmaSelBeta_g);
  Type          sigmaTVsel       = exp(lnsigmaTVsel);

  array<Type>   tau2Obs_pg(nP,nG);
  array<Type>   tauObs_pg(nP,nG);
                tau2Obs_pg        = exp(lntau2Obs_pg);
                tauObs_pg         = exp(.5 * lntau2Obs_pg);
  vector<Type>  tauObs_g          = exp(0.5*lntau2Obs_g);
  vector<Type>  tau2Obs_g         = exp(lntau2Obs_g);
                
  array<Type>   SDProbPosIdx_pg(nP,nG);
                SDProbPosIdx_pg   = exp(lnSDProbPosIdx_pg);
  
  
  vector<Type>  scaleSel_g(nG);
  scaleSel_g.fill(1);

  // Transformed parameters in arrays require loops to fill them
  // Initial N multipliers
  int nInitAges = 0;
  if( initMethod == "surv" )
    nInitAges = nA - 1;
  if( initMethod == "nums" )
    nInitAges = nA;
  
  array<Type> initN_mult_ap(nInitAges,nP);
  // Time-varying selectivity
  array<Type> SelAlpha_pgt(nP,nG,nT);
  array<Type> SelBeta_pgt(nP,nG,nT);
  array<Type> SelAlpha_pg(nP,nG);
  array<Type> SelBeta_pg(nP,nG);
  array<Type> dSelAlpha_pgt(nP,nG,nT);
  array<Type> dSelBeta_pgt(nP,nG,nT);
  array<Type> dSelAlpha_pg(nP,nG);
  array<Type> dSelBeta_pg(nP,nG);

  // Loop and fill
  for( int g = 0; g < nG; g++)
  {
    if( hierSel == 1 )
    {
      SelAlpha_pg.col(g)          = SelAlpha_g(g) * exp( sigmaSelAlpha_g(g) * epsSelAlpha_pg.col(g) );
      SelBeta_pg.col(g)           = SelBeta_g(g) * exp( sigmaSelBeta_g(g) * epsSelBeta_pg.col(g) );

      dSelAlpha_pg.col(g)         = dSelAlpha_g(g) * exp( sigmaSelAlpha_g(g) * epsdSelAlpha_pg.col(g) );
      dSelBeta_pg.col(g)          = dSelBeta_g(g) * exp( sigmaSelBeta_g(g) * epsdSelBeta_pg.col(g) );
    }

    if( hierSel == 0 )
    {
      SelAlpha_pg.col(g)          = exp( epsSelAlpha_pg.col(g));
      SelBeta_pg.col(g)           = exp( epsSelBeta_pg.col(g));

      dSelAlpha_pg.col(g)         = exp( epsdSelAlpha_pg.col(g));
      dSelBeta_pg.col(g)          = exp( epsdSelBeta_pg.col(g));

    }

    scaleSel_g(g) = scaleSel_gt.transpose().col(g).mean();
  }

  SelAlpha_pgt.col(0)  = SelAlpha_pg;
  SelBeta_pgt.col(0)   = SelBeta_pg;

  dSelAlpha_pgt.col(0) = dSelAlpha_pg;
  dSelBeta_pgt.col(0)  = dSelBeta_pg;

  // initial rec devs for fished init
  initN_mult_ap.fill(1);
  for( int p = 0; p < nP; p++)
  {
    if( initCode_p(p) == 1)
      initN_mult_ap.col(p) = exp(sigmaR *fDevs_ap.col(p));

    if( initFcode_p(p) > 0 )
      Finit_p(p) = exp(lnFinit_p(p));

  }

  // Loop over time steps and gears, create a matrix
  // of deviations in selectivity
  array<Type> epsSelAlpha_pgt(nP,nG,nT);
  array<Type> epsSelBeta_pgt(nP,nG,nT);
  epsSelAlpha_pgt.setZero();
  epsSelBeta_pgt.setZero();

  array<Type> epsdSelAlpha_pgt(nP,nG,nT);
  array<Type> epsdSelBeta_pgt(nP,nG,nT);
  epsdSelAlpha_pgt.setZero();
  epsdSelBeta_pgt.setZero();


  int epsSelVecIdx = 0;
  for( int g = 0; g < nG; g++ )
  {
    // Scale fleet average selectivity pars by scaleSel_g
    SelAlpha_g(g) *= scaleSel_g(g);
    SelBeta_g(g) *= scaleSel_g(g);

    dSelAlpha_g(g) *= scaleSel_g(g);
    dSelBeta_g(g) *= scaleSel_g(g);

    for( int t = 0; t < nT; t++ )
    { 
      // Bring previous time-step selecivity
      // forward
      if( t > 0 )
      {
        SelAlpha_pgt.col(t).col(g) = SelAlpha_pgt.col(t-1).col(g);
        SelBeta_pgt.col(t).col(g) = SelBeta_pgt.col(t-1).col(g);

        dSelAlpha_pgt.col(t).col(g) = dSelAlpha_pgt.col(t-1).col(g);
        dSelBeta_pgt.col(t).col(g) = dSelBeta_pgt.col(t-1).col(g);
      
        // Update the epsilon array
        for( int p = 0; p < nP; p++ )
          if( (tvSelFleets(g) == 1) & (A_apgt(0,p,g,t) >= 0) )
          {
            epsSelAlpha_pgt(p,g,t) += sigmaTVsel * epsSelAlpha_vec(epsSelVecIdx);
            epsSelBeta_pgt(p,g,t) += sigmaTVsel * epsSelBeta_vec(epsSelVecIdx);

            epsdSelAlpha_pgt(p,g,t) += sigmaTVsel * epsdSelAlpha_vec(epsSelVecIdx);
            epsdSelBeta_pgt(p,g,t) += sigmaTVsel * epsdSelBeta_vec(epsSelVecIdx);
            epsSelVecIdx++;
          }    
        // Now update SelAlpha_gt and SelBeta_gt - can happen
        // at first time step, since we're using the gear specific
        // selectivity curve when we don't age observations
        SelAlpha_pgt.col(t).col(g)  *= exp(epsSelAlpha_pgt.col(t).col(g));
        SelBeta_pgt.col(t).col(g)   *= exp(epsSelBeta_pgt.col(t).col(g));

        dSelAlpha_pgt.col(t).col(g) *= exp(epsdSelAlpha_pgt.col(t).col(g));
        dSelBeta_pgt.col(t).col(g)  *= exp(epsdSelBeta_pgt.col(t).col(g));
      }

    }
  }

  for( int p = 0; p < nP; p++)
  {
    for( int g = 0; g < nG; g++ )
      for( int t = 0; t < nT; t++ )
      {
        // Now scale by externally supplied scalar
        SelAlpha_pgt(p,g,t)  *= scaleSel_gt(g,t);
        SelBeta_pgt(p,g,t)   *= scaleSel_gt(g,t);

        dSelAlpha_pgt(p,g,t) *= scaleSel_gt(g,t);
        dSelBeta_pgt(p,g,t)  *= scaleSel_gt(g,t);
      }
  }


  // Stock recruit model
  vector<Type> reca_p(nP);              // BH a par
  vector<Type> recb_p(nP);              // BH b par
  vector<Type> phi_p(nP);               // ssbpr
  vector<Type> R0_p(nP);                // unfished eqbm recruitment
  reca_p.setZero();
  recb_p.setZero();
  phi_p.setZero();
  R0_p.setZero();


  // DDM model
  vector<Type> totB0_p(nP);                                   // total start year biomass
  Type m1 = exp(lnm1);                               // mean DDM rate
  vector<Type> m1_p(nP);
  for( int p = 0; p < nP; p++)
    m1_p(p) = exp(lnm1 + sigmaStockM * epslnm1_p(p));    // DD rate for M_p
  totB0_p.setZero();

  // Recruitment devs
  array<Type> omegaR_pt(nP,nT+1);
  omegaR_pt.setZero();
  // int recVecIdx = 0;
  for( int p = 0; p < nP; p++ )
    for(int t = firstRecDev_p(p); t < lastRecDev_p(p); t++)
    {
      omegaR_pt(p,t) = -1*omegaRradius + 2*omegaRradius / (1 + exp(-recDevs_pt(p,t-1)));
    }

  // Optimisation and modeling quantities
  Type objFun   = 0.;                         // Objective function
  Type posPen   = 0.;                         // posFun penalty
  Type rec_nlp  = 0.;                         // recruitment deviations NLL
  Type init_nlp = 0.;                         // Initial rec deviations for fished init
  Type mort_nlp = 0.;                         // lnM prior
  Type Mdev_nlp = 0.;                         // Hierarchical M deviations prior
  Type lnm1_nlp = 0.;                         // lnm1 mean and dev prior
  Type h_nlp    = 0.;                         // Steepness NLP
  Type hDev_nlp = 0.;                         // Hierarchical steepenss deviations prior
  Type tvMnlp   = 0.;                         // Mt deviations correlated MVN NLP
  Type SRnlp    = 0.;                         // SR deviations correlated MVN NLP
  array<Type>   obsIdxNLL_pg(nP,nG);          // observation indices NLL - LN part
  array<Type>   obsIdxDeltaNLL_pg(nP,nG);     // observation indices NLL - Delta part
  vector<Type>  obsCombIdxDeltaNLL_p(nP);     // convex combo indices NLL - Delta part
  array<Type>   ageObsNLL_pg(nP,nG);          // Age observation NLL
  array<Type>   tau2Age_pg(nP,nG);            // Age observation variance
  array<Type>   etaSumSq_pg(nP,nG);           // Sum of squared age logist resids
  array<Type>   nResids_pg(nP,nG);            // Number of residuals
  array<Type>   intrnlAgeLikCpt_pg(nP,nG);    // Internally calculated contribution to the composition obs likelihood
  array<Type>   nObsAge_pg(nP,nG);            // Number of discrete years with observations
  array<Type>   nlptau2idx_pg(nP,nG);         // Observation error var prior
  array<Type>   ageResids_apgt(nA,nP,nG,nT);  // age observation residuals
  array<Type>   probPosIdx_pgt(nP,nG,nT);     // Probability of a positive spawn idx observation
  array<Type>   probPosCombIdx_pt(nP,nT);     // Probability of a positive convex combo spawn idx observation

  // Mixed data
  vector<Type>  obsMixedIdxNLL_g(nG);         // NLL for mixed indices
  vector<Type>  obsCombIdxNLL_p(nP);          // NLL for area specific combined survey idx
  vector<Type>  tau2Age_g(nG);                // Age observation variance for mixed ages
  array<Type>   ageResids_agt(nA,nG,nT);      // age observation resids for mixed ages
  array<Type>   intrnlAgeLikCpt_g(nG);        // Internally calculated contribution to the composition obs likelihood
  array<Type>   ageObsNLL_g(nG);              // Mixed age observations NLL
  vector<Type>  nlptau2idx_g(nG);             // Observation error var prior density (mixed obs)

  // Catch observation likelihood (mortType == 1)
  array<Type>   catObsNLL_pg(nP,nG);
  array<Type>   priorqF_pg(nP,nG);
  array<Type>   catRes_pgt(nP,nG,nT);

  // Age comp debugging info
  array<Type>   ageWt_pgt(nP,nG,nT);          // stock-specific sample size based weights for age comps in each year
  array<Type>   ageWt_gt(nG,nT);              // mixed age comp sample size based weights in each year 
  array<Type>   tcComps_apgt(nA,nP,nG,nT);    // tail compressed stock specific age comps for debugging likelihood
  array<Type>   tcComps_agt(nA,nG,nT);        // tail compressed mixed age comps for debugging likelihood
  array<Type>   tcPred_apgt(nA,nP,nG,nT);     // tail compressed stock specific age comps for debugging likelihood
  array<Type>   tcPred_agt(nA,nG,nT);         // tail compressed mixed age comps for debugging likelihood
  array<Type>   gmObs_pgt(nP,nG,nT);          // Geometric mean of observations
  array<Type>   gmPred_pgt(nP,nG,nT);         // Geometric mean of predictions
  array<Type>   logdetV_pgt(nP,nG,nT);        // Geometric mean of predictions
  array<int>    nBins_pgt(nP,nG,nT);          // Number of observations in each year
  array<Type>   Vchk_aapgt(nA,nA,nP,nG,nT);   // Internally calculated matrix for likelihood comps
  array<Type>   Corr_aapgt(nA,nA,nP,nG,nT);   // Internally calculated matrix for likelihood comps
  array<Type>   K_aapgt(nA,nA,nP,nG,nT);      // Internally calculated matrix for likelihood comps
  array<Type>   H_aapgt(nA,nA,nP,nG,nT);      // Internally calculated matrix for likelihood comps
  array<Type>   F_aapgt(nA,nA,nP,nG,nT);      // Internally calculated matrix for likelihood comps
  array<Type>   Gamma_aapgt(nA,nA,nP,nG,nT);  // Internally calculated matrix for likelihood comps

  // Initialise vectors at 0
  catObsNLL_pg.setZero();
  catRes_pgt.setZero();
  priorqF_pg.setZero();
  obsIdxNLL_pg.fill(0.0);
  obsIdxDeltaNLL_pg.fill(0.0);
  ageObsNLL_pg.fill(0.0);
  tau2Age_pg.fill(0.0);
  nlptau2idx_pg.fill(0.0);
  nlptau2idx_g.fill(0.0);
  ageResids_apgt.setZero();
  obsMixedIdxNLL_g.setZero();
  tau2Age_g.setZero();
  ageResids_agt.setZero();
  ageObsNLL_g.setZero();
  intrnlAgeLikCpt_g.setZero();
  intrnlAgeLikCpt_pg.setZero();
  ageWt_pgt.setZero();
  ageWt_gt.setZero();
  tcComps_apgt.fill(-1);
  tcComps_agt.fill(-1);
  tcPred_apgt.fill(-1);
  tcPred_agt.fill(-1);
  gmObs_pgt.setZero();
  gmPred_pgt.setZero();
  nBins_pgt.setZero();
  logdetV_pgt.setZero();
  Corr_aapgt.setZero();
  K_aapgt.setZero();
  H_aapgt.setZero();
  F_aapgt.setZero();
  probPosIdx_pgt.setZero();
  probPosCombIdx_pt.setZero();
  obsCombIdxDeltaNLL_p.setZero();
  obsCombIdxNLL_p.setZero();


  // Life history schedules
  vector<Type>  age_a(nA);
  vector<Type>  wt_a(nA);
  vector<Type>  len_a(nA);
  array<Type>   surv_ap(nA,nP);
  array<Type>   initSurv_ap(nA,nP);

  // Selectivity values
  array<Type>   sel_ag(nA,nG);  
  sel_ag.setZero();
  // Time-varying selectivity
  array<Type>   sel_apgt(nA,nP,nG,nT);
  sel_apgt.setZero();
  


  // State variables
  array<Type>   N_apt(nA,nP,nT+1);          // Numbers at age
  array<Type>   tmpN_apt(nA,nP,nT+1);       // Numbers at age
  array<Type>   endN_apt(nA,nP,nT+1);       // Numbers at age
  array<Type>   spawnN_apt(nA,nP,nT+1);     // Numbers at age at spawn timing
  array<Type>   movN_apt(nA,nP,nT+1);       // Numbers at age after movement
  array<Type>   initN_apt(nA,nP,nT+1);      // Numbers at age after movement
  array<Type>   termN_apt(nA,nP,nT+1);      // Numbers at age after movement
  array<Type>   B_apt(nA,nP,nT+1);          // Total biomass at age (beginning of year)
  array<Type>   B_pt(nP,nT+1);              // Total biomass (beginning of year)
  array<Type>   vulnB_apgt(nA,nP,nG,nT);    // Vulnerable biomass at age by gear
  array<Type>   vulnB_pgt(nP,nG,nT);        // Vulnerable biomass by gear
  array<Type>   vulnN_apgt(nA,nP,nG,nT);    // Vulnerable numbers at age by gear
  array<Type>   vulnN_pgt(nP,nG,nT);        // Vulnerable numbers by gear
  array<Type>   SB_pt(nP,nT+1);             // Spawning biomass
  array<Type>   Eggs_pt(nP,nT+1);           // Eggs by stock and time
  array<Type>   bhR_pt(nP,nT+1);            // Expected BH recruitments at age-1
  array<Type>   R_pt(nP,nT+1);              // Actualy age-1 recruitments
  array<Type>   M_apt(nA,nP,nT+1);          // Natural mortality array (age, stock, time)
  array<Type>   Z_apt(nA,nP,nT+1);          // Total mortality array (age, stock, time)
  array<Type>   appZ_apt(nA,nP,nT+1);       // Total mortality array (age, stock, time)
  array<Type>   U_pgt(nP,nG,nT);            // Fleet harvest rate
  array<Type>   U_apgt(nA,nP,nG,nT);        // Fleet harvest rate
  array<Type>   F_pgt(nP,nG,nT);            // Instantaneous fishing mortality rate (mortType == 1)
  array<Type>   qF_pgt(nP,nG,nT);           // Fishing/predation mortality per unit effort/predator (mortType == 1)
  array<Type>   lnqF_pgt(nP,nG,nT);         // log-Fishing/predation mortality per unit effort/predator (mortType == 1)
  array<Type>   totSurv_apt(nA,nP,nT);      // Total survival after M and fishing.
  array<Type>   uAge_apgt(nA,nP,nG,nT);     // Vuln proportion at age in each gear at each time
  array<Type>   catAge_apgt(nA,nP,nG,nT);   // Catch at age in each gear at each time step
  // array<Type>   catAge_apgt(nA,nP,nG,nT);   // Catch at age in each gear at each time step
  array<Type>   predPA_apgt(nA,nP,nG,nT);   // Predicted proportion at age in each gear at each time step
  array<Type>   totC_pgt(nP,nG,nT);         // Total removals for a time/gear/pop, adding mixed catches
  array<Type>   expC_pgt(nP,nG,nT);         // Expected total removals for a time/gear/pop (including mixed catches), for mortType_g == 1
  array<Type>   pondC_apgt(nA,nP,nG,nT);    // Ponded fish (numbers) at age/pop/gear/time
  array<Type>   pondC_pgt(nP,nG,nT);        // total ponded fish (biomass) by pop/gear/time
  array<Type>   pondC_gt(nG,nT);            // total mixed ponded fish (biomass) by gear/time
  array<Type>   combIhat_gt(nG,nT);         // Predicted combined state variables being indexed by mixed indices
  array<Type>   psi_pgt(nP,nG,nT);          // Ponded biomas <-> SOK factor by pop/gear/time
  array<Type>   psi_gt(nG,nT);              // Mixed ponded biomas <-> mixed SOK factor by gear/time
  vector<Type>  movFracM_t(nT);             // movement time fraction of M

  // Diluted consumption by predators
  array<Type>   vulnBdiluted_apgt(nA,nP,nG,nT);     // Vulnerable biomass at age by gear
  array<Type>   vulnBdiluted_pgt(nP,nG,nT);         // Vulnerable biomass by gear
  array<Type>   vulnNdiluted_apgt(nA,nP,nG,nT);     // Vulnerable numbers at age by gear
  array<Type>   vulnNdiluted_pgt(nP,nG,nT);         // Vulnerable numbers by gear
  array<Type>   totCundiluted_pgt(nP,nG,nT);        // Total catch from undiluted stock
  array<Type>   totCdiluted_pgt(nP,nG,nT);          // Total catch from diluted stock
  array<Type>   uAgeDiluted_apgt(nA,nP,nG,nT);      // Vuln proportion at age in each gear at each time
  array<Type>   catAgeUndiluted_apgt(nA,nP,nG,nT);  // Catch-at-age from undiluted portion
  array<Type>   catAgeDiluted_apgt(nA,nP,nG,nT);    // catch-at-age from diluted portion
  array<Type>   propCatAgeDiluted_apgt(nA,nP,nG,nT);// Proportion of diluted catch-at-age that comes from main stock 

  // Diluted catch-at-age given a mixing stock
  array<Type>   totCatAge_apgt(nA,nP,nG,nT);
  array<Type>   totCpredilute_pgt(nP,nG,nT);

  // We need to weight proportion caught at age
  // by the proportion of vulnerable fish
  // in each stock later, so we need arrays to
  // hold that info
  array<Type>   mixedVulnN_agt(nA,nG,nT);   // Proportion of vulnerable numbers in each age class
  array<Type>   mixedVulnN_gt(nG,nT);       // Proportion vulnerable in each stock
  array<Type>   mixedVulnB_agt(nA,nG,nT);   // Proportion of vulnerable numbers in each age class
  array<Type>   mixedVulnB_gt(nG,nT);       // Proportion vulnerable in each stock
  array<Type>   splitC_pgt(nP,nG,nT);       // Split value of mixed catch
  array<Type>   propMat_pgt(nP,nG,nT);      // Stock proportion mature estimated from vuln numbers
  array<Type>   propMat_gt(nG,nT);          // Mixed proportion mature estimated from vuln numbers
  vector<Type>  pEff_t(nT);                 // Mixed proportion effective


  // Management quantities
  vector<Type>  termDep_p(nP);     // Terminal depletion
  // Type          projSpawnBio;   // Projected spawning biomass
  // vector<Type>  projExpBio;     // Projected exploitable biomass



  // Observation model quantities
  array<Type>   zSum_pg(nP,nG);     // sum of log(It/Bt)
  array<int>    validObs_pg(nP,nG); // Number of observations
  array<Type>   z_pgt(nP,nG,nT);    // array to hold residuals for NLL calc
  array<Type>   zComb_pt(nP,nT);    // array to hold residuals for combined idx NLL calc
  array<Type>   SSR_pg(nP,nG);      // sum of squared resids
  vector<Type>  SSR_g(nG);          // sum of squared resids

  // Conditional MLEs for survey observations
  array<Type>   lnqhat_pg(nP,nG);    // log scale q_g
  array<Type>   qhat_pg(nP,nG);      // natural scale
  array<Type>   qComb_pg(nP,nG);     // natural scale catchability for combined index components
  array<Type>   qComb_pt(nP,nT);     // natural scale convex combination of survey qs
  array<Type>   tauComb_pt(nP,nT);   // natural scale convex combination of survey qs
  array<Type>   tauComb_pg(nP,nG);   // natural scale obs error SD for combined index

  // Count model years for each pop, based on tInitModel
  vector<int> nModelYrs_p(nP);
  nModelYrs_p.fill(nT);
  nModelYrs_p -= tInitModel_p;

  // Initialise all arrays at 0
  N_apt.fill(0.0);
  movN_apt.fill(0.0);
  endN_apt.fill(0.0);
  spawnN_apt.fill(0.0);
  tmpN_apt.fill(0.0);
  initN_apt.fill(0.0);
  termN_apt.fill(0.0);
  B_apt.fill(0.0);
  B_pt.fill(0.0);
  vulnB_apgt.fill(0.0);
  vulnB_pgt.fill(0.0);
  vulnN_apgt.fill(0.0);
  vulnN_pgt.fill(0.0);
  SB_pt.fill(0.0);
  bhR_pt.fill(0.0);
  R_pt.fill(0.0);
  M_apt.fill(0.0);
  Z_apt.fill(0.0);
  appZ_apt.fill(0.0);
  U_pgt.fill(0.0);
  U_apgt.fill(0.0);
  F_pgt.fill(0.0);
  qF_pgt.fill(0.0);
  lnqF_pgt.fill(0.0);
  totSurv_apt.fill(0.0);
  uAge_apgt.fill(0.0);
  catAge_apgt.fill(0.0);
  predPA_apgt.fill(-1);
  totC_pgt.setZero();
  pondC_apgt.setZero();
  pondC_pgt.setZero();
  pondC_gt.setZero();
  combIhat_gt.setZero();
  psi_pgt.setZero();
  psi_gt.setZero();
  propMat_pgt.setZero();
  propMat_gt.setZero();
  pEff_t.fill(0);
  movFracM_t.setZero();
  qComb_pt.fill(1.0);
  qComb_pg.setZero();
  tauComb_pg.setZero();
  propCatAgeDiluted_apgt.setZero();
  totCatAge_apgt.setZero();
  // vulnerable state variables when diluted with external fish
  vulnBdiluted_apgt.setZero();
  vulnBdiluted_pgt.setZero();
  vulnNdiluted_apgt.setZero();
  vulnNdiluted_pgt.setZero();
  totCundiluted_pgt.setZero();
  totCdiluted_pgt.setZero();
  uAgeDiluted_apgt.setZero();
  catAgeUndiluted_apgt.setZero();
  catAgeDiluted_apgt.setZero();
  propCatAgeDiluted_apgt.setZero();
  totCpredilute_pgt.setZero();
  SSR_pg.setZero();
  SSR_g.setZero();

  // Now fill in qComb_pg
  qComb_pg = exp(lnqComb_pg);
  tauComb_pg = exp(lntauObsComb_pg);


  /*\/\/\/\/\ PROCEDURE SECTION /\/\/\/\*/
  // Calculate life schedules - growth, maturity
  // Fill ages vector
  for( int a = 0; a < nA; a++)
  {
    age_a(a) = a+1;
    // if( a < juveMage )
    //   M_apt.transpose().col(a) += Mjuve;
  }
  // Calculate length curve
  len_a = Linf + (L1-Linf)*exp(-vonK*(age_a-1.));
  len_a(0) = inputL1;
  // Weight
  wt_a  = lenWt(0)*pow(len_a,lenWt(1));  

  // Compute mean and projection weight at age
  array<Type> meanWt_agp(nA,nG,nP);
  array<Type> meanWt_ap(nA,nP);
  array<Type> projWt_agp(nA,nG,nP);
  array<Type> projWt_ap(nA,nP);
  meanWt_agp.setZero();
  projWt_agp.setZero();
  meanWt_ap.setZero();
  projWt_ap.setZero();

  // Loop and calculate mean and projection
  // weights by gear and for the stock average
  for( int p = 0; p < nP; p++)
      for( int a = 0; a < nA; a++ )
      {
        for( int t = tInitModel_p(p); t < nT; t++)
        {
          // Gear
          for( int g = 0; g < nG; g++ )
            meanWt_agp(a,g,p) += W_apgt(a,p,g,t)/nT;

          // Stock average
          meanWt_ap(a,p) += W_apt(a,p,t)/nT;          
        }

        for( int t = nT - nYearsProjAve; t < nT; t++ )
        {
          // Gear average
          for( int g = 0; g < nG; g++ )
            projWt_agp(a,g,p) += W_apgt(a,p,g,t) / nYearsProjAve;

          // Stock average
          projWt_ap(a,p) += W_apt(a,p,t) / nYearsProjAve;
        }
      }

  // Calculate Mortality time series
  // First year uses the initial M
  // Calculate meanM for ssbpr calc
  vector<Type> meanM_p(nP);
  meanM_p.setZero();
  if( densityDepM == 0 )
  {
    for( int p = 0; p < nP; p++ )
    {
      M_apt.col(0).col(p).segment(juveMage,nA - juveMage) += M_p(p);
      for( int t = 1; t < nT; t++ )
      {
        M_apt.col(t).col(p).segment(juveMage,nA - juveMage) = 
            M_apt.col(t-1).col(p).segment(juveMage,nA - juveMage) *
            exp(sigmaM * omegaM_pt(p,t));
      }
    }
  }

  // This is a trick, as the final year isn't included yet, so 
  // we can just sum and divide by nT
  // Overwrite Mjuve_p with eqbm M if DDM is being used
  if(densityDepM == 1)
  {
    M0_p = M_p + exp( - m1_p );
    if( juveMsource == 0 )
      Mjuve_p = M0_p;   

  }

  for( int p = 0; p < nP; p++ )
  {
    if( densityDepM == 0 )
    {
      for( int t = tInitModel_p(p); t < nT; t++)
        meanM_p(p) += M_apt(juveMage,p,t)/nModelYrs_p(p);
      // Now compute the projected year's M
      for( int t = nT - nYearsProjAve; t < nT; t++)
        M_apt.col(nT).col(p).segment(juveMage,nA - juveMage) += M_apt(juveMage,p,t)/nYearsProjAve;


      M0_p(p) = meanM_p(p);
    }

    for( int t = tInitModel_p(p); t <= nT; t++  )
    {
      if( densityDepM == 0 )
      {
        if( juveMsource == 0 )
          M_apt.col(t).col(p).segment(0,juveMage) += meanM_p(p);

        if( juveMsource == 1 )
          M_apt.col(t).col(p).segment(0,juveMage) += Mjuve_p(p);
      }
      if( densityDepM == 1 )
        M_apt.col(t).col(p).segment(0,juveMage) += Mjuve_p(p);        
    }

    // Save meanM as Mjuve for plotting later
    if( juveMsource == 0 & densityDepM == 0 )
      Mjuve_p(p) += meanM_p(p);

  }

  


  // Calculate selectivity
  for( int g = 0; g < nG; g++ )
  {
    Type selX = 0.;
    Type maxSel = 1e-6;
    // Check seltype switch (asymptotic vs dome)
    if( selType_g(g) == 0)
    {
      // asymptotic
      Type xSel50   = SelAlpha_g(g);
      Type xSel95   = SelAlpha_g(g) + SelBeta_g(g);
      Type tmp      = log(19.)/( xSel95 - xSel50 );
      for( int a = 0; a < nA; a++ )
      {
        if( selX_g(g) == 1 )
          selX = len_a(a);
        if( selX_g(g) == 0)
          selX = age_a(a);

        sel_ag(a,g) = 1./( 1. + exp(-tmp*( selX - xSel50 ) ) );

        if( sel_ag(a,g) > maxSel )
          maxSel = sel_ag(a,g);
      }
    }

    if( selType_g(g) == 1)
    {
      // domed Normal selectivity
      Type nSelMean = SelAlpha_g(g);
      Type nSelSD   = SelBeta_g(g);
      for( int a = 0; a < nA; a++)
      {
        // Check if length or age based selectivity
        if( selX_g(g) == 1 )
          selX = len_a(a);
        if( selX_g(g) == 0)
          selX = age_a(a);
        // Compute selectivity function
        sel_ag(a,g) = exp(-1. * square((selX - nSelMean)/nSelSD) );

        if( sel_ag(a,g) > maxSel )
          maxSel = sel_ag(a,g);
      }
    }

    if( selType_g(g) == 2)
    {
      // asymptotic
      Type xSel95   = SelAlpha_g(g);
      Type xSel50   = SelAlpha_g(g) + SelBeta_g(g);
      Type tmp      = log(19.)/( xSel95 - xSel50 );
      for( int a = 0; a < nA; a++ )
      {
        if( selX_g(g) == 1 )
          selX = len_a(a);
        if( selX_g(g) == 0)
          selX = age_a(a);

        sel_ag(a,g) = 1./( 1. + exp(-tmp*( selX - xSel50 ) ) );

        if( sel_ag(a,g) > maxSel )
          maxSel = sel_ag(a,g);
      }
    }

    if( selType_g(g) == 3) // asymptotic domed
    {
      // Get scalar for this type
      // asymptotic
      Type uSel50   = SelAlpha_g(g);
      Type uSel95   = (SelAlpha_g(g) + SelBeta_g(g));

      Type dSel95   = dSelBeta_g(g);
      Type dSel50   = (dSelAlpha_g(g) + dSelBeta_g(g));

      Type utmp      = log(19.)/( uSel95 - uSel50 );
      Type dtmp      = log(19.)/( dSel95 - dSel50 );

      for( int a = 0; a < nA; a++ )
      {
        if( selX_g(g) == 1 )
          selX = len_a(a);
        if( selX_g(g) == 0)
          selX = age_a(a);

        Type tmpSelu = 1./( 1. + exp(-utmp*( selX - uSel50 ) ) );
        Type tmpSeld = 1./( 1. + exp(-dtmp*( selX - dSel50 ) ) );


        sel_ag(a,g) = tmpSelu * tmpSeld;

        if( sel_ag(a,g) > maxSel)
          maxSel = sel_ag(a,g);
      }
    }

    sel_ag.col(g) /= maxSel;

    // Force zero selectivity at age-1 (change to less than minAge)
    if(minAge_g(g) > 1 & selType_g(g) != 2)
      sel_ag.col(g).segment(0,minAge_g(g)-1) = 0;

    // Now do time-varying selectivity by stock
    for( int p = 0; p < nP; p++)
      for( int t = 0; t < nT; t++)
      {
        Type selX = 0.;
        Type maxSel = 0.;
        // Check seltype switch (asymptotic vs dome)
        if( selType_g(g) == 0)
        {
          // asymptotic
          Type xSel50   = SelAlpha_pgt(p,g,t);
          Type xSel95   = SelAlpha_pgt(p,g,t) + SelBeta_pgt(p,g,t);
          Type tmp      = log(19.)/( xSel95 - xSel50 );
          for( int a = 0; a < nA; a++ )
          {
            if( selX_g(g) == 1 )
              selX = len_a(a);
            if( selX_g(g) == 0)
              selX = age_a(a);

            sel_apgt(a,p,g,t) = 1./( 1. + exp(-tmp*( selX - xSel50 ) ) );

            if( sel_apgt(a,p,g,t) > maxSel )
              maxSel = sel_apgt(a,p,g,t);
          }
        }

        // Check seltype switch (asymptotic vs dome)
        if( selType_g(g) == 2)
        {
          // asymptotic
          Type xSel95   = SelAlpha_pgt(p,g,t);
          Type xSel50   = SelAlpha_pgt(p,g,t) + SelBeta_pgt(p,g,t);
          Type tmp      = log(19.)/( xSel95 - xSel50 );
          for( int a = 0; a < nA; a++ )
          {
            if( selX_g(g) == 1 )
              selX = len_a(a);
            if( selX_g(g) == 0)
              selX = age_a(a);

            sel_apgt(a,p,g,t) = 1./( 1. + exp(-tmp*( selX - xSel50 ) ) );

            if( sel_apgt(a,p,g,t) > maxSel )
              maxSel = sel_apgt(a,p,g,t);
          }
        }

        if( selType_g(g) == 1)
        {
          // domed Normal selectivity
          Type nSelMean = SelAlpha_pgt(p,g,t);
          Type nSelSD   = SelBeta_pgt(p,g,t);
          for( int a = 0; a < nA; a++)
          {
            // Check if length or age based selectivity
            if( selX_g(g) == 1 )
              selX = len_a(a);
            if( selX_g(g) == 0)
              selX = age_a(a);
            // Compute selectivity function
            sel_apgt(a,p,g,t) = exp(-1. * square((selX - nSelMean)/nSelSD) );

            if( sel_apgt(a,p,g,t) > maxSel )
              maxSel = sel_apgt(a,p,g,t);
          }

        }  

        if( selType_g(g) == 3) // asymptotic domed
        {
          // asymptotic
          Type uSel50   = SelAlpha_pgt(p,g,t);
          Type uSel95   = SelAlpha_pgt(p,g,t) + SelBeta_pgt(p,g,t);

          Type dSel95   = dSelBeta_pgt(p,g,t);
          Type dSel50   = dSelAlpha_pgt(p,g,t) + dSelBeta_pgt(p,g,t);

          Type utmp     = log(19.)/( uSel95 - uSel50 );
          Type dtmp     = log(19.)/( dSel95 - dSel50 );
          for( int a = 0; a < nA; a++ )
          {
            if( selX_g(g) == 1 )
              selX = len_a(a);
            if( selX_g(g) == 0)
              selX = age_a(a);

            Type tmpSelu = 1./( 1. + exp(-utmp*( selX - uSel50 ) ) );
            Type tmpSeld = 1./( 1. + exp(-dtmp*( selX - dSel50 ) ) );


            sel_apgt(a,p,g,t) = tmpSelu * tmpSeld;

            if( sel_apgt(a,p,g,t) > maxSel)
              maxSel = sel_apgt(a,p,g,t);
          }
        }

        sel_apgt.col(t).col(g).col(p) /= maxSel;


        // Force zero sel at age 1 (change to < minAge)
        if(minAge_g(g) > 1 & selType_g(g) != 2)
          sel_apgt.col(t).col(g).col(p).segment(0,minAge_g(g)-1) = 0;
      }
  }

  // We're gonna update this to an array later
  vector<int> chronIdx_g(nG);
  chronIdx_g.setZero();
  

  // Calculate SSBpr, recruitment parameters
  // Calculate equilibrium unfished survivorship.
  surv_ap.fill(1.0);
  initSurv_ap.fill(1.0);
  array<Type> initZ_ap(nA,nP);
  initZ_ap.setZero();
  for( int p = 0; p < nP; p++ )
  {
    initZ_ap.col(p).segment(0,juveMage) += Mjuve_p(p);
    initZ_ap.col(p).segment(juveMage,nA - juveMage) += M0_p(p);

    if( initFcode_p(p) > 0 )
      initZ_ap.col(p) += sel_apgt.col(tInitModel_p(p)).col(initFcode_p(p) - 1).col(p) * Finit_p(p);
    
    for(int a = 1; a < nA; a++ ) 
    {
      if( a < juveMage )
        surv_ap(a,p) = surv_ap(a-1,p)*exp(-Mjuve_p(p));
      if( a >= juveMage)
        surv_ap(a,p) = surv_ap(a-1,p)*exp(-M0_p(p));

      initSurv_ap(a,p) = initSurv_ap(a-1,p)*exp(-initZ_ap(a-1,p));
      if( initMethod == "surv" )
        initSurv_ap(a,p) *= initN_mult_ap(a-1,p);
    }

    surv_ap(nA-1,p) /= (1.0 - exp(-M0_p(p)));  
    initSurv_ap(nA-1,p) /= (1.0 - exp(-initZ_ap(nA-1,p)));  
  }

  // SSBpr
  vector<Type> totbpr_p(nP);
  totbpr_p.setZero();
  for( int p = 0; p < nP; p++)
    for( int a = 0; a < nA; a ++)
    {
      if( a < juveMage)
        phi_p(p) += ( mat_a(a) * meanWt_ap(a,p) * surv_ap(a,p) * exp(-spawnTiming*Mjuve_p(p) ) );

      if( a >= juveMage )
      {
        phi_p(p) += ( mat_a(a) * meanWt_ap(a,p) * surv_ap(a,p) * exp(-spawnTiming*M0_p(p) ) );
        totbpr_p(p) += meanWt_ap(a,p) * surv_ap(a,p);
      }

    }


  // Now compute R0 from B0
  R0_p = B0_p/phi_p;

  for( int p = 0; p < nP; p++ )
    for( int a = juveMage; a < nA; a++ )
      totB0_p(p) += meanWt_ap(a,p) * surv_ap(a,p) * R0_p(p);

  // Need to calculate some helper indices
  int nDiscFleets = 0;
  int nCtsFleets = 0;
  
  for( int g = 0; g < nG; g++)
  {
    if(mortType_g(g)==0)
      nDiscFleets += 1;

    if(mortType_g(g)==1)
      nCtsFleets += 1;
  }
  // Now a vector of indices of each
  vector<int> ctsFleetIdx(nCtsFleets);
  vector<int> discFleetIdx(nDiscFleets);
  int c = 0;
  int d = 0;
  for( int g = 0; g < nG; g++)
  {
    if(mortType_g(g)==0)
    {
      discFleetIdx(d) = g;
      d++;
    }

    if(mortType_g(g)==1)
    {
      ctsFleetIdx(c) = g;
      c++;
    }
  }


  // Beverton-Holt a parameter.
  Type eggBscalar = 1;
  if( SRindVar == 2)
    eggBscalar = fec * 1e3 * 0.5;
  reca_p = 4.*rSteepness_p*R0_p/(B0_p*(1.-rSteepness_p) * eggBscalar);
  // Beverton-Holt b parameter.
  recb_p = (5.*rSteepness_p-1.)/(B0_p*(1.-rSteepness_p) * eggBscalar);
  
  int propEffVecIdx = 0;
  
  // Loop over time steps, run pop dynamics
  for( int t = 0; t < nT; t++ )
  {
    // Check for calculating proportion effective
    if( calcPropEff_t(t) == 1 )
    {
      pEff_t(t) = propEffBounds(0) + 
                  ( propEffBounds(1) - propEffBounds(0) ) /
                  ( 1 + exp( -logitPropEff_vec( propEffVecIdx ) ) );
      propEffVecIdx++;
    }

    // Initialise
    // Block to initialise population dynamics
    for( int p = 0; p < nP; p++ )
    {
      if( t == tInitModel_p(p))
      {
        if( initRcode_p(p) == 0)
          Rinit_p(p) = R0_p(p);
        if( initRcode_p(p) == 1)
          Rinit_p(p) = exp(lnRinit_p(p));

        N_apt.col(t).col(p) =  Rinit_p(p) * initSurv_ap.col(p);

        // Initialise qF_pgt
        for( int g = 0; g < nG; g++ )
        {
          if( mortType_g(g) == 1 | catSeriesType_g(g) == 1 )
          {
            lnqF_pgt(p,g,t) += lnqFinit_pg(p,g);
          }
        }
        

        if( initMethod == "nums" )
          N_apt.col(t).col(p) *= initN_mult_ap.col(p);

        // Calc total biomass and rec in first year
        B_apt.col(t).col(p) = N_apt.col(t).col(p) * W_apt.col(t).col(p);
        B_pt(p,t) = B_apt.col(t).col(p).segment(juveMage,nA-juveMage).sum();

        if( densityDepM == 1 )
          M_apt.col(t).col(p).segment(juveMage,nA - juveMage) = ( M_p(p) + exp( -m1_p(p) * B_pt(p,t)/totB0_p(p) ) ) * exp( sigmaM * omegaM_pt(p,t));

        R_pt(p,t)  = N_apt(0,p,t);
      }

      // Otherwise extend the RW
      if(t > tInitModel_p(p))
      {
        for( int g = 0; g < nG; g++)
          if(mortType_g(g) == 1 | catSeriesType_g(g) == 1 )
            lnqF_pgt(p,g,t) = lnqF_pgt(p,g,t-1) + sddeltaqF_g(g) * deltaqF_pgt(p,g,t);
      }


      // Calculate cts fishing mortality rate
      if( t >= tInitModel_p(p) )
      {
        Z_apt.col(t).col(p) += M_apt.col(t).col(p);
        for( int g = 0; g < nG; g++)
          if(mortType_g(g) == 1 | catSeriesType_g(g) == 1 )
          {
            qF_pgt(p,g,t)         = exp(lnqF_pgt(p,g,t));
            if(mortType_g(g) == 1)
            {
              if(C_pgt(p,g,t) > 0)
                F_pgt(p,g,t)          = qF_pgt(p,g,t) * E_pgt(p,g,t);

              Z_apt.col(t).col(p)  += sel_apgt.col(t).col(g).col(p) * F_pgt(p,g,t);
            }
          }
      }


    }

    // Now loop over fleets and take
    // catch as necessary
    Type prevTime = 0.;
    // First, get this year's beginning of year numbers
    array<Type>  tmpN_ap(nA,nP);
    array<Type>  wtAge_ap(nA,nP);
    array<Type>  wtAge_apg(nA,nP,nG);
    tmpN_ap.fill(0);
    wtAge_ap.fill(0);
    wtAge_apg.fill(0);

    tmpN_ap   = N_apt.col(t);
    wtAge_ap  = W_apt.col(t);
    wtAge_apg = W_apgt.col(t);

    // Save tmpN_ap
    tmpN_apt.col(t) = tmpN_ap;

    // Apply movement effectively between time-steps 
    // after recruitment - might need to adjust later for
    // recruitment to be a separate process 
    if( useMovement == 1 )
    {
      array<Type> tmpMovN_ap(nA,nP);
      array<Type> tmpInitN_ap(nA,nP);
      array<Type> tmpTermN_ap(nA,nP);
      tmpMovN_ap.setZero();
      tmpInitN_ap.setZero();
      tmpTermN_ap.setZero();
      array<Type> tmpM_ap(nA,nP);
      tmpM_ap = M_apt.col(t);
      tmpN_apt.col(t) = tmpN_ap;
      tmpMovN_ap = applyMovement( tmpM_ap,
                                  tmpN_ap,
                                  mov_ppa,
                                  tmpInitN_ap,
                                  tmpTermN_ap,
                                  prevTime,
                                  Type(0) );

      initN_apt.col(t) = tmpInitN_ap;
      termN_apt.col(t) = tmpTermN_ap;
      movN_apt.col(t) = tmpMovN_ap;
      tmpN_ap = tmpMovN_ap;
    }

    // here is where we would re-calculate the
    // fleet order
    calcFleetOrder( fleetTiming_g,
                    chronIdx_g,
                    mortType_g);
    
    for( int cIdx = 0; cIdx < nDiscFleets; cIdx ++ )
    {
      // Get actual fleet number
      int gIdx = discFleetIdx(cIdx);

      if(mortType_g(gIdx) == 0)
      {
        // Loop over stocks to compute vuln numbers for this fleet/time
        for( int p = 0; p < nP; p++ )
        {  
          // First, vulnerable numbers is found by reducing
          // numbers by fracM
          tmpN_ap.col(p) = N_apt.col(t).col(p) * exp( - fleetTiming_g(gIdx) * Z_apt.col(t).col(p) );

          if( cIdx > 0)
            for( int a = 0; a < nA; a++)
            {
              for( int cc = 0; cc < cIdx; cc++ )  
              {
                int gg = chronIdx_g(cc);
                // Deplete tmpN
                tmpN_ap(a,p) *= (1 - U_pgt(p,gg,t) * sel_apgt(a,p,gg,t));
                
              }
              // apply posFun
              tmpN_ap(a,p) = posfun(tmpN_ap(a,p),Type(1e-6), posPen);
            }

          // refactored vulnerability calcs to remove a loop
          vulnN_apgt.col(t).col(gIdx).col(p) = tmpN_ap.col(p) * sel_apgt.col(t).col(gIdx).col(p);
          vulnB_apgt.col(t).col(gIdx).col(p) = vulnN_apgt.col(t).col(gIdx).col(p)*wtAge_apg.col(gIdx).col(p);  

          // Now add dilution
          vulnNdiluted_apgt.col(t).col(gIdx).col(p) += vulnN_apgt.col(t).col(gIdx).col(p);
          vulnNdiluted_apgt.col(t).col(gIdx).col(p).segment(minAgeDilute,nA-minAgeDilute) += alphaDilute*(diluteN_apt.col(t).col(p) * sel_apgt.col(t).col(gIdx).col(p)).segment(minAgeDilute,nA-minAgeDilute);
          // Consider adding SOG weight-at-age for dilution as well!!
          vulnBdiluted_apgt.col(t).col(gIdx).col(p) = vulnNdiluted_apgt.col(t).col(gIdx).col(p)*wtAge_apg.col(gIdx).col(p);

          // Now sum for annual totals
          vulnB_pgt(p,gIdx,t)         = vulnB_apgt.col(t).col(gIdx).col(p).sum();
          vulnN_pgt(p,gIdx,t)         = vulnN_apgt.col(t).col(gIdx).col(p).sum();
          vulnNdiluted_pgt(p,gIdx,t)  = vulnNdiluted_apgt.col(t).col(gIdx).col(p).sum();
          vulnBdiluted_pgt(p,gIdx,t)  = vulnBdiluted_apgt.col(t).col(gIdx).col(p).sum();

          // Mixed population vulnerable states
          mixedVulnN_agt.col(t).col(gIdx) += vulnN_apgt.col(t).col(gIdx).col(p);
          mixedVulnB_agt.col(t).col(gIdx) += vulnB_apgt.col(t).col(gIdx).col(p);
          mixedVulnN_gt(gIdx,t) += vulnN_pgt(p,gIdx,t);
          mixedVulnB_gt(gIdx,t) += vulnB_pgt(p,gIdx,t);
        }

        // Add fleet catch to totC if a commercial
        // fishery
        if( fleetType_g(gIdx) == 1 )
          totC_pgt.col(t).col(gIdx) += C_pgt.col(t).col(gIdx);

        if( catSeriesType_g(gIdx) == 1)
          totC_pgt.col(t).col(gIdx) *= qF_pgt.col(t).col(gIdx);

        if( fleetType_g(gIdx) == 2 | fleetType_g(gIdx) == 3 )
        {
          // Apply SOK conversion here to substock specific
          // catch
          for( int p = 0; p < nP; p++)
          {
            if( C_pgt(p,gIdx,t) > 0 )
            {
              // Calculate biomass of ponded fish, take
              // as removals

              Type tmpPsi = 0;
              Type tmpEff = 0;
              Type tmpMat = 0;
              Type tmpCat = C_pgt(p,gIdx,t);
              
              if( catSeriesType_g(gIdx) == 1)
                tmpCat *= qF_pgt(p,gIdx,t);

              tmpEff = pEff_t(t);

              vector<Type> vulnN_a(nA);
              vector<Type> sel_a(nA);

              vulnN_a = vulnN_apgt.col(t).col(gIdx).col(p);
              sel_a   = sel_apgt.col(t).col(gIdx).col(p);

              Type tmpPondC = 0;

              
              tmpPondC = convertSOK(  tmpCat,
                                      mat_a, 
                                      sel_a,
                                      tmpEff,
                                      pFem,
                                      vulnN_a,
                                      tmpMat,
                                      fec,
                                      wtAge_apg.col(gIdx).col(p),
                                      gamma_g(gIdx),
                                      tmpPsi,
                                      sokInitF );

              pondC_pgt(p,gIdx,t) = tmpPondC;

              psi_pgt(p,gIdx,t) = tmpPsi;
              propMat_pgt(p,gIdx,t) = tmpMat;
            }

            totC_pgt(p,gIdx,t) += pondC_pgt(p,gIdx,t);
          }


        }


        // Now we need to split the mixed catch if necessary
        if( mC_gt(gIdx,t) > 0 )
        {
          Type mixedCat = 0;
          // Compute mixed catch
          if( fleetType_g(gIdx) == 1 )
             mixedCat = mC_gt(gIdx,t);

          if( fleetType_g(gIdx) == 2 | fleetType_g(gIdx) == 3 )
          {

            Type tmpPsi = 0;
            Type tmpEff = pEff_t(t);
            Type tmpMat = 0;
            Type tmpCat = mC_gt(gIdx,t);

            if( catSeriesType_g(gIdx) == 1)
              tmpCat *= qF_pgt(0,gIdx,t);


            vector<Type> vulnN_a(nA);
            vector<Type> sel_a(nA);

            vulnN_a = mixedVulnN_agt.col(t).col(gIdx);
            sel_a   = sel_apgt.col(t).col(gIdx).col(0);

            // Apply SOK conversion if mC_gt > 0
            Type tmpPond = 0;
            tmpPond = convertSOK( tmpCat,
                                  mat_a, 
                                  sel_a,
                                  tmpEff,
                                  pFem,
                                  vulnN_a,
                                  tmpMat,
                                  fec,
                                  wtAge_ap.col(0),
                                  gamma_g(gIdx),
                                  tmpPsi,
                                  sokInitF );

            pondC_gt(gIdx,t) = tmpPond;

            psi_gt(gIdx,t) = tmpPsi;
            propMat_gt(gIdx,t) = tmpMat;

            mixedCat = pondC_gt(gIdx,t);
          }

          // Now split by vulnerable biomass
          vector<Type> vulnBioRatio = vulnB_pgt.col(t).col(gIdx);
          vulnBioRatio /= mixedVulnB_gt(gIdx,t);

          // Save split catch
          splitC_pgt.col(t).col(gIdx) = vulnBioRatio * mixedCat;

          // now add the mixed catch to the total removals
          totC_pgt.col(t).col(gIdx) += splitC_pgt.col(t).col(gIdx);

          // Save ponded fish if SOK fleet
          if( fleetType_g(gIdx) == 2 | fleetType_g(gIdx) == 3 )
            pondC_pgt.col(t).col(gIdx) += vulnBioRatio * mixedCat;

        }

        // Now loop over stocks again
        for( int p = 0; p < nP; p++)
        {
          totCpredilute_pgt(p,gIdx,t) = totC_pgt(p,gIdx,t);
          if( t >= tInitModel_p(p) )
          {

            // Add a conditional to check if betaDilute_g(gIdx) > 0

            // Get numbers caught at age
            totCundiluted_pgt(p,gIdx,t) = totCpredilute_pgt(p,gIdx,t) * (1 - betaDilute_g(gIdx));
            totCdiluted_pgt(p,gIdx,t)   = totCpredilute_pgt(p,gIdx,t) * ( betaDilute_g(gIdx));

            // Caclulate proportion-at-age in each fleet's 
            // vuln biomass to convert catch to numbers
            uAge_apgt.col(t).col(gIdx).col(p) = vulnB_apgt.col(t).col(gIdx).col(p) / vulnB_pgt(p,gIdx,t);
            // And for the diluted stock
            if( vulnBdiluted_pgt(p,gIdx,t) > 0)
              uAgeDiluted_apgt.col(t).col(gIdx).col(p) = vulnBdiluted_apgt.col(t).col(gIdx).col(p) / (vulnBdiluted_pgt(p,gIdx,t) );

            // catch-at-age for diluted and undiluted
            catAgeUndiluted_apgt.col(t).col(gIdx).col(p)  = uAge_apgt.col(t).col(gIdx).col(p) * totCundiluted_pgt(p,gIdx,t) / wtAge_apg.col(gIdx).col(p);
            catAgeDiluted_apgt.col(t).col(gIdx).col(p)    = uAgeDiluted_apgt.col(t).col(gIdx).col(p) * totCdiluted_pgt(p,gIdx,t) / wtAge_apg.col(gIdx).col(p);

            // Calculate the proportion of numbers-at-age that is from the main stock after dilution
            for( int a = 0; a < nA; a++)
              if(sel_apgt(a,p,gIdx,t) > 0 & vulnNdiluted_apgt(a,p,gIdx,t) > 0 )
                propCatAgeDiluted_apgt(a,p,gIdx,t) = vulnN_apgt(a,p,gIdx,t)/vulnNdiluted_apgt(a,p,gIdx,t);

            catAge_apgt.col(t).col(gIdx).col(p) =  catAgeUndiluted_apgt.col(t).col(gIdx).col(p) + propCatAgeDiluted_apgt.col(t).col(gIdx).col(p) * catAgeDiluted_apgt.col(t).col(gIdx).col(p);

            // Save ponded fish at age if SOK fleet
            if( fleetType_g(gIdx) == 2 | fleetType_g(gIdx) == 3 )
              pondC_apgt.col(t).col(gIdx).col(p) += catAge_apgt.col(t).col(gIdx).col(p);

            // Update totC for diluted removals and calculating harvest rates/Fs
            totC_pgt(p,gIdx,t) = (catAge_apgt.col(t).col(gIdx).col(p) * wtAge_apg.col(gIdx).col(p)).sum();
            
            // Calculate F - there might be a better way here...
            // Read MacCall's paper on alternatives to Pope's approx
            // Question is: what do we use for estimating F? Just U? Crank this
            // on paper first.
            // Pope's approx works better within a single age class, or
            // with DD models (F = log(Nt/Nt-1)-M)
            if(fleetType_g(gIdx) == 1)
              U_pgt(p,gIdx,t) = totC_pgt(p,gIdx,t) / vulnB_pgt(p,gIdx,t);
            if( fleetType_g(gIdx) > 1)
              U_pgt(p,gIdx,t) = pondC_pgt(p,gIdx,t) / vulnB_pgt(p,gIdx,t);
            

            // need to revisit what to do here... U is an important
            // variable now that drives the intra-annual population
            // dynamics, so we should add another variable to keep
            // trace of dead ponded fish (essentially dead discards...)

            // if( fleetType_g(gIdx) == 2 | fleetType_g(gIdx) == 3 )
            //   U_pgt(p,gIdx,t) *= (1 - exp(-postPondM_g(gIdx)));

            // Add a penalty if U is greater than .9 or smaller than 1e-9
            if( C_pgt(p,gIdx,t) > 0 | mC_gt(gIdx,t) > 0 )
            {
              U_pgt(p,gIdx,t) = 1 - posfun( 1. - U_pgt(p,gIdx,t), Type(.05), posPen );
              U_pgt(p,gIdx,t) = posfun( U_pgt(p,gIdx,t), Type(1e-9), posPen );
            }
            
            U_apgt.col(t).col(gIdx).col(p) = U_pgt(p,gIdx,t) * sel_apgt.col(t).col(gIdx).col(p);
            
          }

        }
      } // END if mortType_g == 0
    } // END loop over chronIdx


    // Now loop over continuous fleets, 
    for( int cIdx = 0; cIdx < nCtsFleets; cIdx ++ )
    {
      // estimate vuln bio (start of year)
      for( int p = 0; p < nP; p++ )
      {
        int gIdx = ctsFleetIdx(cIdx);
        // refactored vulnerability calcs to remove a loop
        vulnN_apgt.col(t).col(gIdx).col(p) = N_apt.col(t).col(p) * sel_apgt.col(t).col(gIdx).col(p);
        vulnB_apgt.col(t).col(gIdx).col(p) = vulnN_apgt.col(t).col(gIdx).col(p)*wtAge_apg.col(gIdx).col(p);  

        // Now add dilution
        vulnNdiluted_apgt.col(t).col(gIdx).col(p) += vulnN_apgt.col(t).col(gIdx).col(p);
        vulnNdiluted_apgt.col(t).col(gIdx).col(p).segment(minAgeDilute,nA-minAgeDilute) += alphaDilute*(diluteN_apt.col(t).col(p) * sel_apgt.col(t).col(gIdx).col(p)).segment(minAgeDilute,nA-minAgeDilute);
        // Consider adding SOG weight-at-age for dilution as well!!
        vulnBdiluted_apgt.col(t).col(gIdx).col(p) = vulnNdiluted_apgt.col(t).col(gIdx).col(p)*wtAge_apg.col(gIdx).col(p);

        // Now sum for annual totals
        vulnB_pgt(p,gIdx,t)         = vulnB_apgt.col(t).col(gIdx).col(p).sum();
        vulnN_pgt(p,gIdx,t)         = vulnN_apgt.col(t).col(gIdx).col(p).sum();
        vulnNdiluted_pgt(p,gIdx,t)  = vulnNdiluted_apgt.col(t).col(gIdx).col(p).sum();
        vulnBdiluted_pgt(p,gIdx,t)  = vulnBdiluted_apgt.col(t).col(gIdx).col(p).sum();

        // Mixed population vulnerable states
        mixedVulnN_agt.col(t).col(gIdx) += vulnN_apgt.col(t).col(gIdx).col(p);
        mixedVulnB_agt.col(t).col(gIdx) += vulnB_apgt.col(t).col(gIdx).col(p);
        mixedVulnN_gt(gIdx,t) += vulnN_pgt(p,gIdx,t);
        mixedVulnB_gt(gIdx,t) += vulnB_pgt(p,gIdx,t);


        // Calculate proportion biomass at age
        uAge_apgt.col(t).col(gIdx).col(p) = vulnB_apgt.col(t).col(gIdx).col(p) * exp(-fleetTiming_g(gIdx) *Z_apt.col(t).col(p)) / (vulnB_apgt.col(t).col(gIdx).col(p) * exp(-fleetTiming_g(gIdx) *Z_apt.col(t).col(p))).sum();

        // Calculate total observed catch
        totC_pgt(p,gIdx,t) = C_pgt(p,gIdx,t);

        // While we're here, split mixed catch in proportion to
        // vulnerable biomass
        // Now split by vulnerable biomass
        Type vulnBioRatio = vulnB_pgt(p,gIdx,t)/mixedVulnB_gt(gIdx,t);
      
        // Save split catch
        splitC_pgt(p,gIdx,t) = vulnBioRatio * mC_gt(gIdx,t);


        // now add the mixed catch to the total removals
        totC_pgt(p,gIdx,t) += splitC_pgt(p,gIdx,t);

        // Now dilute - not currently working as expected
        // // And for the diluted stock
        // if( vulnBdiluted_pgt(p,gIdx,t) > 0)
        //   uAgeDiluted_apgt.col(t).col(gIdx).col(p) = vulnBdiluted_apgt.col(t).col(gIdx).col(p) / (vulnBdiluted_pgt(p,gIdx,t) );

        // // 

        // // Get numbers caught at age
        // totCundiluted_pgt(p,gIdx,t) = totC_pgt(p,gIdx,t) * (1 - betaDilute_g(gIdx));
        // totCdiluted_pgt(p,gIdx,t)   = totC_pgt(p,gIdx,t) * ( betaDilute_g(gIdx));

        
        // // catch-at-age for diluted and undiluted
        // catAgeUndiluted_apgt.col(t).col(gIdx).col(p)  = uAge_apgt.col(t).col(gIdx).col(p) * totCundiluted_pgt(p,gIdx,t) / wtAge_apg.col(gIdx).col(p);
        // catAgeDiluted_apgt.col(t).col(gIdx).col(p)    = uAgeDiluted_apgt.col(t).col(gIdx).col(p) * totCdiluted_pgt(p,gIdx,t) / wtAge_apg.col(gIdx).col(p);

        // // Calculate the proportion of numbers-at-age that is from the main stock after dilution
        // for( int a = 0; a < nA; a++)
        //   if(sel_apgt(a,p,gIdx,t) > 0 & vulnNdiluted_apgt(a,p,gIdx,t) > 0 )
        //     propCatAgeDiluted_apgt(a,p,gIdx,t) = vulnN_apgt(a,p,gIdx,t)/vulnNdiluted_apgt(a,p,gIdx,t);

        // catAge_apgt.col(t).col(gIdx).col(p) =  catAgeUndiluted_apgt.col(t).col(gIdx).col(p) + propCatAgeDiluted_apgt.col(t).col(gIdx).col(p) * catAgeDiluted_apgt.col(t).col(gIdx).col(p);

        // // // Save ponded fish at age if SOK fleet
        // // if( fleetType_g(gIdx) == 2 | fleetType_g(gIdx) == 3 )
        // //   pondC_apgt.col(t).col(gIdx).col(p) += catAge_apgt.col(t).col(gIdx).col(p);

        // // Update totC for diluted removals and calculating harvest rates/Fs
        // totC_pgt(p,gIdx,t) -= ((1 - propCatAgeDiluted_apgt.col(t).col(gIdx).col(p)) * catAgeDiluted_apgt.col(t).col(gIdx).col(p) * wtAge_apg.col(gIdx).col(p)).sum();
        
        // Loop over discrete fleets, remove discFleet catch and calculate
        // the expected cts fleet catch
        Type prevTime = 0;
        Type fracZ = 0;
        // Initialise temporary vulnerable biomass
        array<Type> tmpVulnB_apg(nA,nP,nG);
        tmpVulnB_apg = vulnB_apgt.col(t);
        for( int dIdx = 0; dIdx < nDiscFleets; dIdx ++ )
        {
          int ggIdx = discFleetIdx(dIdx);
          fracZ = fleetTiming_g(ggIdx) - prevTime;

          for( int a = 0; a < nA; a++ )
          {
            // Add catch at each age
            expC_pgt(p,gIdx,t) += tmpVulnB_apg(a,p,gIdx) * (1 - exp(-1 * fracZ * Z_apt(a,p,t))) * F_pgt(p,gIdx,t)/Z_apt(a,p,t);

            // Now calculate remaining biomass
            tmpVulnB_apg(a,p,gIdx) = tmpVulnB_apg(a,p,gIdx) * exp(-fracZ * Z_apt(a,p,t)) * (1 - U_apgt(a,p,ggIdx,t)); 
          }

          prevTime = fleetTiming_g(ggIdx);
        }

        // Now calculate last piece of catch
        fracZ = 1-prevTime;
        for( int a = 0; a < nA; a++ )
          expC_pgt(p,gIdx,t) += tmpVulnB_apg(a,p,gIdx) * (1 - exp(-1 * fracZ * Z_apt(a,p,t))) * F_pgt(p,gIdx,t)/Z_apt(a,p,t);

        // totC_pgt(p,gIdx,t) = expC_pgt(p,gIdx,t);

      } // END p loop
      
    } // END cIdx loop (continuous fleets)

    // Finally, advance numbers at age
    // to the following time step by reducing
    // remaining numbers by the remaining mortality.

    // At the same time, calculate numbers at age at
    // spawn timing

    for( int p = 0; p < nP; p++ )
    {
      if( t >= tInitModel_p(p) )
      {
        for(int a = 0; a < nA; a++ )
        {
          // Deplete N_apt by all natural mortality
          endN_apt(a,p,t) = N_apt(a,p,t) * exp( - Z_apt(a,p,t));
          totSurv_apt(a,p,t) = exp(-M_apt(a,p,t));

          // Deplete by spawn timing fraction
          spawnN_apt(a,p,t) = N_apt(a,p,t) * exp( - spawnTiming * Z_apt(a,p,t));


          // Loop over gears and deplete by fishing
          for( int g = 0; g < nG; g++ )
          {
            endN_apt(a,p,t) *= (1 - U_pgt(p,g,t) * sel_apgt(a,p,g,t));
            totSurv_apt(a,p,t) *= (1 - U_pgt(p,g,t) * sel_apgt(a,p,g,t));

            if( fleetTiming_g(g) <= spawnTiming )
              spawnN_apt(a,p,t) *= (1 - U_pgt(p,g,t) * sel_apgt(a,p,g,t));
            
            // Add ponded fish back in to endN - but not spawnN
            if( fleetType_g(g) == 2 | fleetType_g(g) == 3)
            {
              Type fleetFrac = 1 - fleetTiming_g(g);
              endN_apt(a,p,t) += pondC_apgt(a,p,g,t) * exp( -postPondM_g(g) - fleetFrac * Z_apt(a,p,t));
            }

            endN_apt(a,p,t)   = posfun( endN_apt(a,p,t), Type(1e-6), posPen );
            spawnN_apt(a,p,t) = posfun( spawnN_apt(a,p,t), Type(1e-6), posPen );
          }

          // Advance age for a > 0 (recruits done below)
          if( a > 0 )
            N_apt(a,p,t+1) = endN_apt(a-1,p,t);

          // And combine last two age classes in plus group
          if( a == nA - 1)
            N_apt(a,p,t+1) += endN_apt(a,p,t);

          // calculate total mortality
          appZ_apt(a,p,t) = -log(endN_apt(a,p,t)/N_apt(a,p,t));

        }

        // Calculate spawning biomass and expected BH recruitment
        calcSpawn(  spawnN_apt.col(t),
                    SB_pt,
                    M_apt.col(t),
                    W_apt.col(t),
                    mat_a,
                    bhR_pt,
                    Eggs_pt,
                    fec,
                    reca_p,
                    recb_p,
                    t,
                    tInitModel_p,
                    SRindVar );

        // Update numbers at age 1 for the following year
        // Average R recruitment
        if(avgRcode_p(p) == 1)
          N_apt(0,p,t+1) =  Rbar_p(p) * exp(sigmaR * omegaR_pt(p,t+1));

        // BH recruitment
        if( avgRcode_p(p) == 0 )
          N_apt(0,p,t+1) =  bhR_pt(p,t+1) * exp( sigmaR * omegaR_pt(p,t+1) );              

        // Record recruits
        R_pt(p,t+1) = N_apt(0,p,t+1);            

      }
      vector<Type> wtAge(nA);
      // Convert to biomass
      if( t < nT-1)
         wtAge = W_apt.col(t+1).col(p);
      else wtAge = projWt_ap.col(p);

      // Calculate start of year biomass
      B_apt.col(t+1).col(p) = N_apt.col(t+1).col(p) * wtAge;
      B_pt(p,t+1) = B_apt.col(t+1).col(p).segment(juveMage,nA - juveMage).sum();

      // Depensatory M if used
      if( densityDepM == 1 & t >= tInitModel_p(p) )
      {
        M_apt.col(t+1).col(p).segment(juveMage,nA-juveMage) = M_p(p) + exp( -m1_p(p) * (B_pt(p,t+1)/totB0_p(p)));

        if( t < nT - 1 )
          M_apt.col(t+1).col(p).segment(juveMage,nA-juveMage) *= exp( sigmaM * omegaM_pt(p,t+1) );
      }

    } // END p loop for advancing population dynamics
    
    // reset previous time to zero, bug keeps it at 1
    prevTime = 0;
  } // End t loop for pop dynamics - could potentially extend this loop to cover observation models 

  // No catch
  for( int p = 0; p < nP; p++)
    SB_pt(p,nT) = (B_apt.col(nT).col(p) * mat_a * exp(-spawnTiming*M_apt.col(nT).col(p)) ).sum();


  // Project spawning biomass to end of next time step under a range
  // of catch conditions
  array<Type> alloc_pg(nP,nG);
  for( int p = 0; p < nP; p++ )
  {
    vector<Type> totC_g(nG);
    totC_g.setZero();
    for( int t = catAllocYrs(0)-1; t < catAllocYrs(1); t++ )
      for( int g = 0; g < nG; g++ )
        totC_g(g) += totC_pgt(p,g,t);

    if( totC_g.sum() > 0 )
      for( int g = 0; g < nG; g++ )
        alloc_pg(p,g) = totC_g(g)/totC_g.sum(); 
  }

  // Make catch table for projections
  int nV = C_v.size();
  array<Type> catTable_pvk(nP,nV,4);

  array<Type> projN_ap(nA,nP);
  array<Type> projM_ap(nA,nP);
  projN_ap = N_apt.col(nT);
  projM_ap = M_apt.col(nT);

  projCatchBio( projN_ap,
                projM_ap,
                projWt_ap,
                projWt_agp,
                mat_a,
                fleetTiming_g,
                fleetType_g,
                chronIdx_g,
                spawnTiming,
                sel_apgt.col(nT-1),
                alloc_pg,
                C_v,
                postPondM_g,
                B0_p,
                catTable_pvk,
                nT);

  
  // Calculate likelihood functions
  // Observation models //
  
  // Initialise lnqhat and lntauObs at 0 
  lnqhat_pg.fill(0.0);
  z_pgt.fill(0.0);
  zComb_pt.fill(0.0);
  tauComb_pt.fill(0.0);
  zSum_pg.fill(0.0);
  validObs_pg.fill(0);
  etaSumSq_pg.setZero();
  nResids_pg.setZero();
  nObsAge_pg.setZero();
  tau2Age_pg.setZero();
  array<Type> tau2Obshat_pg(nP,nG);
  tau2Obshat_pg.setZero();
  // Loop over stocks and gear types
  for( int p = 0; p < nP; p++)
  {
    for( int gIdx = 0; gIdx < nG; gIdx++ )
    {
      validObs_pg(p,gIdx) = int(0);
      // Stock indices (concentrate obs error var and lnq)
      // Check if this is a survey
      if( calcIndex_g(gIdx) == 1)
      {
        // Loop over time steps
        for( int t = tInitModel_p(p); t < nT; t++)
        {
          Type idxState = 0.;
          if( I_pgt(p,gIdx,t) >= 0.0 )
          {
            // Recover numbers or biomass
            if( survType_g(gIdx) == 1 )
              idxState = vulnN_apgt.col(t).col(gIdx).col(p).sum();
            
            if( survType_g(gIdx) == 0 )
              idxState = vulnB_pgt(p,gIdx,t);
            
            if( survType_g(gIdx) == 2 )
            {
              idxState = (spawnN_apt.col(t).col(p)*mat_a * W_apt.col(t).col(p)).sum();
            }


            // First, calculate probability
            if( deltaIdx_pg(p,gIdx) == 1 )
            {
              Type tmpScalar = 1;
              if( deltaIndVar == 2)
                tmpScalar = B0_p(p);
              Type tmpLogitProb = meanProbPosIdx_pg(p,gIdx) + SDProbPosIdx_pg(p,gIdx) * (idxState/tmpScalar);
              probPosIdx_pgt(p,gIdx,t) = 1 / (1 + exp(-tmpLogitProb ) );
              // Then check if idx is 0 or positive
              if( (I_pgt(p,gIdx,t) == 0) & (probPosIdx_pgt(p,gIdx,t) < 1) )
              {
                // Add bernoulli likelihood here
                obsIdxDeltaNLL_pg(p,gIdx) -= log(1 - probPosIdx_pgt(p,gIdx,t));
              }
            }
            if( deltaIdx_pg(p,gIdx) == 0 )
              probPosIdx_pgt(p,gIdx,t) = 1;
            
            if( I_pgt(p,gIdx,t) > 0)
            {
              // Bernoulli part
              if( deltaIdx_pg(p,gIdx) == 1 )
                obsIdxDeltaNLL_pg(p,gIdx) -= log(probPosIdx_pgt(p,gIdx,t)); 

              // Now apply probability of positive index
              // to expected state variable
              idxState = probPosIdx_pgt(p,gIdx,t) * idxState;
              idxState *= rI_pgt(p,gIdx,t);
              // Calculate residual
              z_pgt(p,gIdx,t) = log(I_pgt(p,gIdx,t)) - log(idxState);
              // Add to sum of residuals
              zSum_pg(p,gIdx) += z_pgt(p,gIdx,t);
              validObs_pg(p,gIdx) += int(1);
            }
          }
        }
        SSR_pg(p,gIdx) = 0.;
        // Calculate conditional MLE of q
        // if a relative index
        if( indexType_g(gIdx) == 0 & validObs_pg(p,gIdx) > 0 )
        {
          // mean residual is lnq (intercept)
          if( qPrior_g(gIdx) == 1 )
            lnqhat_pg(p,gIdx) = (log( mq(gIdx) )/square(sdq(gIdx)) + (zSum_pg(p,gIdx))/tau2Obs_pg(p,gIdx) )/(1/square(sdq(gIdx)) + validObs_pg(p,gIdx)/tau2Obs_pg(p,gIdx));

          if( qPrior_g(gIdx) == 2 )
            lnqhat_pg(p,gIdx) = ( mlnq_g(gIdx)/square(sdlnq_g(gIdx)) + (zSum_pg(p,gIdx))/tau2Obs_pg(p,gIdx) )/(1/square(sdlnq_g(gIdx)) + validObs_pg(p,gIdx)/tau2Obs_pg(p,gIdx));
          
          if( qPrior_g(gIdx) == 0 )
            lnqhat_pg(p,gIdx) = zSum_pg(p,gIdx)/validObs_pg(p,gIdx);
          // Subtract mean from residuals for
          // inclusion in likelihood
          for(int t = 0; t < nT; t++)
            if( I_pgt(p,gIdx,t) > 0.0)
            {
              z_pgt(p,gIdx,t) -= lnqhat_pg(p,gIdx);
            }
            
        }
        // Sum squared resids
        SSR_pg(p,gIdx) += square(z_pgt.transpose().col(p).col(gIdx)).sum();    

        if( validObs_pg(p,gIdx) > 0)
        {
          // Concentrated conditional MLE of observation error
          if( condMLEtauObs == 1 )
          {
            tau2Obs_pg(p,gIdx)    = SSR_pg(p,gIdx) / validObs_pg(p,gIdx);
            tauObs_pg(p,gIdx)     = sqrt(tau2Obs_pg(p,gIdx));
            obsIdxNLL_pg(p,gIdx)  = idxLikeWeight_g(gIdx)*0.5*( validObs_pg(p,gIdx) * log(tau2Obs_pg(p,gIdx)) + validObs_pg(p,gIdx) );
          }

          if( condMLEtauObs == 0 )
          {
            obsIdxNLL_pg(p,gIdx)  = idxLikeWeight_g(gIdx)*0.5*( validObs_pg(p,gIdx) * lntau2Obs_pg(p,gIdx) + SSR_pg(p,gIdx)/tau2Obs_pg(p,gIdx));
          }
        }
      }

      // Catch observations
      if( mortType_g(gIdx) == 1)
      {
        // Need expected catch
        for( int t = 0; t < nT; t++)
          if( C_pgt(p,gIdx,t) > 0 | mC_gt(gIdx,t) > 0 )
          {
            catRes_pgt(p,gIdx,t) = (log(totC_pgt(p,gIdx,t)) - log(expC_pgt(p,gIdx,t)))/resCatCV_g(gIdx);

            // catObsNLL_pg(p,gIdx) -= catLikeWeight_g(gIdx) * dnorm( catRes_pgt(p,gIdx,t),Type(0),Type(1),true);
            catObsNLL_pg(p,gIdx) += pow(catRes_pgt(p,gIdx,t),2);
            if( t == tInitModel_p(p))
              priorqF_pg(p,gIdx) -= dnorm(lnqFinit_pg(p,gIdx),mlnqF_g(gIdx), sdlnqF_g(gIdx),true);
            if( t > tInitModel_p(p))
              priorqF_pg(p,gIdx) -= dnorm(deltaqF_pgt(p,gIdx,t),Type(0), Type(1),true);
          }
      }


      // Age observations //
      // Loop over time steps
      for( int t = tInitModel_p(p); t < nT; t++ )
      { 
        // Check that age observations
        // exist by checking that the plus
        // group has observations
        if( A_apgt(nA-1,p,gIdx,t) >= 0 )
        {
          Type         sumPropAge = 0.;
          // Now estimate predicted catch-at-age
          // This was done already if catch > 0, but 
          // not for all fleets (fishery indep. surveys
          // would be missed)
          // First, calculate prop at age in each fleet
          for(int a = 0; a < nA; a++ )
          {
            predPA_apgt(a,p,gIdx,t)   = uAge_apgt(a,p,gIdx,t);
            // Convert to numbers
            predPA_apgt(a,p,gIdx,t)   /= W_apgt(a,p,gIdx,t);  
            // if( survType_g(gIdx) == 2)
            //   predPA_apgt(a,p,gIdx,t) *= mat(a);

            sumPropAge                += predPA_apgt(a,p,gIdx,t);
          }

          int minAge = minAge_g(gIdx);

          // Save to array and renormalise
          predPA_apgt.col(t).col(gIdx).col(p) /= sumPropAge;

          vector<Type> obsAge = A_apgt.col(t).col(gIdx).col(p).segment(minAge-1,nA-minAge+1);
          vector<Type> predAge = predPA_apgt.col(t).col(gIdx).col(p).segment(minAge-1,nA-minAge+1);

          // Calculate logistic normal likelihood components
          // if age observations exist this year
          if( A_apgt(0,p,gIdx,t) >= 0)
          {
            // ageResids_apgt.col(t).col(gIdx).col(p) = 
            //     calcLogistNormLikelihood(   obsAge, 
            //                                 predAge,
            //                                 minPropAge,
            //                                 etaSumSq_pg(p,gIdx),
            //                                 nResids_pg(p,gIdx) );

            matrix<Type> tmpCorr_aa = Corr_gaa.transpose().col(gIdx).matrix();

            

            vector<Type> tmptcComps(nA - minAge + 1);
            tmptcComps.setZero();
            vector<Type> tmptcPred(nA - minAge + 1);
            tmptcPred.setZero();

            Type tmpgmObs = 0;
            Type tmpgmPred = 0;
            int tmpnBins = 0;
            Type tmplogdetV = 0;
            
            matrix<Type> tmpVchk(nA,nA);
            tmpVchk.setZero();

            matrix<Type> tmpCorr(nA,nA);
            tmpCorr.setZero();

            matrix<Type> tmpKmat(nA,nA);
            tmpKmat.setZero();

            matrix<Type> tmpHmat(nA,nA);
            tmpHmat.setZero();

            matrix<Type> tmpFmat(nA,nA);
            tmpFmat.setZero();

            matrix<Type> tmpGamma(nA,nA);
            tmpGamma.setZero();

            ageResids_apgt.col(t).col(gIdx).col(p).segment(minAge - 1,nA - minAge + 1) =   
                calcCorrLogistNormLikelihood( obsAge, 
                                              predAge,
                                              minPropAge,
                                              etaSumSq_pg(p,gIdx),
                                              nResids_pg(p,gIdx),
                                              meanSampSize_pg(p,gIdx),
                                              compLikeFun,
                                              tmpCorr_aa,
                                              intrnlAgeLikCpt_pg(p,gIdx),
                                              ageWt_pgt(p,gIdx,t),
                                              tmptcComps,
                                              tmptcPred,
                                              tmpgmObs,
                                              tmpgmPred,
                                              tmpnBins,
                                              tmplogdetV,
                                              tmpVchk,
                                              tmpCorr,
                                              tmpKmat,
                                              tmpHmat,
                                              tmpFmat,
                                              tmpGamma  );

            tcComps_apgt.col(t).col(gIdx).col(p).segment(minAge_g(gIdx) - 1,nA - minAge_g(gIdx) + 1) = tmptcComps;
            tcPred_apgt.col(t).col(gIdx).col(p).segment(minAge_g(gIdx) - 1,nA - minAge_g(gIdx) + 1) = tmptcPred;
            if(minAge_g(gIdx) > 1)
            {
              tcPred_apgt.col(t).col(gIdx).col(p).segment(0,minAge_g(gIdx)-1) = 0;
              tcComps_apgt.col(t).col(gIdx).col(p).segment(0,minAge_g(gIdx)-1) = 0;
            }

            nBins_pgt(p,gIdx,t) = tmpnBins;
            logdetV_pgt(p,gIdx,t) = tmplogdetV;

            Vchk_aapgt.col(t).col(gIdx).col(p)  = tmpVchk.array();
            Corr_aapgt.col(t).col(gIdx).col(p)  = tmpCorr.array();
            K_aapgt.col(t).col(gIdx).col(p)     = tmpKmat.array();
            H_aapgt.col(t).col(gIdx).col(p)     = tmpHmat.array();
            F_aapgt.col(t).col(gIdx).col(p)     = tmpFmat.array();
            Gamma_aapgt.col(t).col(gIdx).col(p) = tmpGamma.array();


            // Save geometric means
            gmObs_pgt(p,gIdx,t) = tmpgmObs;
            gmPred_pgt(p,gIdx,t) = tmpgmPred;

            nObsAge_pg(p,gIdx) += 1;
          }
        }
      }
      // Add contribution to age comps likelihood
      if( nResids_pg(p,gIdx) > 0)
      {
        tau2Age_pg(p,gIdx)    = etaSumSq_pg(p,gIdx) / nResids_pg(p,gIdx);
        ageObsNLL_pg(p,gIdx)  += ageCompWeight_g(gIdx) * ( 
                                0.5 * (nResids_pg(p,gIdx)) * log(tau2Age_pg(p,gIdx)) +
                                intrnlAgeLikCpt_pg(p,gIdx) +
                                0.5 * etaSumSq_pg(p,gIdx) / tau2Age_pg(p,gIdx) ) ;
      }
    }

    // Now do combined survey index
    for( int t = tInitModel_p(p); t < nT; t++ )
    {
      Type expIdx = 0.;
      qComb_pt(p,t) = 0.;
      tauComb_pt(p,t) = 0.;

      if( combI_pt(p,t) >= 0.0 )
      {
        // Add up q contributions
        for( int g = 0; g < nG; g++ )
        {
          qComb_pt(p,t) += whichCombIdx_g(g) * qComb_pg(p,g) * rI_pgt(p,g,t);
          tauComb_pt(p,t) += whichCombIdx_g(g) * tauComb_pg(p,g) * tauComb_pg(p,g) * rI_pgt(p,g,t);

        }
        tauComb_pt(p,t) = sqrt(tauComb_pt(p,t));
        // Take sqrt of tauComb
        // tauComb_pt(p,t) = sqrt(tauComb_pt(p,t));

        // First, calculate probability
        // HACKING TOGETHER RIGHT NOW USING DIVE SURVEY PARS
        if( deltaIdx_pg(p,4) == 1 )
        {
          Type tmpLogitProb = meanProbPosIdx_pg(p,4) + SDProbPosIdx_pg(p,4) * SB_pt(p,t);
          probPosCombIdx_pt(p,t) = 1 / (1 + exp(-tmpLogitProb ) );
          // Then check if idx is 0 or positive
          if( (combI_pt(p,t) == 0) & (probPosCombIdx_pt(p,t) < 1) )
          {
            // Add bernoulli likelihood here
            obsCombIdxDeltaNLL_p(p) -= log(1 - probPosCombIdx_pt(p,t));
          }
        }
        if( deltaIdx_pg(p,4) == 0 )
          probPosCombIdx_pt(p,t) = 1;
        
        if( combI_pt(p,t) > 0 )
        {
          // Bernoulli part
          if( deltaIdx_pg(p,4) == 1 )
            obsCombIdxDeltaNLL_p(p) -= log(probPosCombIdx_pt(p,t)); 

          // Now apply probability of positive index
          // to expected state variable
          expIdx = qComb_pt(p,t) * SB_pt(p,t) * probPosCombIdx_pt(p,t);
          // Calculate residual
          zComb_pt(p,t) += log(expIdx) - log(combI_pt(p,t));
          // Add to likelihood function value
          obsCombIdxNLL_p(p) -= dnorm( zComb_pt(p,t), Type(0), tauComb_pt(p,t), true);

          zComb_pt(p,t) /= tauComb_pt(p,t);

        }

      }
    }
  }
  

  // Now add contribution of mixed indices
  // and compositions to the likelihood
  vector<Type>  lnqhat_g(nG);
  vector<Type>  qhat_g(nG);
  vector<Type>  lntauObs_g(nG);
  array<Type>   z_gt(nT,nG);
  vector<Type>  zSum_g(nG);
  vector<Type>  validObs_g(nG);
  lnqhat_g.fill(0.0);
  z_gt.fill(0.0);
  zSum_g.fill(0.0);
  validObs_g.fill(0);
  
  // Age observations
  vector<Type> etaSumSq_g(nG);
  vector<Type> nResids_g(nG);
  vector<Type> nObsAge_g(nG);
  etaSumSq_g.setZero();
  nResids_g.setZero();
  nObsAge_g.setZero();

  array<Type>   predMixedPA_agt(nA,nG,nT);
  predMixedPA_agt.setZero();

  for( int gIdx = 0; gIdx < nG; gIdx++ )
  {
    validObs_g(gIdx) = int(0);
    // Stock indices (concentrate obs error var and lnq)
    // Check if this is a survey
    if( calcIndex_g(gIdx) == 1)
    {
      // Loop over time steps
      for( int t = 0; t < nT; t++)
      {
        Type idxState = 0.;
        if( mI_gt(gIdx,t) > 0.0 )
        {
          for( int p = 0; p < nP; p++ )
          {
            // Recover numbers or biomass
            if( survType_g(gIdx) == 1 )
              idxState += vulnN_apgt.col(t).col(gIdx).col(p).sum();
            if( survType_g(gIdx) == 0 )
              idxState += vulnB_pgt(p,gIdx,t);
            if( survType_g(gIdx) == 2 )
              idxState += SB_pt(p,t);
          }

          // Calculate residual
          z_gt(gIdx,t) = log(mI_gt(gIdx,t)) - log(idxState);
          // Add to sum of residuals
          zSum_g(gIdx) += z_gt(gIdx,t);
          validObs_g(gIdx) += int(1);
        }
      }
      SSR_g(gIdx) = 0.;
      // Calculate conditional MLE of q
      // if a relative index
      if( indexType_g(gIdx) == 0 & validObs_g(gIdx) > 0 )
      {
        // mean residual is lnq (intercept)
        if( qPrior_g(gIdx) == 1)
          lnqhat_g(gIdx)  = (log( mq(gIdx) )/square(sdq(gIdx)) + zSum_g(gIdx)/tau2Obs_g(gIdx) )/(1/square(sdq(gIdx)) + validObs_g(gIdx)/tau2Obs_g(gIdx));

        if( qPrior_g(gIdx) == 0 )
          lnqhat_g(gIdx) = zSum_g(gIdx)/validObs_g(gIdx);

        qhat_g(gIdx)    = exp(lnqhat_g(gIdx));

        // Subtract mean from residuals for
        // inclusion in likelihood
        for(int t = 0; t < nT; t++)
          if( mI_gt(gIdx,t) > 0.0)
          {
            z_gt(gIdx,t) -= lnqhat_pg(gIdx);
          }
          
      }
      // Sum squared resids
      SSR_g(gIdx) += square(z_gt.col(gIdx)).sum();    

      if( validObs_g(gIdx) > 0)
      {
        // Add concentrated nll value using cond MLE of tauObs
        if( idxLikeWeight_g(gIdx) > 0)
          obsMixedIdxNLL_g(gIdx) += 0.5*( lntau2Obs_g(gIdx) + SSR_g(gIdx)/tau2Obs_g(gIdx));
      }
    }

    // Age observations //
    // Loop over time steps
    for( int t = 0; t < nT; t++ )
    { 
      // Check that age observations
      // exist by checking that the plus
      // group has non-negative value (negative => missing)
      if( mA_agt(nA-1,gIdx,t) >= 0 )
      {
        // Now estimate predicted catch-at-age
        // This was done already if catch > 0, but 
        // not for all fleets (fishery indep. surveys
        // would be missed)

        // First, check if catAge 4 is > 0 for any stocks
        // in this fleet at this time step (choosing age
        // 4 as this is a selected age)
        vector<Type> checkCatAge(nP);
        checkCatAge = catAge_apgt.col(t).col(gIdx).transpose().col(3);

        if( checkCatAge.sum() > 0 )
        {
          // make predicted proportions at age
          // from the sum of the catch-at-age
          for( int p = 0; p < nP; p++ )
            predMixedPA_agt.col(t).col(gIdx) += catAge_apgt.col(t).col(gIdx).col(p);
        }

        // Otherwise, we take the average prop, weighted
        // by vulnerable numbers at age
        if( checkCatAge.sum() == 0 )
        {
          for( int a = 0; a < nA; a++ )
            for( int p = 0; p < nP; p++ )
              predMixedPA_agt(a,gIdx,t) += uAge_apgt(a,p,gIdx,t) * vulnN_apgt(a,p,gIdx,t) / mixedVulnN_agt(a,gIdx,t);

        }

        // Save to array and renormalise
        // if( survType_g(gIdx) == 2 ) // spawn index has mature fish only
        //   predMixedPA_agt.col(t).col(gIdx) *= mat;

        predMixedPA_agt.col(t).col(gIdx) /= predMixedPA_agt.col(t).col(gIdx).sum();  

        int minAge  = minAge_g(gIdx);
        vector<Type> obsAge = mA_agt.col(t).col(gIdx).segment(minAge-1,nA-minAge+1);
        vector<Type> predAge = predMixedPA_agt.col(t).col(gIdx).segment(minAge-1,nA-minAge+1);

        // Calculate logistic normal likelihood components
        // if age observations exist this year
        if( mA_agt(0,gIdx,t) >= 0)
        {
          // ageResids_agt.col(t).col(gIdx) = 
          //     calcLogistNormLikelihood(   obsAge, 
          //                                 predAge,
          //                                 minPropAge,
          //                                 etaSumSq_g(gIdx),
          //                                 nResids_g(gIdx) );

          

          matrix<Type> tmpCorr_aa = Corr_gaa.transpose().col(gIdx).matrix();
          vector<Type> tmptcComps(nA - minAge + 1);
          tmptcComps.setZero();
          vector<Type> tmptcPred(nA - minAge + 1);
          tmptcPred.setZero();

          Type tmpgmObs = 0;
          Type tmpgmPred = 0;
          int tmpnBins = 0;
          Type tmplogdetV = 0;

          matrix<Type> tmpVchk(nA,nA);
          tmpVchk.setZero();

          matrix<Type> tmpCorr(nA,nA);
          tmpCorr.setZero();

          matrix<Type> tmpKmat(nA,nA);
          tmpKmat.setZero();

          matrix<Type> tmpHmat(nA,nA);
          tmpHmat.setZero();

          matrix<Type> tmpFmat(nA,nA);
          tmpFmat.setZero();

          matrix<Type> tmpGamma(nA,nA);
          tmpGamma.setZero();

          ageResids_agt.col(t).col(gIdx).segment(minAge - 1,nA - minAge + 1) =
              calcCorrLogistNormLikelihood( obsAge, 
                                            predAge,
                                            minPropAge,
                                            etaSumSq_g(gIdx),
                                            nResids_g(gIdx),
                                            meanSampSize_g(gIdx),
                                            compLikeFun,
                                            tmpCorr_aa,
                                            intrnlAgeLikCpt_g(gIdx),
                                            ageWt_gt(gIdx,t),
                                            tmptcComps,
                                            tmptcPred,
                                            tmpgmObs,
                                            tmpgmPred,
                                            tmpnBins,
                                            tmplogdetV,
                                            tmpVchk,
                                            tmpCorr,
                                            tmpKmat,
                                            tmpHmat,
                                            tmpFmat,
                                            tmpGamma  );

          tcComps_agt.col(t).col(gIdx).segment(minAge - 1,nA - minAge + 1) = tmptcComps;
          
          tcPred_agt.col(t).col(gIdx).segment(minAge - 1,nA - minAge + 1) = tmptcPred;
          if( minAge > 1)
          {
            tcComps_agt.col(t).col(gIdx).segment(0,minAge-1) = 0;
            tcPred_agt.col(t).col(gIdx).segment(0,minAge-1) = 0;
          }
          
          nObsAge_g(gIdx) += 1;
        }
      }
    }
    // Add contribution to age comps likelihood
    if( nResids_g(gIdx) > 0)
    {
      tau2Age_g(gIdx)    += etaSumSq_g(gIdx) / nResids_g(gIdx);
      ageObsNLL_g(gIdx)  += ageCompWeight_g(gIdx) * ( 
                                0.5 * (nResids_g(gIdx)) * log(tau2Age_g(gIdx)) +
                                intrnlAgeLikCpt_g(gIdx) +
                                0.5 * nResids_g(gIdx) ) ;
    }
  }


  // transform stuff
  qhat_g = exp(lnqhat_g);
  qhat_pg = exp(lnqhat_pg);

  // Process error priors
  // Recruitment priors
  // Add recruitment deviations to rec NLL; sigmaR is estimated
  array<Type> SRdevs_pt(nP,nT);     // devs off SR model
  array<Type> Rdevs_pt(nP,nT);      // random effects - may be either Rbar or bhR devs
  vector<Type> meanRdev_p(nP);
  SRdevs_pt.setZero();
  Rdevs_pt.setZero();
  meanRdev_p.setZero();
  

  
  for(int p = 0; p < nP; p++)
  {
    int nDev = 0;
    if( initCode_p(p) == 1)
    {
      vector<Type> x = fDevs_ap.col(p);
      init_nlp -= dnorm( x, Type(0), Type(1), true).sum();
    }
    
    for( int t = tInitModel_p(p)+1; t < nT; t++ )
    {
      if( avgRcode_p(p) == 1 )
      {
        SRdevs_pt(p,t) = (log(R_pt(p,t)) - log(bhR_pt(p,t)))/sigmaR;
        Rdevs_pt(p,t) = omegaR_pt(p,t);

        if( t >= firstRecDev_p(p) & t < lastRecDev_p(p) )
        {
          rec_nlp -= dnorm(Rdevs_pt(p,t), Type(0), Type(1),true);
          rec_nlp -= log(2*omegaRradius) + 2 * log(1 + exp(-recDevs_pt(p,t-1) ));
        }
      }
      if( avgRcode_p(p) == 0 )
      {
        SRdevs_pt(p,t) = omegaR_pt(p,t);
        Rdevs_pt(p,t) = omegaR_pt(p,t);
      }

      if( t >= firstRecDev_p(p) & t < lastRecDev_p(p) )
      {
        nDev++;
        meanRdev_p(p) += SRdevs_pt(p,t);
      }
    }
    meanRdev_p(p) /= nDev;
  }

  // Now apply correlated mvn likelihoods for process error devs
  MVNORM_t<Type> SR_mvn(sigmaR_pp);
  MVNORM_t<Type> tvM_mvn(corrM_pp);
  vector<Type> indInitPop_p(nP);
  vector<Type> sigmaM_p(nP);
  vector<Type> sigmaR_p(nP);
  sigmaM_p.fill(sigmaM);
  sigmaR_p.fill(sigmaR);
  for( int t=0; t < nT; t++ )
  {
    indInitPop_p.fill(0);
    for( int p = 0; p < nP; p++)
      if( t >= firstRecDev_p(p) & t < lastRecDev_p(p) )
        indInitPop_p(p) = 1;

    if( indInitPop_p.prod() == 1 )
    {
      vector<Type> srParVec = SRdevs_pt.col(t).vec() * sigmaR_p;
      if( corrRdevs == 1 )
        SRnlp += SR_mvn(srParVec);
      
      if( corrRdevs == 0 )
      {
        Type tmpsigR = sigmaR_p(0);
        SRnlp -= dnorm( srParVec, Type(0), tmpsigR, true).sum();
      }
      
      // add sigma - avoiding
      // SRnlp += nP * lnsigmaR;

      if(avgRcode_p.sum() != 0)
        SRnlp -= nP * log(2*omegaRradius) + 2 * log(1 + exp(-recDevs_pt.col(t-1).vec()  )).sum();
      
      if(densityDepM == 0 & t < nT )
      {
        vector<Type> mortParVec = omegaM_pt.col(t).vec();
        if( corrMdevs == 1 )
          tvMnlp += tvM_mvn(mortParVec);

        if( corrMdevs == 0 )
          tvMnlp -= dnorm( mortParVec, Type(0), Type(1), true).sum();

      }
      
    }

    if( indInitPop_p.prod() == 0 )
      for( int p = 0; p < nP; p ++ )
      {
        if( t >= firstRecDev_p(p) )
        {
          SRnlp   -= dnorm( SRdevs_pt(p,t) * sigmaR_p(p), Type(0), sigmaR_p(p), true);
          
          if(avgRcode_p(p) != 0)
            SRnlp   -= log(2*omegaRradius) + 2 * log(1 + exp(-recDevs_pt(p,t-1) ));
        }
          
        if(densityDepM == 0 & t < nT & t > tInitModel_p(p) )
          tvMnlp  -= dnorm( omegaM_pt(p,t), Type(0), Type(1), true);        
        
        
        
      }
  }
  
  // Natural mortality hyper-priors
  if( densityDepM == 0 )
  {
    mort_nlp -= dnorm( lnM, log(Mprior(0)),Mprior(1), true);
  }
  
  Type sigRnlp = 0;
  sigRnlp += log(priorSigR) + 0.5*square(sigmaR)/square(priorSigR);
  
  if( densityDepM == 1 )
  {
    mort_nlp -= dnorm( log( M0_p ), log(Mprior(0)), Mprior(1),true).sum();
    mort_nlp += m1_p.sum() + log(M0_p).sum() - lnm1; 
    
    lnm1_nlp -= dnorm( lnm1, log(m1Prior(0)), m1Prior(1), true);    

    for( int t = 0; t < nT; t ++)
      tvMnlp -= dnorm( omegaM_pt.col(t).vec(), Type(0), Type(1), true).sum();
  }

  if( nP > 1)
  {
    Mdev_nlp -= dnorm( epsM_p, Type(0), Type(1), true).sum();
    
    if( densityDepM == 1)
      lnm1_nlp -= dnorm( epslnm1_p, Type(0), Type(1), true).sum();
  }
  




  // Catchability hyperprior
  Type qnlp = 0;
  for(int g = 0; g < nG; g++)
  {
    if(qPrior_g(g) == 2)
    {
      Type resid = (mlnq_g(g) - log(mq(g)))/sdq(g);
      qnlp -= dnorm(resid, Type(0), Type(1), true);
    }

    if( whichCombIdx_g(g) == 1 )
    {
      // Overwrite the cond MLE of q with the
      // combo q
      qhat_pg.col(g) = qComb_pg.col(g);

      if( qPrior_g(g) == 1 )
      {
        vector<Type> resid_p = (lnqComb_pg.col(g) - log(mq(g)))/sdq(g);
        qnlp -= dnorm( resid_p, Type(0), Type(1), true).sum();
      }

      if( qPrior_g(g) == 2 )
      {
        vector<Type> resid_p = (lnqComb_pg.col(g) - mlnq_g(g) )/sdlnq_g(g);
        qnlp -= dnorm( resid_p, Type(0), Type(1), true).sum();

        Type resid = (mlnq_g(g) - log(mq(g)))/sdq(g);
        qnlp -= dnorm(resid, Type(0), Type(1), true);
      }
    }
  }
      

  // Beta prior on steepness
  h_nlp  = (1 - rSteepBetaPrior(0)) * log(rSteepness) + (1 - rSteepBetaPrior(1)) * log(1-rSteepness) - log(0.78) - 2*log(1 + exp(-logit_ySteepness));
  // Add a prior on the steepness stock effects
  hDev_nlp -= dnorm( epsSteep_p, Type(0), Type(1), true).sum();
  // Add a recruitment var IG prior
  // Type sig2R_nlp = 0*(sig2RPrior(0)+Type(1))*2*lnsigmaR + sig2RPrior(1)/square(sigmaR);

  // SOK-Ponded biomass conversion factor prior
  Type psi_nlp = 0;
  // for( int g = 0; g < nG; g++ )
  //   for( int t = 0; t < nT; t++ )
  //   {
  //     if( psi_gt(g,t) > 0 )
  //       psi_nlp -= dnorm( log(psi_gt(g,t)), log(mPsi), sdPsi, true);

  //     for( int p = 0; p < nP; p++ )
  //       if( psi_pgt(p,g,t) > 0 )
  //         psi_nlp -= dnorm( log(psi_pgt(p,g,t)), log(mPsi), sdPsi, true);
  //   }

  // Add gear-wise priors:
  // Catchability, index obs error variance,
  // and a prior on selectivity pars
  vector<Type> selAlphaNLP_g(nG);
  vector<Type> selBetaNLP_g(nG);
  array<Type>  selAlphaDevNLP_pg(nP,nG);
  array<Type>  selBetaDevNLP_pg(nP,nG);
  Type tvselAlphaDevNLP = 0;
  Type tvselBetaDevNLP = 0;
  selAlphaDevNLP_pg.setZero();
  selBetaDevNLP_pg.setZero();
  selAlphaNLP_g.setZero();
  selBetaNLP_g.setZero();

  Type lnqF_nlp = 0;

  // Time varying selectivity devs
  tvselAlphaDevNLP -= dnorm(epsSelAlpha_vec, Type(0), Type(1),true ).sum();
  tvselBetaDevNLP -= dnorm(epsSelBeta_vec, Type(0), Type(1),true ).sum();
  for( int gIdx = 0; gIdx < nG; gIdx++ )
  {
    
    if( hierSel == 1 )
    {
      selAlphaNLP_g(gIdx) -= dnorm( lnSelAlpha_g(gIdx), mlnSelAlpha_g(gIdx), sdSel_g(gIdx), true);
      selBetaNLP_g(gIdx) -= dnorm( lnSelBeta_g(gIdx), mlnSelBeta_g(gIdx), sdSel_g(gIdx), true);
    }

    for( int p=0; p < nP; p++)
    {
      // IG prior on survey obs err variance
      if( calcIndex_g(gIdx) == 1 & validObs_pg(p,gIdx) > 0 & condMLEtauObs == 0  )
      {
        nlptau2idx_pg(p,gIdx) += (obstau2IGa(gIdx)+Type(1))*lntau2Obs_pg(p,gIdx) + obstau2IGb(gIdx)/square(tauObs_pg(p,gIdx));
        nlptau2idx_pg(p,gIdx) -= log(Type(2)) + lntau2Obs_pg(p,gIdx); 
      }

      // Add stock effects for selectivity
      if( hierSel == 1 )
      {
        selAlphaDevNLP_pg(p,gIdx) -= dnorm( epsSelAlpha_pg(p,gIdx), Type(0), Type(1), true);
        selBetaDevNLP_pg(p,gIdx) -= dnorm( epsSelBeta_pg(p,gIdx), Type(0), Type(1), true);
      }

      if( hierSel == 0 )
      {
        selAlphaNLP_g(gIdx) -= dnorm( epsSelAlpha_pg(p,gIdx), mlnSelAlpha_g(gIdx), sdSel_g(gIdx), true);
        selBetaNLP_g(gIdx)  -= dnorm( epsSelBeta_pg(p,gIdx), mlnSelBeta_g(gIdx), sdSel_g(gIdx), true);
      }

    }
    
    // Now for the obs index variance for mixed indices
    if( calcIndex_g(gIdx) == 1 & validObs_g(gIdx) > 0 & condMLEtauObs == 0  )
    {
      nlptau2idx_g(gIdx) += (obstau2IGa(gIdx)+Type(1))*lntau2Obs_g(gIdx) + obstau2IGb(gIdx)/tau2Obs_g(gIdx);
      nlptau2idx_g(gIdx) -= log(Type(2)) + lntau2Obs_g(gIdx);
    }

  }

  // IG priors for obs error variance on the convex combination of indices
  if( whichCombIdx_g.sum() > 0 )
    for( int p = 0; p < nP; p++ )
      for( int g = 0; g < nG; g++)
        if( whichCombIdx_g(g) > 0)
        {
          nlptau2idx_pg(p,g)  += (obstau2IGa(g)+Type(1))*lntauObsComb_pg(p,g) + obstau2IGb(g)/square(tauComb_pg(p,g));
          nlptau2idx_pg(p,g)  -= log(Type(2)) + 2*lntauObsComb_pg(p,g);
        }

  Type logitProbPosIdx_nlp = 0;
  for( int p = 0; p < nP; p++ )
    for( int g = 0; g < nG; g++ )
      if( deltaIdx_pg(p,g) == 1)
      {
        logitProbPosIdx_nlp -= dnorm(lnSDProbPosIdx_pg(p,g),muSDProbPosIdx_g(g),sigSDProbPosIdx_g(g),true);
        logitProbPosIdx_nlp -= dnorm(meanProbPosIdx_pg(p,g),muMeanProbPosIdx_g(g),sigMeanProbPosIdx_g(g),true);
       }
 
  
  // Add positive function penalty to objFun
  objFun += posPenFactor * posPen;

  Type totLike =  obsIdxNLL_pg.sum() +
                  obsIdxDeltaNLL_pg.sum() +
                  obsCombIdxNLL_p.sum() +
                  obsCombIdxDeltaNLL_p.sum() +
                  ageObsNLL_pg.sum() +
                  obsMixedIdxNLL_g.sum() +
                  ageObsNLL_g.sum()+ 
                  catObsNLL_pg.sum();

  Type totPriorDens = 0;

  totPriorDens =  mort_nlp + 
                  tvMnlp +
                  Mdev_nlp +
                  rec_nlp + 
                  init_nlp +
                  h_nlp +
                  hDev_nlp +
                  qnlp +
                  priorqF_pg.sum() +
                  nlptau2idx_pg.sum() +
                  nlptau2idx_g.sum() +
                  tvselAlphaDevNLP +
                  tvselBetaDevNLP +
                  selAlphaDevNLP_pg.sum() +
                  selBetaDevNLP_pg.sum()+
                  (selAlphaNLP_g * selPriorWt_g).sum() +
                  (selBetaNLP_g * selPriorWt_g).sum() +
                  jeffWtB0 * lnB0_p.sum() +
                  jeffWtB0 * lnRinit_p.sum() +
                  jeffWtB0 * lnRbar_p.sum() +
                  lnm1PriorWt * lnm1_nlp +
                  meanRdevWt * square(meanRdev_p).sum() +
                  psi_nlp +
                  SRnlp +
                  sigRnlp +
                  logitProbPosIdx_nlp;




  if( compLikeFun >= 1 )
    totPriorDens += (logitphi1_g*logitphi1_g).sum();

  if( compLikeFun >= 2 )
    totPriorDens += (logitpsi_g*logitpsi_g).sum();

  totPriorDens += corrParWeight * ( (off_diag_M*off_diag_M).sum() +
                              (off_diag_R*off_diag_R).sum() );


  // Add NLL and NLP contributions to objFun
  objFun += dataLikeWt * totLike +
            priorDensWt * totPriorDens;
  

  // Convert some values to log scale
  // for sd reporting
  array<Type> lnprojSB_pv(nP,nV);
  array<Type> lnprojDep_pv(nP,nV);
  
  // Fill arrays
  for( int v = 0; v < nV; v++)
  {
    lnprojSB_pv.col(v)  = log(catTable_pvk.col(1).col(v));
    lnprojDep_pv.col(v) = log(catTable_pvk.col(2).col(v));
  }



  /*\/\/\/\/\ REPORTING SECTION /\/\/\/\*/
  // Variables we want SEs for
  ADREPORT(lnprojSB_pv);
  ADREPORT(lnprojDep_pv);
  

  // Everything else //

  // Leading pars
  REPORT(B0_p);                   // Unfished biomass
  REPORT(totB0_p);                // Unfished total biomass
  REPORT(Rinit_p);                // Initial rec
  REPORT(Rbar_p);                 // Avg R
  REPORT(M_p);                    // Initial/constant natural mortality
  REPORT(M);                      // Species M
  REPORT(m1_p);                   // DDM scaling factor
  REPORT(rSteepness);             // Species steepness
  REPORT(initN_mult_ap);          // Initial numbers multiplier (fished)
  REPORT(fDevs_ap);
  REPORT(rSteepness_p);           // Steepness
  REPORT(sigmaM);                 // Mt devations sd
  REPORT(sigmaR);                 // Recruitment deviations SD
  REPORT(Finit_p);                // Initialisation F
  REPORT(Mjuve_p);                // Juvenile M (1-juveMage)
  REPORT(juveMage);               // Juvenile M (1-3)


  // Growth and maturity schedules
  REPORT(mat_a);
  REPORT(wt_a);
  REPORT(meanWt_ap);
  REPORT(projWt_ap);
  REPORT(meanWt_agp);
  REPORT(projWt_agp);
  REPORT(meanM_p);
  REPORT(M0_p);
  REPORT(len_a);
  REPORT(age_a);
  REPORT(surv_ap);
  REPORT(initSurv_ap);

  // Model dimensions
  REPORT(nT);
  REPORT(nG);
  REPORT(nA);
  REPORT(nP);
  REPORT(nV);


  // Selectivity
  REPORT(sel_ag);
  REPORT(sel_apgt);
  // First branch selectivity parameters (normal dome, increasing, or decreasing)
  REPORT(SelAlpha_g);
  REPORT(SelBeta_g);
  REPORT(SelAlpha_pg);
  REPORT(SelBeta_pg);
  REPORT(SelAlpha_pgt);
  REPORT(SelBeta_pgt);
  REPORT(epsSelAlpha_pgt);
  REPORT(epsSelBeta_pgt);
  REPORT(epsSelAlpha_pg);
  REPORT(epsSelBeta_pg);

  // Down selection (only used for asymptotic dome)
  REPORT(dSelAlpha_g);
  REPORT(dSelBeta_g);
  REPORT(dSelAlpha_pg);
  REPORT(dSelBeta_pg);
  REPORT(dSelAlpha_pgt);
  REPORT(dSelBeta_pgt);
  REPORT(epsdSelAlpha_pgt);
  REPORT(epsdSelBeta_pgt);
  REPORT(epsdSelAlpha_pg);
  REPORT(epsdSelBeta_pg);

  REPORT(scaleSel_gt);
  REPORT(scaleSel_g);

  REPORT(sigmaSelAlpha_g);
  REPORT(sigmaSelBeta_g);
  REPORT(selType_g);
  REPORT(selX_g);

  // SR Variables
  REPORT(R0_p);
  REPORT(reca_p);
  REPORT(recb_p);
  REPORT(phi_p);
  REPORT(totbpr_p);

  // State variables
  REPORT(B_apt);
  REPORT(B_pt);
  REPORT(N_apt);
  REPORT(endN_apt);
  REPORT(initZ_ap);
  REPORT(Eggs_pt);
  REPORT(movN_apt);
  REPORT(tmpN_apt);
  REPORT(initN_apt);
  REPORT(termN_apt);
  REPORT(spawnN_apt);
  REPORT(vulnB_pgt);
  REPORT(vulnB_apgt);
  REPORT(vulnN_pgt);
  REPORT(vulnN_apgt);
  REPORT(mixedVulnN_agt);
  REPORT(mixedVulnB_agt);
  REPORT(mixedVulnN_gt);
  REPORT(mixedVulnB_gt);
  REPORT(SB_pt);
  REPORT(R_pt);
  REPORT(M_apt);
  REPORT(U_pgt);
  REPORT(U_apgt);
  REPORT(totSurv_apt);
  REPORT(Z_apt);
  REPORT(appZ_apt);
  REPORT(bhR_pt);
  REPORT(SRdevs_pt);

  // Fleet reordering and 
  // Pope's approx
  REPORT(chronIdx_g);
  REPORT(nDiscFleets);
  REPORT(nCtsFleets);
  REPORT(ctsFleetIdx);
  REPORT(discFleetIdx);
  // REPORT(usedFleet);
  REPORT(fleetTiming_g);
  REPORT(fleetType_g);
  REPORT(uAge_apgt);
  REPORT(catAge_apgt);
  REPORT(predPA_apgt);
  // REPORT(fracM);
  // REPORT(lastFrac);
  // REPORT(prevFleetTime);

  // Data switches
  REPORT(survType_g);     // Type of index (0 = vuln bio, 1 = vuln numbers)
  REPORT(indexType_g);    // Type of survey (0 = relative, 1 = absolute)
  REPORT(calcIndex_g);    // Calculate fleet index (0 = no, yes = 1)
  REPORT(selType_g);      // Type of selectivity (0 = asymptotic, 1 = domed (normal), 2 = decreasing asymp)
  REPORT(minPropAge);     // Min prop'n at age to avoid accumulation in L-N age comp likelihood
  REPORT(densityDepM);    // Calculate DDM
  
  // Data
  REPORT(C_pgt);
  REPORT(I_pgt);
  REPORT(A_apgt);
  REPORT(W_apt);
  REPORT(mC_gt);
  REPORT(mI_gt);
  REPORT(mA_agt);

  // Debugging movement model
  REPORT(movFracM_t);

  // total removals (distributing mixed catch among stocks)
  REPORT( totC_pgt );
  REPORT( splitC_pgt );

  // Catch model
  REPORT( qF_pgt);
  REPORT( F_pgt );
  REPORT( lnqF_pgt);
  REPORT( expC_pgt );
  REPORT( deltaqF_pgt);
  REPORT( lnqFinit_pg );
  REPORT( catRes_pgt );

  // SOK
  REPORT( pondC_apgt );
  REPORT( pondC_pgt );
  REPORT( pondC_gt );
  REPORT( psi_gt );
  REPORT( psi_pgt );
  REPORT( propMat_gt );
  REPORT( propMat_pgt );
  REPORT( gamma_g );
  REPORT( fec );
  REPORT( mPsi );
  REPORT( sdPsi );
  REPORT( pEff_t );
  REPORT( postPondM_g );
  REPORT( pFem );

  // Movement matrix
  REPORT( mov_ppa );

  // Random effects
  REPORT(omegaM_pt);
  REPORT(omegaR_pt);
  REPORT(recDevs_pt);
  REPORT(Rdevs_pt);
  REPORT(meanRdev_p);

  
  // Observation model quantities
  REPORT(qhat_pg);
  REPORT(lnqhat_pg);
  REPORT(tauObs_pg);
  REPORT(tau2Obshat_pg);
  REPORT(z_pgt);
  REPORT(zSum_pg);
  REPORT(zComb_pt);
  REPORT(validObs_pg);
  REPORT(SSR_pg);
  REPORT(etaSumSq_pg);
  REPORT(tau2Age_pg);
  REPORT(nResids_pg);
  REPORT(nObsAge_pg);
  REPORT(ageResids_apgt);
  REPORT(probPosIdx_pgt);
  REPORT(probPosCombIdx_pt);
  REPORT(SDProbPosIdx_pg);
  REPORT(meanProbPosIdx_pg);
  REPORT( catObsNLL_pg );
  REPORT(priorqF_pg);

  // Convex combination of surveys pars
  REPORT(qComb_pt);
  REPORT(qComb_pg);
  REPORT(tauComb_pg);
  REPORT(tauComb_pt);


  // Age comp debugging arrays
  REPORT(ageWt_pgt);
  REPORT(ageWt_gt);
  REPORT(tcComps_apgt);
  REPORT(tcComps_agt);
  REPORT(tcPred_apgt);
  REPORT(tcPred_agt);
  REPORT(gmObs_pgt);
  REPORT(gmPred_pgt);
  REPORT(nBins_pgt);
  REPORT(logdetV_pgt);
  REPORT(Vchk_aapgt);
  REPORT(Corr_aapgt);
  REPORT(H_aapgt);
  REPORT(F_aapgt);
  REPORT(K_aapgt);
  REPORT(Gamma_aapgt);

  // Mixed observations
  REPORT(qhat_g);
  REPORT(lnqhat_g);
  REPORT(tauObs_g);
  REPORT(z_gt);
  REPORT(zSum_g);
  REPORT(validObs_g);
  REPORT(SSR_g);
  REPORT(etaSumSq_g);
  REPORT(tau2Age_g);
  REPORT(nResids_g);
  REPORT(nObsAge_g);
  REPORT(ageResids_agt);
  REPORT(predMixedPA_agt);

  // Switches and weights
  REPORT(ageCompWeight_g);
  REPORT(idxLikeWeight_g);
  REPORT(tInitModel_p);


  // Likelihood/prior values
  REPORT(objFun);
  REPORT(totLike);
  REPORT(totPriorDens);
  REPORT(ageObsNLL_pg);
  REPORT(ageObsNLL_g);
  REPORT(intrnlAgeLikCpt_pg);
  REPORT(intrnlAgeLikCpt_g);
  REPORT(obsIdxNLL_pg);
  REPORT(obsIdxDeltaNLL_pg);
  REPORT(obsMixedIdxNLL_g);
  REPORT(obsCombIdxDeltaNLL_p);
  REPORT(obsCombIdxNLL_p);
  REPORT(rec_nlp);
  REPORT(SRnlp);
  REPORT(mort_nlp);
  REPORT(init_nlp);
  REPORT(tvMnlp);
  REPORT(Mdev_nlp);
  REPORT(lnm1_nlp);
  REPORT(qnlp);
  REPORT(psi_nlp);
  REPORT(init_nlp);
  REPORT(sigRnlp);
  REPORT(nlptau2idx_pg);
  REPORT(nlptau2idx_g);
  REPORT(h_nlp);
  REPORT(hDev_nlp);
  REPORT(posPen);
  REPORT(selAlphaDevNLP_pg);
  REPORT(selBetaDevNLP_pg);
  REPORT(selAlphaNLP_g);
  REPORT(selBetaNLP_g);
  REPORT(tvselAlphaDevNLP);
  REPORT(tvselBetaDevNLP);
  REPORT(obstau2IGa);
  REPORT(obstau2IGb);
  REPORT(corrM_pp);
  REPORT(corrR_pp);
  REPORT(cholM_pp);
  REPORT(cholR_pp);
  REPORT(sigmaM_pp);
  REPORT(sigmaR_pp);
  REPORT(SigmaM_pp);
  REPORT(SigmaR_pp);
  REPORT(logitProbPosIdx_nlp);

  // Dilute catch from SOG
  REPORT(propCatAgeDiluted_apgt);
  REPORT(totCatAge_apgt);
  REPORT(diluteN_apt);
  REPORT(totCdiluted_pgt);
  REPORT(totCundiluted_pgt);
  REPORT(uAgeDiluted_apgt);
  REPORT(catAgeDiluted_apgt);
  REPORT(catAgeUndiluted_apgt);
  REPORT(totCpredilute_pgt);

  REPORT(vulnNdiluted_apgt);
  REPORT(vulnNdiluted_pgt);
  REPORT(vulnBdiluted_apgt);
  REPORT(vulnBdiluted_pgt);

  // Catch table
  REPORT( catTable_pvk );
  REPORT( alloc_pg );

  // age obs correlation mtx parameters
  REPORT(Corr_gaa);
  REPORT(corr_ga);
  REPORT(phi1_g);
  REPORT(phi2_g);
  REPORT(psi_g);
  
  return objFun;
} // END objective function
