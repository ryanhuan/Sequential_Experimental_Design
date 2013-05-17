#include "SARM.h"

/*********************************************************************
 * SARM
 *********************************************************************/

SARM::SARM(Controls const &refControls, int* const nIterRef, 
	   int* const storageCtrRef, int* const iterHistoryRef, 
	   double** const XHistoryRef, double* const valHistoryRef, 
	   double* const gradNormHistoryRef)
{
  
  /* Load controls, and pointers for variables from StochasticSearch
   * class. */
  primary = refControls;
  nIter = nIterRef;
  storageCtr = storageCtrRef;
  iterHistory = iterHistoryRef;
  XHistory = XHistoryRef;
  valHistory = valHistoryRef;
  gradNormHistory = gradNormHistoryRef;
  DPCosts = new DPCostsInsideExpectation(refControls);
  
  flagCheckGradientFD = 0;
  if (flagCheckGradientFD && primary.rank == 0)
    cout << "% Warning: SARM checking gradient using FD at every computation. "
	 << endl;

  /* Allocate memory. */
  XCur = new double[primary.nControlsDim];
  gradTemp = new double[primary.nControlsDim];
  gradTemp2 = new double[primary.nControlsDim];
  gradCur = new double[primary.nControlsDim];

  /* Random number generator initializations. */
  rngType = gsl_rng_ranlxs0;
  generatorNoise1 = gsl_rng_alloc(rngType);
  generatorNoise2 = gsl_rng_alloc(rngType);
  generatorNoiseFDCheck = gsl_rng_alloc(rngType);
  gsl_rng_env_setup();
  gsl_rng_set(generatorNoise1, rand() + algParams.rank);
  gsl_rng_set(generatorNoise2, rand() + algParams.rank);
  gsl_rng_set(generatorNoiseFDCheck, rand() + algParams.rank);

}

SARM::~SARM()
{
  
  /* Free memory. */
  delete DPCosts;
  delete [] XCur;
  delete [] gradTemp;
  delete [] gradTemp2;
  delete [] gradCur;
  gsl_rng_free(generatorNoise1);
  gsl_rng_free(generatorNoise2);
  gsl_rng_free(generatorNoiseFDCheck);
  
}

void SARM::initialPosition(double (*innerFcn)(Controls const &, 
					      double const * const), 
			   double const * const state)
{

  /* Set random seed for a new optimization run. Note technically we
   * do not need CRN for this algorithm, but it may be used for
   * debugging purposes. */
  CRNSeed = rand() + primary.rank;
  gsl_rng_set(generatorNoise1, CRNSeed);
  
  /* Load initial position from controls, save initial position. The
   * initial point is assumed to be valid. */
  for (int i = 0; i < primary.nControlsDim; i++)
  {
    XCur[i] = primary.XInitial[i];
    XHistory[0][i] = primary.XInitial[i];
  }

  /* For FD gradient estimate. */
  gsl_rng_memcpy(generatorNoise2, generatorNoise1);

  /* For validating initial gradient against FD. */
  if (primary.checkInitialGradient || flagCheckGradientFD)
    /* The current generator seed stored. */
    gsl_rng_memcpy(generatorNoiseFDCheck, generatorNoise1);
  
  /* The objective evaluation is ONLY required for objectiveFcn = 0
   * (i.e, expected utility), because its gradient is a byproduct of
   * the objective evaluation. */
  /* Then evaluate the gradient. */ 
  
  valCur = 0.0;
  for (int i = 0; i < primary.nControlsDim; i++)
    gradCur[i] = 0.0;
  for (int i = 0; i < primary.nObjMC; i++)
  {
    valTemp = -DPCosts -> sampleOnce(innerFcn, XCur, generatorNoise1);
    
    DPCosts -> sampleOnceGradientFD(innerFcn,
				    stateKRef, XCur, generatorNoise2, 
				    gradTemp);
    // DPCosts -> sampleOnceGradient(coefsJkp1Ref, refTableJkp1Ref, 
    // 				  stateKRef, XCur, gradTemp);
    
    /* Since the sampleOnce computes the positive J_{k+1} and the
     * gradient, we need to negate the gradient as well. */
    for (int i = 0; i < primary.nControlsDim; i++)
      gradTemp[i] = -gradTemp[i];
    
    /* The value is not used, but might as well record it for
     * debugging purposes. */
    valCur += valTemp;
    
    for (int i = 0; i < primary.nControlsDim; i++)
      gradCur[i] += gradTemp[i];

    /* We are not counting the number of gradient evaluations here,
     * although it can be incorporated. Should 1 gradient evaluation
     * (1 realization) be viewed as equivalent to 1 objective
     * evaluation in terms of computational effort quantification?
     * Perhaps the most straightforward testing method is just to look
     * at wall clock time. */
  }
  /* The value is not computed except for if expected utility is
   * used. */
  valCur /= double(primary.nObjMC);
  valHistory[0] = valCur;
  for (int i = 0; i < primary.nControlsDim; i++)
    gradCur[i] /= double(primary.nObjMC);
  gradNormHistory[0] = normOf1DArray(gradCur, primary.nControlsDim, 
				     primary.gradNormChoice);
  
  /* Validate initial gradient against FD. */
  if (primary.checkInitialGradient || flagCheckGradientFD)
  {
    for (int i = 0; i < primary.nControlsDim; i++)
      gradTemp2[i] = 0.0;

    for (int i = 0; i < primary.nObjMC; i++)
    {
      /* The stored generator seed is used. */
      DPCosts -> sampleOnceGradientFD(coefsJkp1Ref, refTableJkp1Ref, 
				      stateKRef, XCur, generatorNoiseFDCheck, 
				      gradTemp);
      
      /* Since the sampleOnce computes the positive J_{k+1} and the
       * gradient, we need to negate the gradient. */
      for (int i = 0; i < primary.nControlsDim; i++)
	gradTemp[i] = -gradTemp[i];

      for (int i = 0; i < primary.nControlsDim; i++)
	gradTemp2[i] += gradTemp[i];
    }

    for (int i = 0; i < primary.nControlsDim; i++)
      gradTemp2[i] /= double(primary.nObjMC);
    
    cout << "% Initial gradient validation:" << endl;
    cout << "%    Analytic Gradient           FD Gradient" << endl;
    for (int i = 0; i < primary.nControlsDim; i++)
      cout << "% " << setw(23) << gradCur[i] << "  " 
	   << setw(23) << gradTemp2[i] << endl;
    cout << endl;
  }

  /* If the SPSA gain sequences are used. Note that all relevant
   * inputs are under the SPSA section. */
  if (primary.SARMGainSeq == 2)
  {
    /* Initialize SPSA controls parameters from input file. */
    SPSAa = primary.SPSAInputa;
    SPSAA = primary.SPSAInputA;
    SPSAalpha = primary.SPSAInputalpha;
    
    /* If parameter detection requested. */
    if (primary.detectSPSAParams)
    {
      /* Set alpha. */
      SPSAalpha = 0.602;
      
      /* Set A. Recommended to be proportional (e.g., 10% of) to the
       * max number of iterations allowed or expected. */
      SPSAA = 0.1 * double(primary.maxOptIters);
  
      /* Set a. Recommended such that a/(A+1)^alpha |\hat{g}_0(x_0)|
       * (does not have to be norm, can be anything that reflects the
       * order-of-magnitude of the gradient entries) is approx equal
       * to the smallest of the desired change magnitudes in the
       * search space. Let the smallest of the desired change
       * magnitudes be 10% of the box diagonal length. */
      /* If the gradient norm is too small, then SPSAa could become
       * too large, causing slow or lack of convergence. Thus, a lower
       * bound value of 0.1 is set in the denominator. */
      double size(0.0);
      for (int i = 0; i < primary.nControlsDim; i++)
	size += (primary.controlsRightBound[i] - primary.controlsLeftBound[i]) 
	  * (primary.controlsRightBound[i] - primary.controlsLeftBound[i]);
      size = sqrt(size);
      SPSAa = 0.1 * size * pow(SPSAA + 1.0, SPSAalpha) / max(gradNormHistory[0], 0.1);
    }
    
  }
    
}

void SARM::iterateOnce(double const * const coefsJkp1Ref, 
		       int const * const * const refTableJkp1Ref, 
		       double const * const stateKRef)
{
  
  /* Update counter. */
  nIter[0]++;
  storageCtr[0]++;
  iterHistory[storageCtr[0]] = nIter[0];
  
  /* Update gain sequence. */
  update_a_k();

  for (int i = 0; i < primary.nControlsDim; i++)
  {
    /* Get to the next point. */
    XCur[i] -= a_k * gradCur[i];
    
    /* Euclidean projection to box constraints if violating. */
    if (XCur[i] < primary.controlsLeftBound[i])
      XCur[i] = primary.controlsLeftBound[i];
    else if (XCur[i] > primary.controlsRightBound[i])
      XCur[i] = primary.controlsRightBound[i];
    
    /* Store the new position. */
    XHistory[storageCtr[0]][i] = XCur[i];
    
  }
  
  /* For FD gradient estimate. */
  gsl_rng_memcpy(generatorNoise2, generatorNoise1);
  
  /* For validating gradient against FD. */
  if (flagCheckGradientFD)
    /* The current generator seed stored. */
    gsl_rng_memcpy(generatorNoiseFDCheck, generatorNoise1);

  /* Compute new gradient estimate. */
  valCur = 0.0;
  valTemp = 0.0;
  for (int i = 0; i < primary.nControlsDim; i++)
    gradCur[i] = 0.0;
  for (int i = 0; i < primary.nObjMC; i++)
  {
    valTemp = -DPCosts -> sampleOnce(coefsJkp1Ref, refTableJkp1Ref, 
				     stateKRef, XCur, generatorNoise1);

    DPCosts -> sampleOnceGradientFD(coefsJkp1Ref, refTableJkp1Ref, 
				    stateKRef, XCur, generatorNoise2, 
				    gradTemp);

    // DPCosts -> sampleOnceGradient(coefsJkp1Ref, refTableJkp1Ref, 
    // 				  stateKRef, XCur, gradTemp);

    /* Since the sampleOnce computes the positive J_{k+1} and the
     * gradient, we need to negate the gradient. */
    for (int i = 0; i < primary.nControlsDim; i++)
      gradTemp[i] = -gradTemp[i];
    
    /* The value is not used, but might as well record it for
     * debugging purposes. */
    valCur += valTemp;
    
    for (int i = 0; i < primary.nControlsDim; i++)
      gradCur[i] += gradTemp[i];
    
    /* We are not counting the number of gradient evaluations here,
     * although it can be incorporated. Should 1 gradient evaluation
     * (1 realization) be viewed as equivalent to 1 objective
     * evaluation in terms of computational effort quantification?
     * Perhaps the most straightforward testing method is just to look
     * at wall clock time. */
  }
  /* The value is not computed except for if expected utility is
   * used. */
  valCur /= double(primary.nObjMC);
  valHistory[storageCtr[0]] = valCur;
  for (int i = 0; i < primary.nControlsDim; i++)
    gradCur[i] /= double(primary.nObjMC);
  gradNormHistory[storageCtr[0]] = normOf1DArray(gradCur, primary.nControlsDim, 
						 primary.gradNormChoice);
  
  /* Validate gradient against FD. */
  if (flagCheckGradientFD)
  {
    for (int i = 0; i < primary.nControlsDim; i++)
      gradTemp2[i] = 0.0;
    
    for (int i = 0; i < primary.nObjMC; i++)
    {
      /* The stored generator seed is used. */
      DPCosts -> sampleOnceGradientFD(coefsJkp1Ref, refTableJkp1Ref, 
				      stateKRef, XCur, generatorNoiseFDCheck, 
				      gradTemp);
      
      /* Since the sampleOnce computes the positive J_{k+1} and the
       * gradient, we need to negate the gradient. */
      for (int i = 0; i < primary.nControlsDim; i++)
	gradTemp[i] = -gradTemp[i];
      
      for (int i = 0; i < primary.nControlsDim; i++)
	gradTemp2[i] += gradTemp[i];
    }

    for (int i = 0; i < primary.nControlsDim; i++)
      gradTemp2[i] /= double(primary.nObjMC);

    cout << "% Gradient validation:" << endl;
    cout << "%    Analytic Gradient           FD Gradient" << endl;
    for (int i = 0; i < primary.nControlsDim; i++)
      cout << "% " << setw(23) << gradCur[i] << "  " 
	   << setw(23) << gradTemp2[i] << endl;
    cout << endl;
  }
  
}

void SARM::update_a_k()
{

  /* Update gain sequences. */
  switch (primary.SARMGainSeq)
  {
  case 1:
    /* 1/k. */
    a_k = primary.gainMultiplier / double(nIter[0]);
    break;
    
  case 2:
    /* SPSA gain sequence. */
    a_k = SPSAa / pow(SPSAA + nIter[0] + 1.0, SPSAalpha);
    break;
    
  default:
    cout << "Error: SARMGainSeq " << primary.SARMGainSeq
	 << " not available." << endl;
    exit(1);
  }
    
}



/*********************************************************************
 *********************************************************************/
