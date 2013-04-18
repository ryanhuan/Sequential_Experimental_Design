#include "DPCostsInsideExpectation.h"

/*********************************************************************
 * DPCostsInsideExpectation
 *********************************************************************/

DPCostsInsideExpectation::DPCostsInsideExpectation(Controls const &refControls)
{
  
  /* Load controls, and pointers for external variables. */
  primary = refControls;
  
  /* Allocate memory. */
  stateKp1 = new double[primary.nStatesDim];
  noiseK = new double[primary.nNoiseDim];
  thetaSample = new double[primary.nParamsDim];
  forwardModelOutputs = new double[primary.nNoiseDim];
  noiseStdDev = new double[primary.nNoiseDim];
  z = new double[primary.nNoiseDim];
  zetas = new double[primary.nStatesDim];
  psi = new double*[primary.nStatesDim];
  psi[0] = new double[primary.nStatesDim * (primary.pOrder + 1)];
  dpsi = new double*[primary.nStatesDim];
  dpsi[0] = new double[primary.nStatesDim * (primary.pOrder + 1)];
  for (int i = 1; i < primary.nStatesDim; i++)
  {
    psi[i] = psi[i - 1] + primary.pOrder + 1;
    dpsi[i] = dpsi[i - 1] + primary.pOrder + 1;
  }
  dValuedStateKp1Storage = new double[primary.nStatesDim];
  dStateKp1dControl = new double*[primary.nStatesDim];
  dStateKp1dControl[0] = new double[primary.nStatesDim * primary.nControlsDim];
  for (int i = 1; i < primary.nStatesDim; i++)
    dStateKp1dControl[i] = dStateKp1dControl[i - 1] + primary.nControlsDim;
  
  psiProd = 0.0;
  
}

DPCostsInsideExpectation::~DPCostsInsideExpectation()
{
  
  /* Free memory. */
  delete [] stateKp1;
  delete [] noiseK;
  delete [] thetaSample;
  delete [] forwardModelOutputs;
  delete [] noiseStdDev;
  delete [] z;
  delete [] zetas;
  delete [] psi[0];
  delete [] psi;
  delete [] dpsi[0];
  delete [] dpsi;
  delete [] dValuedStateKp1Storage;
  delete [] dStateKp1dControl[0];
  delete [] dStateKp1dControl;
  
}

void DPCostsInsideExpectation::setNoiseStdDev()
{

  /* Compute the noise standard deviation. */
  for (int i = 0; i < primary.nNoiseDim; i++)
    noiseStdDev[i] = primary.noiseStdDevConstant;

}  

double DPCostsInsideExpectation::sampleOnce(double const * const coefsJkp1Ref, 
					    int const * const * const refTableJkp1Ref, 
					    double const * const stateKRef, 
					    double const * const controlKRef, 
					    gsl_rng* generatorAll)
{

  /* Set noise standard deviation. */
  setNoiseStdDev();

  /* Generate noiseK. */
  thetaSample[0] = stateKRef[0] 
    + sqrt(stateKRef[1]) * gsl_ran_gaussian(generatorAll, 1.0);
  
  /* Evaluate forward model to obtain output. */
  /* F = m * g + e. */
  z[0] = gsl_ran_gaussian(generatorAll, 1.0);
  noiseK[0] = controlKRef[0] * thetaSample[0] + noiseStdDev[0] * z[0];

  /* Evaluate system equation to compute state at stage k+1. */
  evalSystemEquation(stateKRef, controlKRef, noiseK);

  /* Evaluate J_{k+1}. */
  return evalJkp1(coefsJkp1Ref, refTableJkp1Ref, stateKp1)
    - primary.stageCostQuadControlsWeight * controlKRef[0] * controlKRef[0];
      
}


void DPCostsInsideExpectation::evalSystemEquation
(double const * const stateKRef, double const * const controlKRef, 
 double const * const noiseKRef)
{
  
  /* For now, the system equation is Bayes' Theorem on the 2D Gaussian
   * states. */
  
  /* Update mean and standard deviation. Note the second variable of
   * state is the variance, not standard deviation. */
  stateKp1[0] = (noiseKRef[0] * controlKRef[0] * stateKRef[1] + stateKRef[0] 
		 * noiseStdDev[0] * noiseStdDev[0]) 
    / (controlKRef[0] * controlKRef[0] * stateKRef[1] 
       + noiseStdDev[0] * noiseStdDev[0]);
  
  stateKp1[1] = stateKRef[1] * noiseStdDev[0] * noiseStdDev[0] 
    / (controlKRef[0] * controlKRef[0] * stateKRef[1] 
       + noiseStdDev[0] * noiseStdDev[0]);
  
}

double DPCostsInsideExpectation::evalJkp1(double const * const coefsJkp1Ref, 
					  int const * const * const refTableJkp1Ref, 
					  double const * const stateKp1Ref)
{

  double output = 0.5 * (stateKp1Ref[1] / primary.initialState[1]
			 + (stateKp1Ref[0] - primary.initialState[0]) 
			 * (stateKp1Ref[0] - primary.initialState[0]) 
			 / primary.initialState[1]
			 - log(stateKp1Ref[1] / primary.initialState[1]) - 1.0);

  
  // /* Transform the state back to the (normalized) input variables (the
  //  * xis or zetas) to the PC expansion for J_{k+1}. */
  // for (int i = 0; i < primary.nStatesGaussianPCWeight; i++)
  //   zetas[i] = (stateKp1Ref[i] - primary.gaussianPCWeightParams[i * 2]) 
  //     / sqrt(primary.gaussianPCWeightParams[i * 2 + 1]);

  // for (int i = 0; i < primary.nStatesUniformPCWeight; i++)
  //   zetas[i + primary.nStatesGaussianPCWeight] = 
  //     (stateKp1Ref[i + primary.nStatesGaussianPCWeight] 
  //      - primary.uniformPCWeightParams[i * 2]) * 2.0 / 
  //     (primary.uniformPCWeightParams[i * 2 + 1] - primary.uniformPCWeightParams[i * 2])
  //     - 1.0;

  // /* Evaluate the PC expansion (J_{k+1}) at these dimensionless state
  //  * values. */
  // /* Compute the polynomial basis values (psi's) for this zeta. These
  //  * have been verified manually. */
  // for (int i = 0; i < primary.nStatesGaussianPCWeight; i++)
  //   threeTermHermite(zetas[i], primary.pOrder, psi[i]);
  // for (int i = 0; i < primary.nStatesUniformPCWeight; i++)
  //   threeTermLegendre(zetas[i + primary.nStatesGaussianPCWeight], 
  // 		      primary.pOrder, psi[i + primary.nStatesGaussianPCWeight]);
  
  // double output;
  
  // /* Compute coefficient sum. */
  // for (int i = 0; i < primary.nTotalPCETerms; i++)
  // {
    
  //   /* Initialization and product. */
  //   psiProd = 1.0;
    
  //   /* Contribution. */
  //   for (int j = 0; j < primary.nStatesDim; j++)
  //     psiProd *= psi[j][refTableJkp1Ref[i][j]];
    
  //   if (i == 0)
  //   {
  //     /* First term. */
  //     output = coefsJkp1Ref[i] * psiProd;
  //     kC = 0.0;
  //   }
  //   else
  //   {
  //     /* Subsequent terms. */
  //     kY = coefsJkp1Ref[i] * psiProd - kC;
  //     kT = output + kY;
  //     kC = (kT - output) - kY;
  //     output = kT;
  //   }
    
  // }
  
  return output;  
  
}

void DPCostsInsideExpectation::sampleOnceGradient
(double const * const coefsJkp1Ref, 
 int const * const * const refTableJkp1Ref, 
 double const * const stateKRef, double const * const controlKRef,
 double * const grad)
{

  /* Compute the gradient of the value function with respect to the
   * state at k+1. Here we assume the state at k+1 is already computed
   * from the previous evalJkp1(), thus it is not recomputed. */
  dValuedStateKp1(coefsJkp1Ref, refTableJkp1Ref);

  /* Compute the Jacobian of the state at k+1 (i.e., the system
   * equation output) with respect to the control variables. */
  dSystemEquationdControl(stateKRef, controlKRef, noiseK);
  
  /* Construct the gradient by using the total derivative rule. Do not
   * forget the derivative from the stage cost as well. */
  for (int i = 0; i < primary.nControlsDim; i++)
  {
    grad[i] = - 2.0 * primary.stageCostQuadControlsWeight * controlKRef[0];
    for (int j = 0; j < primary.nStatesDim; j++)
      grad[i] += dValuedStateKp1Storage[j] * dStateKp1dControl[j][i];
  }

}

void DPCostsInsideExpectation::dValuedStateKp1
(double const * const coefsJkp1Ref, 
 int const * const * const refTableJkp1Ref)
{

  /* This analytic derivative hass been verified with finite
   * difference. */

  /* The polynomial basis values (psi's) for this zeta are already
   * computed from the preious evalJkp1() call. Make sure this
   * function is ONLY called right after its corresponding evalJkp1()
   * call. */
    
  /* Compute the polynomial basis DERIVATIVE values (dpsi's) for this
   * zeta. These have been verified manually. */
  for (int i = 0; i < primary.nStatesGaussianPCWeight; i++)
    threeTermHermiteDerivatives(zetas[i], primary.pOrder, psi[i], dpsi[i]);
  for (int i = 0; i < primary.nStatesUniformPCWeight; i++)
    threeTermLegendreDerivatives(zetas[i + primary.nStatesGaussianPCWeight], 
				 primary.pOrder, 
				 psi[i + primary.nStatesGaussianPCWeight],
				 dpsi[i + primary.nStatesGaussianPCWeight]);
  
  /* Loop over all components of the state variable (what we are
   * differentiating with respect to here). */
  double dXidState(0.0);
  for (int a = 0; a < primary.nStatesDim; a++)
  {
    
    if (a < primary.nStatesGaussianPCWeight)
      dXidState = 1.0 / sqrt(primary.gaussianPCWeightParams[a * 2 + 1]);
    else if (a >= primary.nStatesGaussianPCWeight)
      dXidState = 2.0 
	/ (primary.uniformPCWeightParams[(a - primary.nStatesGaussianPCWeight) * 2 + 1] 
	   - primary.uniformPCWeightParams[(a - primary.nStatesGaussianPCWeight) * 2]);
    
    /* Compute coefficient sum. */
    for (int i = 0; i < primary.nTotalPCETerms; i++)
    {
      
      /* Initialization and product. */
      psiProd = 1.0;
      
      for (int j = 0; j < primary.nStatesDim; j++)
	if (j == a)
	  psiProd *= dpsi[j][refTableJkp1Ref[i][j]] * dXidState;
	else
	  psiProd *= psi[j][refTableJkp1Ref[i][j]];
      
      if (i == 0)
      {
	/* First term. */
	dValuedStateKp1Storage[a] = coefsJkp1Ref[i] * psiProd;
	kC = 0.0;
      }
      else
      {
	/* Subsequent terms. */
	kY = coefsJkp1Ref[i] * psiProd - kC;
	kT = dValuedStateKp1Storage[a] + kY;
	kC = (kT - dValuedStateKp1Storage[a]) - kY;
	dValuedStateKp1Storage[a] = kT;
      }

    }
    
  }

//   //!!! FD:
//   /* Loop over all components of the state variable (what we are
//    * differentiating with respect to here). */
//   double stateP(0.0), stateM(0.0);
//   double fp(0.0), fm(0.0);
//   double** psiP = new double*[primary.nStatesDim];
//   psiP[0] = new double[primary.nStatesDim * (primary.pOrder + 1)];
//   double** psiM = new double*[primary.nStatesDim];
//   psiM[0] = new double[primary.nStatesDim * (primary.pOrder + 1)];
//   for (int i = 1; i < primary.nStatesDim; i++)
//   {
//     psiP[i] = psiP[i - 1] + primary.pOrder + 1;
//     psiM[i] = psiM[i - 1] + primary.pOrder + 1;
//   }

//   stateP = zetas[1] + sqrt(EPS);
//   stateM = zetas[1] - sqrt(EPS);
// //  cout << "z " << zetas[0] << "  " << sqrt(EPS) << endl;

//   // for (int i = 0; i < primary.nStatesGaussianPCWeight; i++)
//   //   threeTermHermite(stateP, primary.pOrder, psiP[i]);
//   for (int i = 0; i < primary.nStatesGaussianPCWeight; i++)
//     threeTermHermite(zetas[i], primary.pOrder, psiP[i]);
//   for (int i = 0; i < primary.nStatesUniformPCWeight; i++)
//     threeTermLegendre(stateP, 
// 		      primary.pOrder, psiP[i + primary.nStatesGaussianPCWeight]);
//   // for (int i = 0; i < primary.nStatesUniformPCWeight; i++)
//   //   threeTermLegendre(zetas[i + primary.nStatesGaussianPCWeight], 
//   // 		      primary.pOrder, psiP[i + primary.nStatesGaussianPCWeight]);

//   // for (int i = 0; i < primary.nStatesGaussianPCWeight; i++)
//   //   threeTermHermite(stateM, primary.pOrder, psiM[i]);
//   for (int i = 0; i < primary.nStatesGaussianPCWeight; i++)
//     threeTermHermite(zetas[i], primary.pOrder, psiM[i]);
//   for (int i = 0; i < primary.nStatesUniformPCWeight; i++)
//     threeTermLegendre(stateM, 
// 		      primary.pOrder, psiM[i + primary.nStatesGaussianPCWeight]);
//   // for (int i = 0; i < primary.nStatesUniformPCWeight; i++)
//   //   threeTermLegendre(zetas[i + primary.nStatesGaussianPCWeight], 
//   // 		      primary.pOrder, psiM[i + primary.nStatesGaussianPCWeight]);
    
//     /* Compute coefficient sum. */
//     for (int i = 0; i < primary.nTotalPCETerms; i++)
//     {
      
//       /* Initialization and product. */
//       psiProd = 1.0;
      
//       for (int j = 0; j < primary.nStatesDim; j++)
// 	psiProd *= psiP[j][refTableJkp1Ref[i][j]];
      
//       if (i == 0)
//       {
// 	/* First term. */
// 	fp = coefsJkp1Ref[i] * psiProd;
// 	kC = 0.0;
//       }
//       else
//       {
// 	/* Subsequent terms. */
// 	kY = coefsJkp1Ref[i] * psiProd - kC;
// 	kT = fp + kY;
// 	kC = (kT - fp) - kY;
// 	fp = kT;
//       }

//     }

//     /* Compute coefficient sum. */
//     for (int i = 0; i < primary.nTotalPCETerms; i++)
//     {
      
//       /* Initialization and product. */
//       psiProd = 1.0;
      
//       for (int j = 0; j < primary.nStatesDim; j++)
// 	psiProd *= psiM[j][refTableJkp1Ref[i][j]];
      
//       if (i == 0)
//       {
// 	/* First term. */
// 	fm = coefsJkp1Ref[i] * psiProd;
// 	kC = 0.0;
//       }
//       else
//       {
// 	/* Subsequent terms. */
// 	kY = coefsJkp1Ref[i] * psiProd - kC;
// 	kT = fm + kY;
// 	kC = (kT - fm) - kY;
// 	fm = kT;
//       }

//     }
//     double fd = (fp-fm) / (2.0 * sqrt(EPS));
//     cout << dValuedStateKp1Storage[1] << "  " << fd << endl;
  
}

void DPCostsInsideExpectation::dSystemEquationdControl
(double const * const stateKRef, double const * const controlKRef, 
 double const * const noiseKRef)
{
  
  /* For now, the system equation is Bayes' Theorem on the 2D Gaussian
   * states. We take the derivative of each component with respect to
   * each component of the control. This Jacobian is verified with
   * finite difference. */
  
  /* Compute Jacobian. */
  dStateKp1dControl[0][0] = noiseKRef[0] * stateKRef[1]
    / (controlKRef[0] * controlKRef[0] * stateKRef[1] 
       + noiseStdDev[0] * noiseStdDev[0])
    - (noiseKRef[0] * controlKRef[0] * stateKRef[1] + stateKRef[0] 
       * noiseStdDev[0] * noiseStdDev[0]) * 2.0 * controlKRef[0] * stateKRef[1]
    / ((controlKRef[0] * controlKRef[0] * stateKRef[1] 
	+ noiseStdDev[0] * noiseStdDev[0]) 
       * (controlKRef[0] * controlKRef[0] * stateKRef[1] 
	  + noiseStdDev[0] * noiseStdDev[0]))
    + controlKRef[0] * stateKRef[1] * thetaSample[0] 
    / (controlKRef[0] * controlKRef[0] * stateKRef[1] 
       + noiseStdDev[0] * noiseStdDev[0]);
  
  dStateKp1dControl[1][0] = -(stateKRef[1] * noiseStdDev[0] * noiseStdDev[0])
    * 2.0 * controlKRef[0] * stateKRef[1]
    / ((controlKRef[0] * controlKRef[0] * stateKRef[1] 
	+ noiseStdDev[0] * noiseStdDev[0]) 
       * (controlKRef[0] * controlKRef[0] * stateKRef[1] 
	  + noiseStdDev[0] * noiseStdDev[0]));

  // //!!! FD
  // double* stateKp1P = new double[primary.nStatesDim];
  // double* stateKp1M = new double[primary.nStatesDim];
  // double controlP(0.0), controlM(0.0);
  // double* noiseKP = new double[primary.nNoiseDim];
  // double* noiseKM = new double[primary.nNoiseDim];
  // controlP = controlKRef[0] + sqrt(EPS);
  // controlM = controlKRef[0] - sqrt(EPS);
  
  // noiseKP[0] = controlP * thetaSample[0] + noiseStdDev[0] * z[0];
  // noiseKM[0] = controlM * thetaSample[0] + noiseStdDev[0] * z[0];

  // stateKp1P[0] = (noiseKP[0] * controlP * stateKRef[1] + stateKRef[0] 
  // 		 * noiseStdDev[0] * noiseStdDev[0]) 
  //   / (controlP * controlP * stateKRef[1] 
  //      + noiseStdDev[0] * noiseStdDev[0]);
  
  // stateKp1P[1] = stateKRef[1] * noiseStdDev[0] * noiseStdDev[0] 
  //   / (controlP * controlP * stateKRef[1] 
  //      + noiseStdDev[0] * noiseStdDev[0]);


  // stateKp1M[0] = (noiseKM[0] * controlM * stateKRef[1] + stateKRef[0] 
  // 		 * noiseStdDev[0] * noiseStdDev[0]) 
  //   / (controlM * controlM * stateKRef[1] 
  //      + noiseStdDev[0] * noiseStdDev[0]);
  
  // stateKp1M[1] = stateKRef[1] * noiseStdDev[0] * noiseStdDev[0] 
  //   / (controlM * controlM * stateKRef[1] 
  //      + noiseStdDev[0] * noiseStdDev[0]);

  // double fd1 = (stateKp1P[0]-stateKp1M[0]) / (2.0 * sqrt(EPS));
  // double fd2 = (stateKp1P[1]-stateKp1M[1]) / (2.0 * sqrt(EPS));
  // cout << dStateKp1dControl[0][0] << "  " << fd1 << endl;
  // cout << dStateKp1dControl[1][0] << "  " << fd2 << endl;
  
}

void DPCostsInsideExpectation::sampleOnceGradientFD
(double const * const coefsJkp1Ref, 
 int const * const * const refTableJkp1Ref, 
 double const * const stateKRef, 
 double const * const controlKRef, 
 gsl_rng* generator, double * const grad)
{
  
  /* Note we assume the generator is already initialized with the
   * appropriate seed when passed into this function. */
  static gsl_rng_type const * rngType = gsl_rng_ranlxs0;
  static gsl_rng* generatorRef = gsl_rng_alloc(rngType);
  gsl_rng_memcpy(generatorRef, generator);
  
  /* Initializations. */
  double perturb(0.0), fp(0.0), fm(0.0);
  static double* controlKPert = new double[primary.nControlsDim];
  
  for (int i = 0; i < primary.nControlsDim; i++)
    controlKPert[i] = controlKRef[i];
  
  for (int i = 0; i < primary.nControlsDim; i++)
  {
    /* Reset previous perturbation. */
    if (i > 0)
      controlKPert[i - 1] = controlKRef[i - 1];

    /* Positive perturbation. */
    perturb = max<double>(fabs(controlKRef[i]), 1.0) * sqrt(EPS);
    controlKPert[i] = controlKRef[i] + perturb;
    fp = 0.0;
    gsl_rng_memcpy(generator, generatorRef);
    fp = sampleOnce(coefsJkp1Ref, refTableJkp1Ref, 
		    stateKRef, controlKPert, generator);

    /* Negative perturbation. */
    controlKPert[i] = controlKRef[i] - perturb;
    fm = 0.0;
    gsl_rng_memcpy(generator, generatorRef);
    fm = sampleOnce(coefsJkp1Ref, refTableJkp1Ref, 
		    stateKRef, controlKPert, generator);
    
    /* Compute gradient component. */
    grad[i] = (fp - fm) / (2.0 * perturb);
  }
  
}
