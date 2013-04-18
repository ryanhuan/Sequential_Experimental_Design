#include "userFunctions.h"

double systemEquation(Controls const &primary, double const * const state, 
		      double const * const control, double const * const disturbance,
		      double const * const newState)
{
  
  /* Bayesian inference for conjugate linear Gaussian model. */
  newState[0] = (disturbance[0] / control[0]) / (primary.!!!
  
}

double maxExpectation(Controls const &primary, GenericInputs &allInputs, 
		      double const * const state)
{
  
  /* Initialization. */
  double value(0.0);

  /* Initialize GSL random generator. */
  static gsl_rng_type const * rngType(gsl_rng_ranlxs0);
  static gsl_rng* generatorInit = gsl_rng_alloc(rngType);
  gsl_rng_env_setup();
  gsl_rng_set(generatorInit, rand() + primary.rank);
  
  /* Initialization of StochasticSearch (mostly just memory
   * allocation); does not include making the initial position at this
   * point. */
  static StochasticSearch stoch(primary);
  
  double tempValue(0.0);
  
  /* Multiple optimizations and take the maximum high-quality
   * objective value. */
  for (int z = 0; z < primary.nStochOptPerStateK; z++)
  {
    
    /* Note the X here in optimization algorithms are the design
     * variables (i.e., the control, not the state). The state is
     * fixed for an entire optimization run. The location of the state
     * is already stored in the StochasticSearch class when it is
     * initialized below, so all state changes follow
     * automatically. */

    if (primary.randomizeXInitial == 1)
    {
      if (primary.nControlsDim == 1)
      {
	/* Equi-spaced initial points in 1D design space only. */
	if (primary.nStochOptPerStateK == 1)
	  primary.XInitial[0] = (primary.controlsRightBound[0] +
				 + primary.controlsLeftBound[0]) / 2;
	else
	  primary.XInitial[0] = double(z) / double(primary.nStochOptPerStateK - 1) 
	    * (primary.controlsRightBound[0] - primary.controlsLeftBound[0])
	    + primary.controlsLeftBound[0];
      }
      else
	/* Generate random initial position for multi-D design
	 * space. */
	for (int i = 0; i < primary.nControlsDim; i++)
	  primary.XInitial[i] = gsl_rng_uniform(generatorInit)
	    * (primary.controlsRightBound[i] - primary.controlsLeftBound[i])
	    + primary.controlsLeftBound[i];
    }
    
    /* Initialize StochSearch class. */
    stoch.initialPosition(allInputs.innerFcn, state);
    
    /* Optimize. */
    stoch.search(allInputs.innerFcn, state);
    
    /* Extract final stochastic optimization result. Note this value
     * is negative estimated value of J_k, because in optimization we
     * minimize the negative. */
    tempValue = -stoch.exportFinalObjHighQualityMC();

    /* Take the maximum. If want to see distribution, need to store
     * them or output them. */
    if (z == 0)
      value = tempValue;
    else
      value = max<double>(value, tempValue);
  }
    
  return value;
  
}

double evalTerminalRewardFromKm1(Controls const &primary, 
				 double const * const stateKm1, 
				 double const * const controlKm1, 
				 double const * const disturbanceKm1)
{
  
  /* Initialization. */
  double value(0.0);

  /* !!!. */
  
  /* Evaluate the KL divergence between the final state and the
   * initial prior. */
  if (primary.nStatesDim == 2)
    value = 0.5 * (state[1] / primary.initialState[1]
		   + (state[0] - primary.initialState[0]) 
		   * (state[0] - primary.initialState[0]) 
		   / primary.initialState[1]
		   - log(state[1] / primary.initialState[1]) - 1.0);
  else
  {
    cout << "Error: KL divergence evaluation for multi-dimensional PARAMETER space not yet implemented. " << endl;
    exit(1);
  }
  
  return value;
  
}

double evalTerminalRewardFromKm1(Controls const &primary, 
				 double const * const stateKm1, 
				 double const * const controlKm1, 
				 double const * const disturbanceKm1)
{
  
  /* Initialization. */
  double value(0.0);


  
  /* Evaluate the KL divergence between the final state and the
   * initial prior. */
  if (primary.nStatesDim == 2)
    value = 0.5 * (state[1] / primary.initialState[1]
		   + (state[0] - primary.initialState[0]) 
		   * (state[0] - primary.initialState[0]) 
		   / primary.initialState[1]
		   - log(state[1] / primary.initialState[1]) - 1.0);
  else
  {
    cout << "Error: KL divergence evaluation for multi-dimensional PARAMETER space not yet implemented. " << endl;
    exit(1);
  }
  
  return value;
  
}

double evalJk(Controls const &primary, double const * const stateKRef,
	      int const currentStage, double const * const coefsJkp1Ref, 
	      int const * const * const refTableJkp1Ref)
{
  
  /* Initialization. */
  double value(0.0);
//!!!
  // /* Initialize GSL random generator. */
  // static gsl_rng_type const * rngType(gsl_rng_ranlxs0);
  // static gsl_rng* generatorInit = gsl_rng_alloc(rngType);
  
  // /* Initialization of StochasticSearch (mostly just memory
  //  * allocation); does not include making the initial position at this
  //  * point. */
  // static StochasticSearch stoch(primary);
  
  // // if (currentStage == primary.nStages)
  // // {
  // //   /* Evaluate the KL divergence between the final state and the
  // //    * initial prior. */
  // //   if (primary.nStatesDim == 2)
  // //     value = 0.5 * (stateKRef[1] / primary.initialState[1]
  // // 		     + (stateKRef[0] - primary.initialState[0]) 
  // // 		     * (stateKRef[0] - primary.initialState[0]) 
  // // 		     / primary.initialState[1]
  // // 		     - log(stateKRef[1] / primary.initialState[1]) - 1.0);
  // //   else
  // //   {
  // //     cout << "Error: KL divergence evaluation for multi-dimensional PARAMETER space not yet implemented. " << endl;
  // //     exit(1);
  // //   }
  // // }
  // // else
  // // {
  //   /* Else we will compute a Monte Carlo estimate to the Jk
  //    * function. */
    
  //   /* Initialize GSL random generator. */
  //   gsl_rng_env_setup();
  //   gsl_rng_set(generatorInit, rand() + primary.rank);

  //   double tempValue(0.0);
    
  //   /* Multiple optimizations and take the maximum high-quality
  //    * objective value. */
  //   for (int z = 0; z < primary.nStochOptPerStateK; z++)
  //   {
    
  //     /* Note the X here in optimization algorithms are the design
  //      * variables (i.e., the control, not the state). The state is
  //      * fixed for an entire optimization run. The location of the state
  //      * is already stored in the StochasticSearch class when it was
  //      * first initialized, so all state changes follow
  //      * automatically. */

  //     if (primary.randomizeXInitial == 1)
  //     {
  // 	if (primary.nControlsDim == 1)
  // 	{
  // 	  /* Equi-spaced initial points in 1D design space only. */
  // 	  if (primary.nStochOptPerStateK == 1)
  // 	    primary.XInitial[0] = (primary.controlsRightBound[0] +
  // 	      + primary.controlsLeftBound[0]) / 2;
  // 	  else
  // 	    primary.XInitial[0] = double(z) / double(primary.nStochOptPerStateK - 1) 
  // 	      * (primary.controlsRightBound[0] - primary.controlsLeftBound[0])
  // 	      + primary.controlsLeftBound[0];
  // 	}
  // 	else
  // 	  /* Generate random initial position for multi-D design
  // 	   * space. */
  // 	  for (int i = 0; i < primary.nControlsDim; i++)
  // 	    primary.XInitial[i] = gsl_rng_uniform(generatorInit)
  // 	      * (primary.controlsRightBound[i] - primary.controlsLeftBound[i])
  // 	      + primary.controlsLeftBound[i];
  //     }
    
  //     /* Initialize StochSearch class. */
  //     stoch.initialPosition(coefsJkp1Ref, refTableJkp1Ref, stateKRef);
    
  //     /* Optimize. */
  //     stoch.search(coefsJkp1Ref, refTableJkp1Ref, stateKRef);
      
  //     /* Extract final stochastic optimization result. Note this value
  //      * is negative estimated value of J_k, because in optimization
  //      * we minimize the negative. */
  //     tempValue = -stoch.exportFinalObjHighQualityMC();

  //     /* Take the maximum. If want to see distribution, need to store
  //      * them or output them. */
  //     if (z == 0)
  // 	value = tempValue;
  //     else
  // 	value = max<double>(value, tempValue);
  //   }
  // // }
  
  return value;
  
}
