#include "userFunctions.h"

void setNoiseStdDev(InputParams const &algParams, vector<double> &noiseStdDev)
{
  
  /* Set the Gaussian noise standard deviation. */
  switch (algParams.noiseModel)
  {
  case 1:
    /* Constant standard deviation. */
    for (unsigned int i = 0; i < noiseStdDev.size(); i++)
      noiseStdDev[i] = 1.0;

    break;
    
  default:
    cout << "Error: Noise model " << algParams.noiseModel
	 << " not available." << endl;
    exit(1);
  }
  
}

void forwardModel(InputParams const &algParams, vector<double> const &theta, 
		  vector<double> const &state, vector<double> const &control, 
		  vector<double> &modelOutputs)
{

  /* Evaluate forward model (no noise). */
  switch (algParams.forwardModel)
  {
  case 1:
    /* G = d * theta. */
    if (modelOutputs.size() != 1)
    {
      cout << "Error: forwardModel " << algParams.forwardModel << " for 1D model output only. " 
	   << endl;
      exit(1);
    }
    modelOutputs[0] = control[0] * theta[0];

    break;
    
  default:
    cout << "Error: Forward model " << algParams.forwardModel
	 << " not available." << endl;
    exit(1);
  }

}

void generateDisturbance(InputParams const &algParams, vector<double> const &state, 
			 vector<double> const &control, gsl_rng* generator,
			 vector<double> &disturbance)
{

  /* Generate model parameter sample. */
  //!!! in future generalize with a sample-generating subrouting
  //!!! depending on the state we choose
  static vector<double> thetaSample(algParams.nInferenceParamsDim, 0.0);
  for (unsigned int i = 0; i < thetaSample.size(); i++)
    thetaSample[i] = state[0] + sqrt(state[1]) * gsl_ran_gaussian(generator, 1.0);
  
  /* Set noise standard deviation and evaluate forward model. */
  static vector<double> noiseStdDev(algParams.nDisturbanceDim, 0.0);
  setNoiseStdDev(algParams, noiseStdDev);
  static vector<double> modelOutputs(algParams.nDisturbanceDim, 0.0);
  forwardModel(algParams, thetaSample, state, control, modelOutputs);
  
  /* Assume additive noise. */
  for (unsigned int i = 0; i < disturbance.size(); i++)
    disturbance[i] = modelOutputs[i] 
      + noiseStdDev[i] * gsl_ran_gaussian(generator, 1.0);
  
}

void systemEquation(InputParams const &algParams, vector<double> const &state, 
		    vector<double> const &control, vector<double> const &disturbance,
		    vector<double> &newState)
{
  
  /* Compute the Gaussian noise standard deviation. */
  static vector<double> noiseStdDev(algParams.nDisturbanceDim, 0.0);
  setNoiseStdDev(algParams, noiseStdDev);
  
  //!!!Bayesian inference for conjugate 1D linear Gaussian
  //!!!model. Generalization needed for future for different state choices.
  newState[0] = (disturbance[0] * control[0] * state[1] + state[0] 
		 * noiseStdDev[0] * noiseStdDev[0]) 
    / (control[0] * control[0] * state[1] + noiseStdDev[0] * noiseStdDev[0]);
  newState[1] = state[1] * noiseStdDev[0] * noiseStdDev[0] 
    / (control[0] * control[0] * state[1] + noiseStdDev[0] * noiseStdDev[0]);
  
}

double stageCost(InputParams const &algParams, vector<double> const &state, 
		 vector<double> const &control, vector<double> const &disturbance)
{
  
  double cost(0.0);
  
  /* Compute the stage cost. */
  //!!! quadratic control-cost for now.
  for (unsigned int i = 0; i < control.size(); i++)
    cost += 0.01 * control[i] * control[i];
  
  return cost;

}

double maxExpectation(InputParams const &algParams, GenericInputs const &allInputs, 
		      vector<double> const &state)
{
  
  static vector<double> dummy(algParams.nControlsDim, 0.0);
  return maxExpectation(algParams, allInputs, state, dummy);
  
}

double maxExpectation(InputParams const &algParams, GenericInputs const &allInputs, 
		      vector<double> const &state, 
		      vector<double> &finalPositionExternal)
{
  
  /* Initialization. */
  double value(0.0), tempValue(0.0);
  static vector<double> tempPosition(finalPositionExternal.size(), 0.0);
  
  /* Initialize GSL random generator. */
  static gsl_rng_type const * rngType(gsl_rng_ranlxs0);
  static gsl_rng* generatorInit = gsl_rng_alloc(rngType);
  gsl_rng_env_setup();
  gsl_rng_set(generatorInit, rand() + algParams.rank);
  
  /* Initialization of StochasticSearch (mostly just memory
   * allocation); does not include making the initial position at this
   * point. */
  static StochasticSearch stoch(algParams);
  static vector<double> XInitial(algParams.nSearchDim, 0.0);
  
  /* Multiple optimizations and take the maximum high-quality
   * objective value. */
  for (int z = 0; z < algParams.nStochOptPerState; z++)
  {
    
    /* Note the X here in optimization algorithms are the design
     * variables (i.e., the control, not the state). The state is
     * fixed for an entire optimization run. The location of the state
     * is already stored in the StochasticSearch class when it is
     * initialized below, so all state changes follow
     * automatically. */

    if (algParams.randomizeXInitial == 1)
    {
      if (algParams.nSearchDim == 1)
      {
  	/* Equi-spaced initial points in 1D search space only. */
  	if (algParams.nStochOptPerState == 1)
  	  XInitial[0] = (algParams.controlsUpperBounds[0] +
			 + algParams.controlsLowerBounds[0]) / 2;
  	else
  	  XInitial[0] = double(z) / double(algParams.nStochOptPerState - 1) 
  	    * (algParams.controlsUpperBounds[0] - algParams.controlsLowerBounds[0])
  	    + algParams.controlsLowerBounds[0];
      }
      else
  	/* Generate random initial position for multi-D design
  	 * space. */
  	for (unsigned int i = 0; i < XInitial.size(); i++)
  	  XInitial[i] = gsl_rng_uniform(generatorInit)
  	    * (algParams.controlsUpperBounds[i] - algParams.controlsLowerBounds[i])
  	    + algParams.controlsLowerBounds[i];
    }
    else
      for (unsigned int i = 0; i < XInitial.size(); i++)
	XInitial[i] = algParams.userXInitial[i];
    
    /* Initialize StochSearch class. */
    stoch.initializePosition(XInitial);
    
    /* Optimize. */
    stoch.search(state, allInputs);
    
    /* Extract final stochastic optimization result. Note this value
     * is negative estimated value of J_k, because in optimization we
     * minimize the negative. */
    tempValue = -stoch.exportFinalObjHighQualityMC();

    /* Extract final position. */
    stoch.exportFinalPosition(tempPosition);

    /* Take the maximum. If want to see distribution, need to store
     * them or output them. */
    if (z == 0)
    {
      value = tempValue;
      finalPositionExternal = tempPosition;
    }
    else
    {
      if (tempValue > value)
      {
	value = tempValue;
	finalPositionExternal = tempPosition;
      }
    }
  }
    
  return value;
  
}

double evalTerminalReward(InputParams const &algParams, vector<double> const &state)
{
  
  /* Initialization. */
  double value(0.0);

  /* Evaluate the KL divergence between the final state and the
   * initial prior. */
  //!!! for now we will use the analytical formula for the KL between
  //!!! two 1D Gaussians. Generalize in future.
  if (algParams.nStatesDim == 2)
    value = 0.5 * (state[1] / algParams.initialState[1]
		   + (state[0] - algParams.initialState[0]) 
		   * (state[0] - algParams.initialState[0]) 
		   / algParams.initialState[1]
		   - log(state[1] / algParams.initialState[1]) - 1.0);
  else
  {
    cout << "Error: KL divergence evaluation for multi-dimensional PARAMETER space not yet implemented. " << endl;
    exit(1);
  }
  
  return value;
  
}

// double evalJk(Controls const &algParams, double const * const stateKRef,
// 	      int const currentStage, double const * const coefsJkp1Ref, 
// 	      int const * const * const refTableJkp1Ref)
// {
  
//   /* Initialization. */
//   double value(0.0);
// //!!!
//   // /* Initialize GSL random generator. */
//   // static gsl_rng_type const * rngType(gsl_rng_ranlxs0);
//   // static gsl_rng* generatorInit = gsl_rng_alloc(rngType);
  
//   // /* Initialization of StochasticSearch (mostly just memory
//   //  * allocation); does not include making the initial position at this
//   //  * point. */
//   // static StochasticSearch stoch(algParams);
  
//   // // if (currentStage == algParams.nStages)
//   // // {
//   // //   /* Evaluate the KL divergence between the final state and the
//   // //    * initial prior. */
//   // //   if (algParams.nStatesDim == 2)
//   // //     value = 0.5 * (stateKRef[1] / algParams.initialState[1]
//   // // 		     + (stateKRef[0] - algParams.initialState[0]) 
//   // // 		     * (stateKRef[0] - algParams.initialState[0]) 
//   // // 		     / algParams.initialState[1]
//   // // 		     - log(stateKRef[1] / algParams.initialState[1]) - 1.0);
//   // //   else
//   // //   {
//   // //     cout << "Error: KL divergence evaluation for multi-dimensional PARAMETER space not yet implemented. " << endl;
//   // //     exit(1);
//   // //   }
//   // // }
//   // // else
//   // // {
//   //   /* Else we will compute a Monte Carlo estimate to the Jk
//   //    * function. */
    
//   //   /* Initialize GSL random generator. */
//   //   gsl_rng_env_setup();
//   //   gsl_rng_set(generatorInit, rand() + algParams.rank);

//   //   double tempValue(0.0);
    
//   //   /* Multiple optimizations and take the maximum high-quality
//   //    * objective value. */
//   //   for (int z = 0; z < algParams.nStochOptPerStateK; z++)
//   //   {
    
//   //     /* Note the X here in optimization algorithms are the design
//   //      * variables (i.e., the control, not the state). The state is
//   //      * fixed for an entire optimization run. The location of the state
//   //      * is already stored in the StochasticSearch class when it was
//   //      * first initialized, so all state changes follow
//   //      * automatically. */

//   //     if (algParams.randomizeXInitial == 1)
//   //     {
//   // 	if (algParams.nControlsDim == 1)
//   // 	{
//   // 	  /* Equi-spaced initial points in 1D design space only. */
//   // 	  if (algParams.nStochOptPerStateK == 1)
//   // 	    algParams.XInitial[0] = (algParams.controlsRightBound[0] +
//   // 	      + algParams.controlsLeftBound[0]) / 2;
//   // 	  else
//   // 	    algParams.XInitial[0] = double(z) / double(algParams.nStochOptPerStateK - 1) 
//   // 	      * (algParams.controlsRightBound[0] - algParams.controlsLeftBound[0])
//   // 	      + algParams.controlsLeftBound[0];
//   // 	}
//   // 	else
//   // 	  /* Generate random initial position for multi-D design
//   // 	   * space. */
//   // 	  for (int i = 0; i < algParams.nControlsDim; i++)
//   // 	    algParams.XInitial[i] = gsl_rng_uniform(generatorInit)
//   // 	      * (algParams.controlsRightBound[i] - algParams.controlsLeftBound[i])
//   // 	      + algParams.controlsLeftBound[i];
//   //     }
    
//   //     /* Initialize StochSearch class. */
//   //     stoch.initialPosition(coefsJkp1Ref, refTableJkp1Ref, stateKRef);
    
//   //     /* Optimize. */
//   //     stoch.search(coefsJkp1Ref, refTableJkp1Ref, stateKRef);
      
//   //     /* Extract final stochastic optimization result. Note this value
//   //      * is negative estimated value of J_k, because in optimization
//   //      * we minimize the negative. */
//   //     tempValue = -stoch.exportFinalObjHighQualityMC();

//   //     /* Take the maximum. If want to see distribution, need to store
//   //      * them or output them. */
//   //     if (z == 0)
//   // 	value = tempValue;
//   //     else
//   // 	value = max<double>(value, tempValue);
//   //   }
//   // // }
  
//   return value;
  
// }
