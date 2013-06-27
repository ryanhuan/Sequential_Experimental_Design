#include "stochasticSearch.h"

/*********************************************************************
 * StochasticSearch
 *********************************************************************/

StochasticSearch::StochasticSearch()
{
  initialized = 0;
}

StochasticSearch::StochasticSearch(InputParams const &algParamsRef)
{
  initialize(algParamsRef);
}

StochasticSearch::StochasticSearch(StochasticSearch const &other)
{
  
  /* Initializations. */
  algParams = other.algParams;

  nIter = other.nIter;
  storageCtr = other.storageCtr;
  XInitial = other.XInitial;
  dataLabels = other.dataLabels;
  terminationLabel = other.terminationLabel;
  finalObjHighQualityMC = other.finalObjHighQualityMC;
  iterHistory = other.iterHistory;
  XHistory = other.XHistory;
  valHistory = other.valHistory;
  gradNormHistory = other.gradNormHistory;
  terminationCondition = other.terminationCondition;
  nConsecRelXNorm = other.nConsecRelXNorm;
  previousXNorm = other.previousXNorm;
  curXNorm = other.curXNorm;
  isLocalMin = other.isLocalMin;
  
  //!!!
  // DPCostsForHighQualityMC = new DPCostsInsideExpectation(refControls);

  /* Random number generator initializations. */
  rngType = gsl_rng_ranlxs0;
  generatorNoiseForHighQualityMC = gsl_rng_alloc(rngType);
  gsl_rng_env_setup();
  gsl_rng_set(generatorNoiseForHighQualityMC, rand() + algParams.rank);
  
  initialized = other.initialized;
  positionInitialized = other.positionInitialized;

}

StochasticSearch::~StochasticSearch()
{
  
  /* Free memory. */
  if (initialized)
  {
    gsl_rng_free(generatorNoiseForHighQualityMC);
  }
  
}

StochasticSearch& StochasticSearch::operator=(StochasticSearch const &rhs)
{

  /* Protect against invalid self-assignment. */
  if (this != &rhs)
  {
    
    /* Initializations. */
    algParams = rhs.algParams;

    nIter = rhs.nIter;
    storageCtr = rhs.storageCtr;
    XInitial = rhs.XInitial;
    dataLabels = rhs.dataLabels;
    terminationLabel = rhs.terminationLabel;
    finalObjHighQualityMC = rhs.finalObjHighQualityMC;
    iterHistory = rhs.iterHistory;
    XHistory = rhs.XHistory;
    valHistory = rhs.valHistory;
    gradNormHistory = rhs.gradNormHistory;
    terminationCondition = rhs.terminationCondition;
    nConsecRelXNorm = rhs.nConsecRelXNorm;
    previousXNorm = rhs.previousXNorm;
    curXNorm = rhs.curXNorm;
    isLocalMin = rhs.isLocalMin;
  
    //!!!
    // DPCostsForHighQualityMC = new DPCostsInsideExpectation(refControls);

    /* Random number generator initializations. */
    rngType = gsl_rng_ranlxs0;
    generatorNoiseForHighQualityMC = gsl_rng_alloc(rngType);
    gsl_rng_env_setup();
    gsl_rng_set(generatorNoiseForHighQualityMC, rand() + algParams.rank);
  
    initialized = rhs.initialized;
    positionInitialized = rhs.positionInitialized;

  }
  
  /* By convention, always return *this. */
  return *this;
  
}

void StochasticSearch::initialize(InputParams const &algParamsRef)
{

  /* Load input parameters. */
  algParams = algParamsRef;
  
  /* Allocate memory. */
  iterHistory.assign(algParams.maxOptIters + 1, 0);
  XHistory.assign(algParams.maxOptIters + 1, vector<double>(algParams.nSearchDim, 0.0));
  valHistory.assign(algParams.maxOptIters + 1, 0.0);
  XInitial.assign(algParams.nSearchDim, 0.0);

  //!!!
  // DPCostsForHighQualityMC = new DPCostsInsideExpectation(refControls);

  /* Initialize specific search algorithm class. */
  switch (algParams.optMethod)
  {
  case 6:
    /* SARM. */
    gradNormHistory.assign(algParams.maxOptIters + 1, 0.0);
    // robbinsMonro.initialize(primary, nIter, storageCtr, iterHistory, XHistory,
    // 			    valHistory, gradNormHistory);
    break;
    
  default:
    cout << "Error: Optimization method " << algParams.optMethod
	 << " not available." << endl;
    exit(1);
  }

  /* Initializations. Note that we are not initializing position at
   * this point. Do that separately. */
  makeLabels();

  /* Random number generator initializations. */
  rngType = gsl_rng_ranlxs0;
  generatorNoiseForHighQualityMC = gsl_rng_alloc(rngType);
  gsl_rng_env_setup();
  gsl_rng_set(generatorNoiseForHighQualityMC, rand() + algParams.rank);

  initialized = 1;

}

void StochasticSearch::makeLabels()
{
  
  /* Construct progress and text file label. */
  if (algParams.optMethod == 6)
    /* SARM. */
    dataLabels = string("% Iter,    Gradient Norm,       Coordinates [")
      + num2string<int>(algParams.nSearchDim) + string(" columns]");
  else
  {
    cout << "Error: Optimization method " << algParams.optMethod
	 << " not available." << endl;
    exit(1);
  }
  
}

void StochasticSearch::initializePosition(vector<double> const &XInitialRef)
{

  /* Initialize if not initialized. */
  if (!initialized)
  {
    cout << "Attempting to use StochasticSearch without initializing, please check code. " 
	 << endl;
    exit(1);
  }
  
  /* Initializations. */
  for (unsigned int i = 0; i < XInitial.size(); i++)
    XInitial[i] = XInitialRef[i];
  nIter = 0;
  storageCtr = 0;
  nConsecRelXNorm = 0;
  isLocalMin = 0;  
  iterHistory.assign(iterHistory.size(), 0.0);
  
  /* Call additional initialization functions that are specific to the
   * selected algorithm. For example, gradient estimate if a
   * gradient-based algorithm is used. */
  switch (algParams.optMethod)
  {
  case 6:
    /* SARM. */
//    robbinsMonro -> initialPosition(innerFcn);
    break;
  default:
    cout << "Error: Optimization method " << algParams.optMethod
	 << " not available." << endl;
    exit(1);
  }
  
  /* Set output format and display label. */
  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(16);
  
  /* Display initial conditions. */
  if (algParams.displayOptProgress)
  {
    cout << dataLabels << endl;
    cout << setw(4) << iterHistory[storageCtr] << "  ";
    cout << setw(23) << gradNormHistory[storageCtr] << "  ";
    for (int i = 0; i < algParams.nSearchDim; i++)
      cout << setw(23) << XHistory[storageCtr][i] << "  ";
    cout << endl;
  }

  positionInitialized = 1;
  
}

void StochasticSearch::search(vector<double> const &state, GenericInputs allInputs)
{

  if (!positionInitialized)
  {
    cout << "Attempting to use StochasticSearch without initializing, please check code. " 
	 << endl;
    exit(1);
  }

//   // /* Proceed only if at least 1 iteration is requested. */
//   // if (primary.maxOptIters > 0)
//   // {

//   //   /* Initializations. */
//   //   previousXNorm = normOf1DArray(XHistory[0], primary.nSearchDim, 
//   // 				  primary.relXNormTerminateNormChoice);
    
//   //   do
//   //   {
      
//   //     /* If the user wants to store and display progress every n
//   //      * iterations, then a for-loop should be placed inside the
//   //      * switch statement, wrapping the iterations. */
//   //     switch (primary.optMethod)
//   //     {
//   //     case 6:
//   // 	/* SARM. */
//   // 	robbinsMonro -> iterateOnce(coefsJkp1Ref, refTableJkp1Ref, stateKRef);
//   // 	break;
//   //     default:
//   // 	cout << "Error: Optimization method " << primary.optMethod
//   // 	     << " not available." << endl;
//   // 	exit(1);
//   //     }

//   //     /* Display progress. If local minimum is reached, the next point
//   //      * is not computed and thus should not be displayed (it would
//   //      * display gibberish). */
//   //     if (primary.displayOptProgress && !isLocalMin)
//   //     {
//   // 	cout << setw(4) << iterHistory[storageCtr[0]] << "  ";
//   // 	cout << setw(23) << gradNormHistory[storageCtr[0]] << "  ";
//   // 	for (int i = 0; i < primary.nSearchDim; i++)
//   // 	  cout << setw(23) << XHistory[storageCtr[0]][i] << "  ";
//   // 	cout << endl;
//   //     }
      
//   //     /* If termination due to reaching local minimum, summary should
//   //      * be the last computed point (this termination occurs before
//   //      * the next iteration is computed). */
//   //     if (isLocalMin)
//   //     {
//   // 	storageCtr[0]--;
//   // 	terminationCondition = 8;
//   //     }
//   //     else
//   // 	/* Check terminating conditions. */
//   // 	terminationCondition = checkTermination();
      
//   //   }while(terminationCondition == 0);

//   // }

  //!!! Hack final position. */
  XHistory[storageCtr] = XInitial;
  
  /* Perform high quality Monte Carlo evaluation at the final
   * position. Note we still take the negative of the DPCosts to be
   * consistent with what we did in the optimization algorithm due to
   * minimization implementation of the algorithm. */
  
  gsl_rng_set(generatorNoiseForHighQualityMC, rand() + algParams.rank);
  finalObjHighQualityMC = 0.0;

  for (int i = 0; i < algParams.nFinalObjHighQualityMC; i++)
  {
    //!!! here is just a temporary code for evaluating inside the
    //!!! expectation. If we use actual stochastic optimization we
    //!!! will need to use DPCostsInsideExpectation class.
    vector<double> disturbance(algParams.nDisturbanceDim, 0.0);
    generateDisturbance(algParams, state, XInitial, 
			generatorNoiseForHighQualityMC, disturbance);
    
    /* Evaluate system equation to compute state at stage k+1. */
    vector<double> newState(algParams.nStatesDim, 0.0);
    (*(allInputs.systemEqnPtr))(algParams, state, XInitial, disturbance, newState);
    
    /* Evaluate stage and future rewards. */
    if (((*(allInputs.futureFcnPtr)) != NULL) 
    	&& (allInputs.futureLinearArchClassPtr == NULL))
      /* J_{k+1} is terminal reward. */
      finalObjHighQualityMC += -((*(allInputs.futureFcnPtr))(algParams, newState)
      				 - (*(allInputs.stageFcnPtr))(algParams, state, 
      							      XInitial, disturbance));
    else if (((*(allInputs.futureFcnPtr)) == NULL) 
    	     && (allInputs.futureLinearArchClassPtr != NULL))
      /* J_{k+1} is intermediate reward. */
      finalObjHighQualityMC += -((allInputs.futureLinearArchClassPtr->evalArchitecture)
      				 (algParams, newState)
      				 - (*(allInputs.stageFcnPtr))(algParams, state, 
      							      XInitial, disturbance));

    else
    {
      cout << "Error: either both types of future function pointer defined or both NULL. Check function pointer assignments. " << endl;
      exit(1);
    }
  }  
  finalObjHighQualityMC /= double(algParams.nFinalObjHighQualityMC);
  
  /* Display summary. */
  if (algParams.displayOptSummary)
    displaySummary();
  
}

int StochasticSearch::checkTermination()
{

  /* In general, we need NOT to exclude counting changes that are
   * exactly 0 to the stopping rules, because even thought these are
   * usually from rejected steps in an algorithm, which are not an
   * indication of convergence, but such a "stuck" algorithm should be
   * stopped as no more progress is being made. */
  terminationCondition = 0;
  
  /* Norm of current position. This condition is not applied to NMNS,
   * because the change of the minimum point (which is stored in
   * history) is not representative of convergence. It could remain
   * unchanged between iterations even though the simplex has evolved
   * substantially. */
  curXNorm = vectorNorm(XHistory[storageCtr], algParams.normChoice);
  if (fabs(curXNorm - previousXNorm) < algParams.relXNormTerminateTol)
    nConsecRelXNorm++;
  else
    nConsecRelXNorm = 0;
  previousXNorm = curXNorm;
  
  /* Determine the conditions. */
  switch (algParams.optMethod)
  {
  case 6:
    /* SARM. */
    if (nIter >= algParams.maxOptIters)
      terminationCondition = 1;
    else if (nConsecRelXNorm >= algParams.nConsecRelXNormTerminateTol)
      terminationCondition = 4;
    break;

  default:
    cout << "Error: Optimization method " << algParams.optMethod
	 << " not available." << endl;
    exit(1);
  }

  return terminationCondition;

}

void StochasticSearch::displaySummary()
{
  
  /* Initializations. */
  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(16);

  /* Obtain termination label. */
  makeTerminationLabel();
  
  cout << endl
       << terminationLabel << endl
       << "% Total optimization iterations = " << iterHistory[storageCtr] << endl;
  cout << "% Final gradient norm           = " << setw(23) << gradNormHistory[storageCtr] << endl;
  
  cout << "% Final coordinates: ";
  for (unsigned int i = 0; i < XHistory[0].size(); i++)
    cout << "% " << setw(23) << XHistory[storageCtr][i];
  cout << endl;

}

double StochasticSearch::exportFinalObjHighQualityMC()
{
  return finalObjHighQualityMC;
}

void StochasticSearch::exportFinalPosition(vector<double> &finalXExternal)
{
  finalXExternal = XHistory[storageCtr];
}

void StochasticSearch::makeTerminationLabel()
{
  
  /* Initializations. */
  terminationLabel = "% Termination: ";
  
  /* Determine the termination condition. */
  switch (terminationCondition)
  {
  case 1:
    terminationLabel += "Total number of iterations reached.";
    break;
  case 2:
    terminationLabel += "Total number of algorithm-relevant function evaluations reached.";
    break;
  case 3:
    terminationLabel += "Total number of total function evaluations reached.";
    break;
  case 4:
    terminationLabel += "Norm of position converged.";
    break;
  case 5:
    terminationLabel += "Objective value converged.";
    break;
  case 6:
    terminationLabel += "Gradient norm converged.";
    break;
  case 7:
    terminationLabel += "Simplex diameter converged.";
    break;
  case 8:
    terminationLabel += "Local minimum reached.";
    break;
  default:
    cout << "Error: Termination condition " << terminationCondition
	 << " not applicable. Please check algorithm." << endl;
    exit(1);
  }
  
}



/*********************************************************************
 *********************************************************************/
