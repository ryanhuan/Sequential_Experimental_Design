#include "stochasticSearch.h"

/*********************************************************************
 * StochasticSearch
 *********************************************************************/

StochasticSearch::StochasticSearch(Controls const &refControls)
{
  
  /* Load controls. */
  primary = refControls;
  
  /* Allocate memory. */
  nIter = new int[1];
  storageCtr = new int[1];

  iterHistory = new int[primary.maxOptIters + 1];
  XHistory = new double*[primary.maxOptIters + 1];
  XHistory[0] = new double[(primary.maxOptIters + 1) * primary.nControlsDim];
  for (int i = 1; i < primary.maxOptIters + 1; i++)
    XHistory[i] = XHistory[i - 1] + primary.nControlsDim;
  valHistory = new double[primary.maxOptIters + 1];

  DPCostsForHighQualityMC = new DPCostsInsideExpectation(refControls);

  /* Initialize specific search algorithm class. */
  switch (primary.optMethod)
  {
  case 6:
    /* SARM. */
    gradNormHistory = new double[primary.maxOptIters + 1];
    robbinsMonro = new SARM(primary, nIter, storageCtr, iterHistory, XHistory,
    			    valHistory, gradNormHistory);
    break;
    
  default:
    cout << "Error: Optimization method " << primary.optMethod
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

}

StochasticSearch::~StochasticSearch()
{
  
  /* Free memory. */
  delete [] nIter;
  delete [] storageCtr;

  delete [] iterHistory;
  delete [] XHistory[0];
  delete [] XHistory;
  delete [] valHistory;

  delete DPCostsForHighQualityMC;

  switch (primary.optMethod)
  {
  case 6:
    delete robbinsMonro;
    delete [] gradNormHistory;
    break;
  default:
    cout << "Error: Optimization method " << primary.optMethod
	 << " not available." << endl;
    exit(1);
  }

  gsl_rng_free(generatorNoiseForHighQualityMC);

}

void StochasticSearch::makeLabels()
{
  
  /* Construct progress and text file label. */
  if (primary.optMethod == 6)
    /* SARM. */
    dataLabels = string("% Iter,    Gradient Norm,       Coordinates [")
      + num2string<int>(primary.nControlsDim) + string(" columns]");
  else
  {
    cout << "Error: Optimization method " << primary.optMethod
	 << " not available." << endl;
    exit(1);
  }
  
}

void StochasticSearch::initialPosition(double (*innerFcn)(Controls const &, 
							  double const * const), 
				       double const * const state)
{
  
  /* Some storage initializations. */
  nIter[0] = 0;
  storageCtr[0] = 0;
  nConsecRelXNorm = 0;
  isLocalMin = 0;  
  
  /* More initializations. */
  iterHistory[0] = 0;
  
  /* Call additional initialization functions that are specific to the
   * selected algorithm. For example, gradient estimate if a
   * gradient-based algorithm is used. */
  switch (primary.optMethod)
  {
  case 6:
    /* SARM. */
    robbinsMonro -> initialPosition(innerFcn, state);
    break;
  default:
    cout << "Error: Optimization method " << primary.optMethod
	 << " not available." << endl;
    exit(1);
  }
  
  /* Set output format and display label. */
  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(16);
  
  /* Display initial conditions. */
  if (primary.displayOptProgress)
  {
    cout << dataLabels << endl;
    cout << setw(4) << iterHistory[storageCtr[0]] << "  ";
    cout << setw(23) << gradNormHistory[storageCtr[0]] << "  ";
    for (int i = 0; i < primary.nControlsDim; i++)
      cout << setw(23) << XHistory[storageCtr[0]][i] << "  ";
    cout << endl;
  }
  
}

void StochasticSearch::search(double (*innerFcn)(Controls const &, 
						 double const * const), 
			      double const * const state)
{

  // /* Proceed only if at least 1 iteration is requested. */
  // if (primary.maxOptIters > 0)
  // {

  //   /* Initializations. */
  //   previousXNorm = normOf1DArray(XHistory[0], primary.nControlsDim, 
  // 				  primary.relXNormTerminateNormChoice);
    
  //   do
  //   {
      
  //     /* If the user wants to store and display progress every n
  //      * iterations, then a for-loop should be placed inside the
  //      * switch statement, wrapping the iterations. */
  //     switch (primary.optMethod)
  //     {
  //     case 6:
  // 	/* SARM. */
  // 	robbinsMonro -> iterateOnce(coefsJkp1Ref, refTableJkp1Ref, stateKRef);
  // 	break;
  //     default:
  // 	cout << "Error: Optimization method " << primary.optMethod
  // 	     << " not available." << endl;
  // 	exit(1);
  //     }

  //     /* Display progress. If local minimum is reached, the next point
  //      * is not computed and thus should not be displayed (it would
  //      * display gibberish). */
  //     if (primary.displayOptProgress && !isLocalMin)
  //     {
  // 	cout << setw(4) << iterHistory[storageCtr[0]] << "  ";
  // 	cout << setw(23) << gradNormHistory[storageCtr[0]] << "  ";
  // 	for (int i = 0; i < primary.nControlsDim; i++)
  // 	  cout << setw(23) << XHistory[storageCtr[0]][i] << "  ";
  // 	cout << endl;
  //     }
      
  //     /* If termination due to reaching local minimum, summary should
  //      * be the last computed point (this termination occurs before
  //      * the next iteration is computed). */
  //     if (isLocalMin)
  //     {
  // 	storageCtr[0]--;
  // 	terminationCondition = 8;
  //     }
  //     else
  // 	/* Check terminating conditions. */
  // 	terminationCondition = checkTermination();
      
  //   }while(terminationCondition == 0);

  // }

  /* Perform high quality Monte Carlo evaluation at the final
   * position. Note we still take the negative of the DPCosts to be
   * consistent with what we did in the optimization algorithm due to
   * minimization implementation of the algorithm. */

  gsl_rng_set(generatorNoiseForHighQualityMC, rand() + primary.rank);
  finalObjHighQualityMC = 0.0;
  for (int i = 0; i < primary.nFinalObjHighQualityMC; i++)
    finalObjHighQualityMC += -DPCostsForHighQualityMC -> 
      sampleOnce(innerFcn, state, primary.XInitial, generatorNoiseForHighQualityMC);
  finalObjHighQualityMC /= double(primary.nFinalObjHighQualityMC);
  
  /* Display summary. */
  if (primary.displayOptSummary)
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
  curXNorm = normOf1DArray(XHistory[storageCtr[0]], primary.nControlsDim,
			   primary.relXNormTerminateNormChoice);
  if (fabs(curXNorm - previousXNorm) < primary.relXNormTerminateTol)
    nConsecRelXNorm++;
  else
    nConsecRelXNorm = 0;
  previousXNorm = curXNorm;
  
  /* Determine the conditions. */
  switch (primary.optMethod)
  {
  case 6:
    /* SARM. */
    if (nIter[0] >= primary.maxOptIters)
      terminationCondition = 1;
    else if (nConsecRelXNorm >= primary.nConsecRelXNormTerminateTol)
      terminationCondition = 4;
    break;

  default:
    cout << "Error: Optimization method " << primary.optMethod
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
       << "% Total optimization iterations = " << iterHistory[storageCtr[0]] << endl;
  cout << "% Final gradient norm           = " << setw(23) << gradNormHistory[storageCtr[0]] << endl;
  
  cout << "% Final coordinates: ";
  for (int i = 0; i < primary.nControlsDim; i++)
    cout << "% " << setw(23) << XHistory[storageCtr[0]][i];
  cout << endl;

}

double StochasticSearch::exportFinalObjHighQualityMC()
{
  return finalObjHighQualityMC;
}

void StochasticSearch::exportFinalPosition(double * const finalXExternal)
{
  for (int i = 0; i < primary.nControlsDim; i++)
    finalXExternal[i] = XHistory[storageCtr[0]][i];
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
