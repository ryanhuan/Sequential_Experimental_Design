#include "fileRead.h"
#include "inputParams.h"

/*********************************************************************
 * InputParams
 *********************************************************************/

void initializeControls(Controls &primary)
{
  
  /* General. */
  cout.precision(16);
  cout.setf(ios::scientific);
//  srand(time(NULL));
  primary.dCompTol = 1.0e-14;
  
  /* Dynamic programming. */
  primary.nStages = -1;
  primary.nStatesDim = -1;
  primary.initialState = NULL;
  primary.nControlsDim = -1;
  primary.controlsLeftBound = NULL;
  primary.controlsRightBound = NULL;
  primary.nNoiseDim = -1;
  primary.noiseStdDevConstant = 1.0;
  primary.nParamsDim = -1;
  primary.nStochOptPerStateK = 1;
  primary.stageCostQuadControlsWeight = 0.0;
  
  /* ADP value function approximation. */
  primary.featuresChoice = 1;
  primary.coefsConstructionMethod = 1;
  primary.nRegressionSamples = 1;
  primary.pOrder = 0;
  primary.nFeatures = 1;

  /* Stochastic Optimization. */
  primary.optMethod = 1;
  primary.maxOptIters = 0;
  primary.relXNormTerminateNormChoice = 2;
  primary.relXNormTerminateTol = 0.0;
  primary.nConsecRelXNormTerminateTol = 1;
  primary.gradNormChoice = 2;
  primary.nObjMC = 1;
  primary.XInitial = NULL;
  primary.randomizeXInitial = 0;
  primary.nFinalObjHighQualityMC = 1;
  primary.displayOptProgress = 1;
  primary.displayOptSummary = 1;
  
  /* SARM. */
  primary.checkInitialGradient = 0;
  primary.SARMGainSeq = 1;
  primary.gainMultiplier = 1.0;
  primary.detectSPSAParams = 0;
  primary.SPSAInputa = 0.16;
  primary.SPSAInputc = 0.16;
  primary.SPSAInputA = 100.0;
  primary.SPSAInputalpha = 0.602;
  primary.SPSAInputgamma = 0.101;

  /* Read any Controls input from file. */
  readFileControls(primary);

}

void addMPIInfoControls(Controls &primary, int const nTasks, 
			int const rank)
{
  
  /* Store MPI variables. */
  primary.nTasks = nTasks;
  primary.rank = rank;
  
}

void finalizeControls(Controls &primary)
{

  /* Free memory. */
  delete [] primary.initialState;
  delete [] primary.controlsLeftBound;
  delete [] primary.controlsRightBound;
  delete [] primary.XInitial;
  
}
