#include "fileRead.h"
#include "inputParams.h"

/*********************************************************************
 * InputParams
 *********************************************************************/

InputParams::InputParams()
{
  setDefaultValues();
}

InputParams::InputParams(string const fName)
{
  setDefaultValues();
  readInputParamsFile(fName);
}

InputParams::InputParams(InputParams const &other)
{
  
  /* General. */
  doubleCompTol = other.doubleCompTol;
  
  /* MPI. */
  nTasks = other.nTasks;
  rank = other.rank;
  status = other.status;
  
  /* Dynamic programming. */
  nStages = other.nStages;
  nStatesDim = other.nStatesDim;
  nControlsDim = other.nControlsDim;
  nDisturbanceDim = other.nDisturbanceDim;
  nInferenceParamsDim = other.nInferenceParamsDim;
  nSearchDim = other.nSearchDim;
  initialState = other.initialState;
  controlsLowerBounds = other.controlsLowerBounds;
  controlsUpperBounds = other.controlsUpperBounds;
  disturbanceStructure = other.disturbanceStructure;
  nStochOptPerState = other.nStochOptPerState;
  
  /* ADP value function approximation. */
  featuresChoice = other.featuresChoice;
  pOrder = other.pOrder;
  coefsConstructionMethod = other.coefsConstructionMethod;
  nRegressionSamples = other.nRegressionSamples;
  nFeatures = other.nFeatures;
  
  /* Stochastic Optimization. */
   optMethod = other.optMethod;
   maxOptIters = other.maxOptIters;
   relXNormTerminateNormChoice = other.relXNormTerminateNormChoice;
   relXNormTerminateTol = other.relXNormTerminateTol;
   nConsecRelXNormTerminateTol = other.nConsecRelXNormTerminateTol;
   normChoice = other.normChoice;
   nObjMC = other.nObjMC;
   userXInitial = other.userXInitial;
   randomizeXInitial = other.randomizeXInitial;
   nFinalObjHighQualityMC = other.nFinalObjHighQualityMC;
   displayOptProgress = other.displayOptProgress;
   displayOptSummary = other.displayOptSummary;
  
  // /* SARM. */
  //  checkInitialGradient = 0;
  //  SARMGainSeq = 1;
  //  gainMultiplier = 1.0;
  //  detectSPSAParams = 0;
  //  SPSAInputa = 0.16;
  //  SPSAInputc = 0.16;
  //  SPSAInputA = 100.0;
  //  SPSAInputalpha = 0.602;
  //  SPSAInputgamma = 0.101;
  
}

InputParams::~InputParams()
{
}

InputParams& InputParams::operator=(InputParams const &rhs)
{

  /* Protect against invalid self-assignment. */
  if (this != &rhs)
  {

    /* General. */
    doubleCompTol = rhs.doubleCompTol;
    
    /* MPI. */
    nTasks = rhs.nTasks;
    rank = rhs.rank;
    status = rhs.status;
    
    /* Dynamic programming. */
    nStages = rhs.nStages;
    nStatesDim = rhs.nStatesDim;
    nControlsDim = rhs.nControlsDim;
    nDisturbanceDim = rhs.nDisturbanceDim;
    nInferenceParamsDim = rhs.nInferenceParamsDim;
    nSearchDim = rhs.nSearchDim;
    initialState = rhs.initialState;
    controlsLowerBounds = rhs.controlsLowerBounds;
    controlsUpperBounds = rhs.controlsUpperBounds;
    disturbanceStructure = rhs.disturbanceStructure;
    nStochOptPerState = rhs.nStochOptPerState;
  
    /* ADP value function approximation. */
    featuresChoice = rhs.featuresChoice;
    pOrder = rhs.pOrder;
    coefsConstructionMethod = rhs.coefsConstructionMethod;
    nRegressionSamples = rhs.nRegressionSamples;
    nFeatures = rhs.nFeatures;
    
    /* Stochastic Optimization. */
    optMethod = rhs.optMethod;
    maxOptIters = rhs.maxOptIters;
    relXNormTerminateNormChoice = rhs.relXNormTerminateNormChoice;
    relXNormTerminateTol = rhs.relXNormTerminateTol;
    nConsecRelXNormTerminateTol = rhs.nConsecRelXNormTerminateTol;
    normChoice = rhs.normChoice;
    nObjMC = rhs.nObjMC;
    userXInitial = rhs.userXInitial;
    randomizeXInitial = rhs.randomizeXInitial;
    nFinalObjHighQualityMC = rhs.nFinalObjHighQualityMC;
    displayOptProgress = rhs.displayOptProgress;
    displayOptSummary = rhs.displayOptSummary;
    
    // /* SARM. */
    //  checkInitialGradient = 0;
    //  SARMGainSeq = 1;
    //  gainMultiplier = 1.0;
    //  detectSPSAParams = 0;
    //  SPSAInputa = 0.16;
    //  SPSAInputc = 0.16;
    //  SPSAInputA = 100.0;
    //  SPSAInputalpha = 0.602;
    //  SPSAInputgamma = 0.101;
    
  }
  
  /* By convention, always return *this. */
  return *this;
  
}

void InputParams::setDefaultValues()
{

  /* Set default values, no memory allocation has been made at this
   * point. */

  /* General. */
  cout.precision(16);
  cout.setf(ios::scientific);
//  srand(time(NULL));
  doubleCompTol = 1.0e-14;
  
  /* Dynamic programming. */
  nStages = -1;
  nStatesDim = -1;
  nControlsDim = -1;
  nDisturbanceDim = -1;
  nInferenceParamsDim = -1;
  nSearchDim = -1;
  disturbanceStructure = 0;
  nStochOptPerState = 1;
  
  /* ADP value function approximation. */
  featuresChoice = 1;
  pOrder = 0;
  coefsConstructionMethod = 1;
  nRegressionSamples = 1;
  nFeatures = 1;
  
  /* Stochastic Optimization. */
   optMethod = 1;
   maxOptIters = 0;
   relXNormTerminateNormChoice = 2;
   relXNormTerminateTol = 0.0;
   nConsecRelXNormTerminateTol = 1;
   normChoice = 2;
   nObjMC = 1;
   randomizeXInitial = 0;
   nFinalObjHighQualityMC = 1;
   displayOptProgress = 1;
   displayOptSummary = 1;
  
  // /* SARM. */
  //  checkInitialGradient = 0;
  //  SARMGainSeq = 1;
  //  gainMultiplier = 1.0;
  //  detectSPSAParams = 0;
  //  SPSAInputa = 0.16;
  //  SPSAInputc = 0.16;
  //  SPSAInputA = 100.0;
  //  SPSAInputalpha = 0.602;
  //  SPSAInputgamma = 0.101;

}

void InputParams::readInputParamsFile(string const fName)
{
  
  /* Initializations. */
  ifstream f;
  char data[256], name[256], val[256];
  int pos(0);
  
  /* Open file. */
  f.open(fName.c_str());
  
  /* Read valid lines. */
  while(!f.eof())
  {

    f.getline(data, 256);
    
    /* Skip empty lines and comments (starting a line with '%'). */
    while ((strlen(data) == 0 || data[0] == ' ' || data[0] == '%') 
	   && !f.eof())
      f.getline(data, 256);
    if (f.eof())
      break;

    /* Initialization for each line read. */
    pos = 0;

    /* Extract variable name. */
    while (data[pos] != ' ')
    {
      name[pos] = data[pos];
      pos++;
    }
    name[pos] = '\0';

    /* Read and assign variable values. */
    
    /* General. */
    if (strcmp(name, "doubleCompTol") == 0)
    {
      readScalarAsString(val, data, pos);
      doubleCompTol = atof(val);
    }
    
    /* Dynammic programming. */
    else if (strcmp(name, "nStages") == 0)
    {
      readScalarAsString(val, data, pos);
      nStages = atoi(val);
    }
    else if (strcmp(name, "nStatesDim") == 0)
    {
      readScalarAsString(val, data, pos);
      nStatesDim = atoi(val);

      if (nStatesDim != 2)
      {
	cout << "Error: Currently nStatesDim is only configured for 2. "
	     << endl;
	exit(1);
      }
    }
    else if (strcmp(name, "nControlsDim") == 0)
    {
      readScalarAsString(val, data, pos);
      nControlsDim = atoi(val);

      nSearchDim = nControlsDim;
    }
    else if (strcmp(name, "nDisturbanceDim") == 0)
    {
      readScalarAsString(val, data, pos);
      nDisturbanceDim = atoi(val);
    }
    else if (strcmp(name, "nInferenceParamsDim") == 0)
    {
      readScalarAsString(val, data, pos);
      nInferenceParamsDim = atoi(val);
    }
    else if (strcmp(name, "initialState") == 0)
    {
      /* This must come after reading nControlsDim. */
      if (nStatesDim > 0)
      {
	initialState.assign(nStatesDim, 0.0);
	readVector<double>(initialState, DOUBLE_TYPE, data, pos);
      }
      else
      {
	cout << "Error: nStatesDim must be specified before initialState." << endl;
	exit(1);
      }
    }
    else if (strcmp(name, "controlsLowerBounds") == 0)
    {
      /* This must come after reading nControlsDim. */
      if (nControlsDim > 0)
      {
	controlsLowerBounds.assign(nControlsDim, 0.0);
	readVector<double>(controlsLowerBounds, DOUBLE_TYPE, data, pos);
      }
      else
      {
	cout << "Error: nControlsDim must be specified before controlsLeftBound." << endl;
	exit(1);
      }
    }
    else if (strcmp(name, "controlsUpperBounds") == 0)
    {
      /* This must come after reading nControlsDim. */
      if (nControlsDim > 0)
      {
	controlsUpperBounds.assign(nControlsDim, 0.0);
	readVector<double>(controlsUpperBounds, DOUBLE_TYPE, data, pos);
      }
      else
      {
	cout << "Error: nControlsDim must be specified before controlsRightBound." << endl;
	exit(1);
      }
    }
    else if (strcmp(name, "disturbanceStructure") == 0)
    {
      readScalarAsString(val, data, pos);
      disturbanceStructure = atoi(val);
    }
    else if (strcmp(name, "nStochOptPerState") == 0)
    {
      readScalarAsString(val, data, pos);
      nStochOptPerState = atoi(val);
    }
    
    /* ADP value function approximation. */
    else if (strcmp(name, "featuresChoice") == 0)
    {
      readScalarAsString(val, data, pos);
      featuresChoice = atoi(val);
    }
    else if (strcmp(name, "pOrder") == 0)
    {
      readScalarAsString(val, data, pos);
      pOrder = atoi(val);
    }
    else if (strcmp(name, "coefsConstructionMethod") == 0)
    {
      readScalarAsString(val, data, pos);
      coefsConstructionMethod = atoi(val);
    }
    else if (strcmp(name, "nRegressionSamples") == 0)
    {
      readScalarAsString(val, data, pos);
      nRegressionSamples = atoi(val);
    }
    
    /* Stochastic Optimization. */
    else if (strcmp(name, "optMethod") == 0)
    {
      readScalarAsString(val, data, pos);
      optMethod = atoi(val);
    }
    else if (strcmp(name, "maxOptIters") == 0)
    {
      readScalarAsString(val, data, pos);
      maxOptIters = atoi(val);
    }
    else if (strcmp(name, "relXNormTerminateNormChoice") == 0)
    {
      readScalarAsString(val, data, pos);
      relXNormTerminateNormChoice = atoi(val);
    }
    else if (strcmp(name, "relXNormTerminateTol") == 0)
    {
      readScalarAsString(val, data, pos);
      relXNormTerminateTol = atof(val);
    }
    else if (strcmp(name, "nConsecRelXNormTerminateTol") == 0)
    {
      readScalarAsString(val, data, pos);
      nConsecRelXNormTerminateTol = atoi(val);
    }
    else if (strcmp(name, "normChoice") == 0){
      readScalarAsString(val, data, pos);
      normChoice = atoi(val);
    }
    else if (strcmp(name, "nObjMC") == 0)
    {
      readScalarAsString(val, data, pos);
      nObjMC = atoi(val);
    }
    else if (strcmp(name, "userXInitial") == 0)
    {
      /* This must come after assigning nSearchDim. */
      if (nSearchDim > 0)
      {
	userXInitial.assign(nSearchDim, 0.0);
	readVector<double>(userXInitial, DOUBLE_TYPE, data, pos);
      }
      else
      {
    	cout << "Error: nSearchDim must be specified before userXInitial." << endl;
    	exit(1);
      }
    }
    else if (strcmp(name, "randomizeXInitial") == 0)
    {
      readScalarAsString(val, data, pos);
      randomizeXInitial = atoi(val);
    }
    else if (strcmp(name, "nFinalObjHighQualityMC") == 0)
    {
      readScalarAsString(val, data, pos);
      nFinalObjHighQualityMC = atoi(val);
    }
    else if (strcmp(name, "displayOptProgress") == 0)
    {
      readScalarAsString(val, data, pos);
      displayOptProgress = atoi(val);
    }
    else if (strcmp(name, "displayOptSummary") == 0)
    {
      readScalarAsString(val, data, pos);
      displayOptSummary = atoi(val);
    }
    
    // /* SARM. */
    // else if (strcmp(name, "checkInitialGradient") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    checkInitialGradient = atoi(val);
    // }
    // else if (strcmp(name, "SARMGainSeq") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    SARMGainSeq = atoi(val);
    // }
    // else if (strcmp(name, "gainMultiplier") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    gainMultiplier = atof(val);
    // }
    // else if (strcmp(name, "detectSPSAParams") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    detectSPSAParams = atoi(val);
    // }
    // else if (strcmp(name, "SPSAInputa") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    SPSAInputa = atof(val);
    // }
    // else if (strcmp(name, "SPSAInputc") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    SPSAInputc = atof(val);
    // }
    // else if (strcmp(name, "SPSAInputA") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    SPSAInputA = atof(val);
    // }
    // else if (strcmp(name, "SPSAInputalpha") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    SPSAInputalpha = atof(val);
    // }
    // else if (strcmp(name, "SPSAInputgamma") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    SPSAInputgamma = atof(val);
    // }
    
    else
    {
      cout << "Error: Unrecognized input variable " << name << '.' << endl;
      exit(1);
    }
  }
  
  /* Close files. */
  f.close();
  
  /* Compute total number of features. */
  switch (featuresChoice)
  {
  case 1:
    /* Total order polynomial. */
    nFeatures = nCr(pOrder + nStatesDim, pOrder);
    break;
    
  case 2:
    /* Terms from Gaussian KL divergence: 1.0, mean, mean-squared,
     * variance, log-variance. In that order. */
    nFeatures = 5;
    break;
    
  default:
    cout << "Error: Features choice " << featuresChoice
	 << " not available." << endl;
    exit(1);
  }
  
  /* Post processing. */
  // nFinalObjHighQualityMC = max<int>( nFinalObjHighQualityMC,
  // 				     nObjMC);

}

void InputParams::addMPIInfo(int const nTasksRef, int const rankRef)
{
  
  /* Store MPI variables. */
   nTasks = nTasksRef;
   rank = rankRef;
  
}



/*********************************************************************
 *********************************************************************/
