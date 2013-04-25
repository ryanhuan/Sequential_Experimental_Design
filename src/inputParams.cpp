#include "fileRead.h"
#include "inputParams.h"

/*********************************************************************
 * InputParams
 *********************************************************************/

InputParams::InputParams(string const fName)
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
  nStochOptPerState = 1;
  //!!! check
  // systemEqnPtr = systemEquation;
  // stageCostPtr = stageCost;
  
  /* ADP value function approximation. */
  featuresChoice = 1;
  pOrder = 0;
  coefsConstructionMethod = 1;
  nRegressionSamples = 1;
  nFeatures = 1;
  
  // /* Stochastic Optimization. */
  //  optMethod = 1;
  //  maxOptIters = 0;
  //  relXNormTerminateNormChoice = 2;
  //  relXNormTerminateTol = 0.0;
  //  nConsecRelXNormTerminateTol = 1;
  //  gradNormChoice = 2;
  //  nObjMC = 1;
  //  XInitial = NULL;
  //  randomizeXInitial = 0;
  //  nFinalObjHighQualityMC = 1;
  //  displayOptProgress = 1;
  //  displayOptSummary = 1;
  
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

  /* Read input file, all memory allocation is made externally as file
   * is read. */
  readInputParamsFile(fName);

}

InputParams::~InputParams()
{

  /* Free memory. */
  // delete []  XInitial;
  
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

      controlsLowerBounds.assign(nControlsDim, 0.0);
      controlsUpperBounds.assign(nControlsDim, 0.0);
//      XInitial = new double[ nControlsDim];
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
    else if (strcmp(name, "controlsLowerBounds") == 0)
    {
      /* This must come after reading nControlsDim. */
      if (nControlsDim > 0)
	readVector<double>(controlsLowerBounds, DOUBLE_TYPE, data, pos);
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
	readVector<double>(controlsUpperBounds, DOUBLE_TYPE, data, pos);
      else
      {
	cout << "Error: nControlsDim must be specified before controlsRightBound." << endl;
	exit(1);
      }
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
    
    // /* Stochastic Optimization. */
    // else if (strcmp(name, "optMethod") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    optMethod = atoi(val);
    // }
    // else if (strcmp(name, "maxOptIters") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    maxOptIters = atoi(val);
    // }
    // else if (strcmp(name, "relXNormTerminateNormChoice") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    relXNormTerminateNormChoice = atoi(val);
    // }
    // else if (strcmp(name, "relXNormTerminateTol") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    relXNormTerminateTol = atof(val);
    // }
    // else if (strcmp(name, "nConsecRelXNormTerminateTol") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    nConsecRelXNormTerminateTol = atoi(val);
    // }
    // else if (strcmp(name, "gradNormChoice") == 0){
    //   readScalarAsString(val, data, pos);
    //    gradNormChoice = atoi(val);
    // }
    // else if (strcmp(name, "nObjMC") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    nObjMC = atoi(val);
    // }
    // else if (strcmp(name, "XInitial") == 0)
    // {
    //   /* This must come after reading nDim. */
    //   if ( nControlsDim > 0)
    // 	readVector<double>( XInitial, DOUBLE_TYPE,
    // 			    nControlsDim, data, pos);
    //   else
    //   {
    // 	cout << "Error: nControlsDim must be specified before XInitial." << endl;
    // 	exit(1);
    //   }
    // }
    // else if (strcmp(name, "nFinalObjHighQualityMC") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    nFinalObjHighQualityMC = atoi(val);
    // }
    // else if (strcmp(name, "randomizeXInitial") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    randomizeXInitial = atoi(val);
    // }
    // else if (strcmp(name, "displayOptProgress") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    displayOptProgress = atoi(val);
    // }
    // else if (strcmp(name, "displayOptSummary") == 0)
    // {
    //   readScalarAsString(val, data, pos);
    //    displayOptSummary = atoi(val);
    // }
    
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
    //!!!
    // nFeatures = nCr( pOrder +  nStatesDim,  pOrder);
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
