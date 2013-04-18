#include "fileRead.h"

/* Functions with templates are defined in header file. */

void readFileControls(Controls &primary)
{
  
  /* Initializations. */
  string fName("controls.inp");
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
    if (strcmp(name, "dCompTol") == 0)
    {
      readScalar(val, data, pos);
      primary.dCompTol = atof(val);
    }

    /* Dynammic programming. */
    else if (strcmp(name, "nStages") == 0)
    {
      readScalar(val, data, pos);
      primary.nStages = atoi(val);
    }
    else if (strcmp(name, "nStatesDim") == 0)
    {
      readScalar(val, data, pos);
      primary.nStatesDim = atoi(val);
      primary.initialState = new double[primary.nStatesDim];

      if (primary.nStatesDim != 2)
      {
	cout << "Error: Currently nStatesDim is only configured for 2. "
	     << endl;
	exit(1);
      }
    }
    else if (strcmp(name, "initialState") == 0)
    {
      /* This must come after reading nStatesDim. */
      if (primary.nStatesDim > 0)
	readVector<double>(primary.initialState, DOUBLE_TYPE,
			   primary.nStatesDim, data, pos);
      else
      {
	cout << "Error: nStatesDim must be specified before initialState." << endl;
	exit(1);
      }
    }
    else if (strcmp(name, "nControlsDim") == 0)
    {
      readScalar(val, data, pos);
      primary.nControlsDim = atoi(val);
      primary.controlsLeftBound = new double[primary.nControlsDim];
      primary.controlsRightBound = new double[primary.nControlsDim];
      primary.XInitial = new double[primary.nControlsDim];
    }
    else if (strcmp(name, "controlsLeftBound") == 0)
    {
      /* This must come after reading nControlsDim. */
      if (primary.nControlsDim > 0)
	readVector<double>(primary.controlsLeftBound, DOUBLE_TYPE,
			   primary.nControlsDim, data, pos);
      else
      {
	cout << "Error: nControlsDim must be specified before controlsLeftBound." << endl;
	exit(1);
      }
    }
    else if (strcmp(name, "controlsRightBound") == 0)
    {
      /* This must come after reading nControlsDim. */
      if (primary.nControlsDim > 0)
	readVector<double>(primary.controlsRightBound, DOUBLE_TYPE,
			   primary.nControlsDim, data, pos);
      else
      {
	cout << "Error: nControlsDim must be specified before controlsRightBound." << endl;
	exit(1);
      }
    }
    else if (strcmp(name, "nNoiseDim") == 0)
    {
      readScalar(val, data, pos);
      primary.nNoiseDim = atoi(val);
    }
    else if (strcmp(name, "noiseStdDevConstant") == 0)
    {
      readScalar(val, data, pos);
      primary.noiseStdDevConstant = atof(val);
    }
    else if (strcmp(name, "nParamsDim") == 0)
    {
      readScalar(val, data, pos);
      primary.nParamsDim = atoi(val);
    }
    else if (strcmp(name, "nStochOptPerStateK") == 0)
    {
      readScalar(val, data, pos);
      primary.nStochOptPerStateK = atoi(val);
    }
    else if (strcmp(name, "stageCostQuadControlsWeight") == 0)
    {
      readScalar(val, data, pos);
      primary.stageCostQuadControlsWeight = atof(val);
    }
    
    /* ADP value function approximation. */
    else if (strcmp(name, "featuresChoice") == 0)
    {
      readScalar(val, data, pos);
      primary.featuresChoice = atoi(val);
    }
    else if (strcmp(name, "coefsConstructionMethod") == 0)
    {
      readScalar(val, data, pos);
      primary.coefsConstructionMethod = atoi(val);
    }
    else if (strcmp(name, "nRegressionSamples") == 0)
    {
      readScalar(val, data, pos);
      primary.nRegressionSamples = atoi(val);
    }
    else if (strcmp(name, "pOrder") == 0)
    {
      readScalar(val, data, pos);
      primary.pOrder = atoi(val);
    }
    
    /* Stochastic Optimization. */
    else if (strcmp(name, "optMethod") == 0)
    {
      readScalar(val, data, pos);
      primary.optMethod = atoi(val);
    }
    else if (strcmp(name, "maxOptIters") == 0)
    {
      readScalar(val, data, pos);
      primary.maxOptIters = atoi(val);
    }
    else if (strcmp(name, "relXNormTerminateNormChoice") == 0)
    {
      readScalar(val, data, pos);
      primary.relXNormTerminateNormChoice = atoi(val);
    }
    else if (strcmp(name, "relXNormTerminateTol") == 0)
    {
      readScalar(val, data, pos);
      primary.relXNormTerminateTol = atof(val);
    }
    else if (strcmp(name, "nConsecRelXNormTerminateTol") == 0)
    {
      readScalar(val, data, pos);
      primary.nConsecRelXNormTerminateTol = atoi(val);
    }
    else if (strcmp(name, "gradNormChoice") == 0){
      readScalar(val, data, pos);
      primary.gradNormChoice = atoi(val);
    }
    else if (strcmp(name, "nObjMC") == 0)
    {
      readScalar(val, data, pos);
      primary.nObjMC = atoi(val);
    }
    else if (strcmp(name, "XInitial") == 0)
    {
      /* This must come after reading nDim. */
      if (primary.nControlsDim > 0)
    	readVector<double>(primary.XInitial, DOUBLE_TYPE,
    			   primary.nControlsDim, data, pos);
      else
      {
    	cout << "Error: nControlsDim must be specified before XInitial." << endl;
    	exit(1);
      }
    }
    else if (strcmp(name, "nFinalObjHighQualityMC") == 0)
    {
      readScalar(val, data, pos);
      primary.nFinalObjHighQualityMC = atoi(val);
    }
    else if (strcmp(name, "randomizeXInitial") == 0)
    {
      readScalar(val, data, pos);
      primary.randomizeXInitial = atoi(val);
    }
    else if (strcmp(name, "displayOptProgress") == 0)
    {
      readScalar(val, data, pos);
      primary.displayOptProgress = atoi(val);
    }
    else if (strcmp(name, "displayOptSummary") == 0)
    {
      readScalar(val, data, pos);
      primary.displayOptSummary = atoi(val);
    }
    
    /* SARM. */
    else if (strcmp(name, "checkInitialGradient") == 0)
    {
      readScalar(val, data, pos);
      primary.checkInitialGradient = atoi(val);
    }
    else if (strcmp(name, "SARMGainSeq") == 0)
    {
      readScalar(val, data, pos);
      primary.SARMGainSeq = atoi(val);
    }
    else if (strcmp(name, "gainMultiplier") == 0)
    {
      readScalar(val, data, pos);
      primary.gainMultiplier = atof(val);
    }
    else if (strcmp(name, "detectSPSAParams") == 0)
    {
      readScalar(val, data, pos);
      primary.detectSPSAParams = atoi(val);
    }
    else if (strcmp(name, "SPSAInputa") == 0)
    {
      readScalar(val, data, pos);
      primary.SPSAInputa = atof(val);
    }
    else if (strcmp(name, "SPSAInputc") == 0)
    {
      readScalar(val, data, pos);
      primary.SPSAInputc = atof(val);
    }
    else if (strcmp(name, "SPSAInputA") == 0)
    {
      readScalar(val, data, pos);
      primary.SPSAInputA = atof(val);
    }
    else if (strcmp(name, "SPSAInputalpha") == 0)
    {
      readScalar(val, data, pos);
      primary.SPSAInputalpha = atof(val);
    }
    else if (strcmp(name, "SPSAInputgamma") == 0)
    {
      readScalar(val, data, pos);
      primary.SPSAInputgamma = atof(val);
    }
    
    else
    {
      cout << "Error: Unrecognized input variable " << name << '.' << endl;
      exit(1);
    }
  }
  
  /* Close files. */
  f.close();
  
  /* Compute total number of features. */
  switch (primary.featuresChoice)
  {
  case 1:
    /* Total order polynomial. */
    primary.nFeatures = nCr(primary.pOrder + primary.nStatesDim, primary.pOrder);
    break;

  case 2:
    /* Terms from Gaussian KL divergence: 1.0, mean, mean-squared,
     * variance, log-variance. In that order. */
    primary.nFeatures = 5;
    break;
    
  default:
    cout << "Error: Features choice " << primary.featuresChoice
	 << " not available." << endl;
    exit(1);
  }

  /* Post processing. */
  primary.nFinalObjHighQualityMC = max<int>(primary.nFinalObjHighQualityMC,
					    primary.nObjMC);

}

void readScalar(char* const val, char const * const data, int &pos)
{

  /* Initializations. */
  int ctrVal(0);
  
  /* Find starting point for numeric value. */
  while (data[pos] == ' ' || data[pos] == '=')
    pos++;
  
  /* Record numeric value. */
  while (data[pos] != ' ' && data[pos] != '\0')
  {
    val[ctrVal] = data[pos];
    ctrVal++;
    pos++;
  }
  val[ctrVal] = '\0';
  
}
