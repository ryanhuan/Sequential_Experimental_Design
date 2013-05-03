#include "linearArch.h"

/*********************************************************************
 * LinearArch
 *********************************************************************/

LinearArch::LinearArch()
{
  initialized = 0;
}

LinearArch::LinearArch(InputParams const &refAlgParams)
{
  initialize(refAlgParams);
}

LinearArch::~LinearArch()
{

  if (initialized)
  {
    /* Free memory. */
    gsl_rng_free(generator);
  }
  
}

LinearArch::LinearArch(LinearArch const &other)
{

  /* Initializations. */
  algParams = other.algParams;

  featuresCoefs = other.featuresCoefs;
  featuresVals = other.featuresVals;
  stateSample = other.stateSample;

  /* Total order polynomial features. */
  if (algParams.featuresChoice == 1)
    refTable = other.refTable;
  
  /* Construct coefficients via construction. */
  if (algParams.coefsConstructionMethod == 1)
  {
    /* Note ATrans is stored in its transposed form from the algebraic
     * form of the matrix to accomodate LAPACK (which uses its
     * transpose). Note ATrans is overwritten by right singular
     * vectors upon output of linear least squares computation. LHS
     * (in non-transposed form) will be used to store the original
     * matrix. */
    ATrans = other.ATrans;
    LHS = other.LHS;
    B = other.B;
    RHS = other.RHS;
    /* Size of soln array required to be at least
     * algParams.nRegressionSamples for LAPACK. */
    soln = other.soln;
    
    singularValues = other.singularValues;
    lwork = other.lwork;
    work = other.work;
  }

  /* Random number generator initialization. */
  rngType = gsl_rng_ranlxs0;
  generator = gsl_rng_alloc(rngType);
  gsl_rng_env_setup();  

  initialized = other.initialized;

}

LinearArch& LinearArch::operator=(LinearArch const &rhs)
{

  /* Protect against invalid self-assignment. */
  if (this != &rhs)
  {

    /* Initializations. */
    algParams = rhs.algParams;
    
    featuresCoefs = rhs.featuresCoefs;
    featuresVals = rhs.featuresVals;
    stateSample = rhs.stateSample;
    
    /* Total order polynomial features. */
    if (algParams.featuresChoice == 1)
      refTable = rhs.refTable;
    
    /* Construct coefficients via construction. */
    if (algParams.coefsConstructionMethod == 1)
    {
      /* Note ATrans is stored in its transposed form from the
       * algebraic form of the matrix to accomodate LAPACK (which uses
       * its transpose). Note ATrans is overwritten by right singular
       * vectors upon output of linear least squares computation. LHS
       * (in non-transposed form) will be used to store the original
       * matrix. */
      ATrans = rhs.ATrans;
      LHS = rhs.LHS;
      B = rhs.B;
      RHS = rhs.RHS;
      /* Size of soln array required to be at least
       * algParams.nRegressionSamples for LAPACK. */
      soln = rhs.soln;
      
      singularValues = rhs.singularValues;
      lwork = rhs.lwork;
      work = rhs.work;
    }
    
    /* Random number generator initialization. */
    rngType = gsl_rng_ranlxs0;
    generator = gsl_rng_alloc(rngType);
    gsl_rng_env_setup();  

    initialized = rhs.initialized;

  }
  
  /* By convention, always return *this. */
  return *this;
  
}

void LinearArch::initialize(InputParams const &refAlgParams)
{
  
  /* Initializations. */
  algParams = refAlgParams;

  featuresCoefs.assign(algParams.nFeatures, 0.0);
  featuresVals.assign(algParams.nFeatures, 0.0);
  stateSample.assign(algParams.nStatesDim, 0.0);

  /* Total order polynomial features. */
  if (algParams.featuresChoice == 1)
  {
    /* Make refTable. */
    refTable.assign(algParams.nFeatures, vector<int>(algParams.nStatesDim, 0));
    int ctrCoords(0);
    vector<int> tempCoords(algParams.nStatesDim, 0);
    permPolyOrders(algParams.nStatesDim, algParams.pOrder, ctrCoords, tempCoords);
  }
  
  /* Construct coefficients via construction. */
  if (algParams.coefsConstructionMethod == 1)
  {
    /* Memory allocation for regression method. Note ATrans is stored
     * in its transposed form from the algebraic form of the matrix to
     * accomodate LAPACK (which uses its transpose). Note ATrans is
     * overwritten by right singular vectors upon output of linear
     * least squares computation. LHS (in non-transposed form) will be
     * used to store the original matrix. Note matrices in the form of
     * vector of vectors will NOT have contiguous memory
     * allocation. */
    ATrans.assign(algParams.nFeatures, vector<double>(algParams.nRegressionSamples, 0.0));
    LHS.assign(algParams.nRegressionSamples, vector<double>(algParams.nFeatures, 0.0));
    B.assign(algParams.nRegressionSamples, 0.0);
    RHS.assign(algParams.nRegressionSamples, 0.0);
    /* Size of soln array required to be at least
     * algParams.nRegressionSamples for LAPACK. */
    soln.assign(algParams.nFeatures, 0.0);

    singularValues.assign(algParams.nFeatures, 0.0);
    lwork = 3 * min<int>(algParams.nFeatures, algParams.nRegressionSamples) 
      + max<int>(2 * min<int>(algParams.nFeatures, algParams.nRegressionSamples),
		 max<int>(algParams.nFeatures, algParams.nRegressionSamples));
    work.assign(lwork, 0.0);
  }
  
  /* Random number generator initialization. */
  rngType = gsl_rng_ranlxs0;
  generator = gsl_rng_alloc(rngType);
  gsl_rng_env_setup();

  initialized = 1;
  
}

void LinearArch::permPolyOrders(int zetaRemain, int const upperLimit, 
				int &ctrCoords, vector<int> &tempCoords)
{

  /* Base case. */
  if (zetaRemain == 1)
  {
    for (int z = 0; z < upperLimit + 1; z++)
    {
      tempCoords[algParams.nStatesDim - 1] = z;
      for (int i = 0; i < algParams.nStatesDim; i++)
	refTable[ctrCoords][i] = tempCoords[i];
      ctrCoords++;
    }
  }

  /* Recusrive case. */
  else
  {
    zetaRemain--;

    for (int z = 0; z < upperLimit + 1; z++)
    {
      tempCoords[algParams.nStatesDim - 1 - zetaRemain] = z;
      permPolyOrders(zetaRemain, upperLimit - z, ctrCoords, tempCoords);
    }
  }

}

void LinearArch::makeCoefs(double (*trueFcnRef)(InputParams const &, 
						GenericInputs const &, 
						vector<double> const &),
			   GenericInputs const &trueFcnInputsRef)
{

  /* Check for initialization. */
  if (!initialized)
  {
    cout << "Attempting to use LinearArch without initializing, please check code. " 
	 << endl;
    exit(1);
  }
  
  /* Set pointer to true function. */
  trueFcn = trueFcnRef;
  trueFcnInputs = trueFcnInputsRef;

  switch (algParams.coefsConstructionMethod)
  {
  case 1:
    /* Regression. */

    /* Sample regression data points. */
    for (int p = 0; p < algParams.nRegressionSamples; p++)
    {

      /* !!! Currently sample according to prior only, in future want
       * to sample according to the approximated state measure,
       * perhaps adaptively. */
      stateSample[0] = gsl_ran_gaussian(generator, 1.0)
	* sqrt(algParams.initialState[1]) + algParams.initialState[0];
      stateSample[1] = gsl_rng_uniform(generator)
      	* (algParams.initialState[1] - 1.0e-5) + 1.0e-5;

      /* Evaluate the features on this state, and construct the LHS
       * matrix. */
      evalAllFeatures(stateSample, LHS[p]); 
      for (unsigned int i = 0; i < ATrans.size(); i++)
	ATrans[i][p] = LHS[p][i];
      
      /* Evaluate the true function values, and construct the RHS
       * vector. */
      B[p] = trueFcn(algParams, trueFcnInputs, stateSample);
      RHS[p] = B[p];
    }

    // cout << "LHS = " << endl;
    // for (unsigned int i = 0; i < LHS.size(); i++)
    // {
    //   for (unsigned int j = 0; j < LHS[0].size(); j++)
    // 	cout << LHS[i][j] << "  ";
    //   cout << endl;
    // }
    // cout << endl;

    // for (unsigned int i = 0; i < RHS.size(); i++)
    // {
    //   RHS[i] = i;
    //   B[i] = RHS[i];
    // }

    // cout << "RHS = " << endl;
    // for (unsigned int i = 0; i < RHS.size(); i++)
    //   cout << B[i] << endl;
    // cout << endl;

    /* Solve linear least squares problem. Linear least square solve
     * verified with MATLAB. */
    linearLeastSquares(ATrans, singularValues, B, work, lwork, soln);
    
    /* Store features coefficients. */
    for (unsigned int i = 0; i < soln.size(); i++)
      featuresCoefs[i] = soln[i];

    // cout << "soln = " << endl;
    // for (unsigned int i = 0; i < soln.size(); i++)
    //   cout << i << "  " << B[i] << endl;
    
    break;
    
  default:
    cout << "Error: Approximation function coefficient construction method " 
	 << algParams.coefsConstructionMethod << " not available." << endl;
    exit(1);
  }
  
}

void LinearArch::evalAllFeatures(vector<double> const &inputVar, 
				 vector<double> &storage)
{
  
  switch (algParams.featuresChoice)
  {
  case 1:
    /* Total order polynomial. */
    for (unsigned int i = 0; i < storage.size(); i++)
      storage[i] = pow(inputVar[0], double(refTable[i][0])) 
	* pow(inputVar[1], double(refTable[i][1]));

    break;

  case 2:
    /* Terms from Gaussian KL divergence: 1.0, mean, mean-squared,
     * variance, log-variance. In that order. */
    storage[0] = 1.0;
    storage[1] = inputVar[0];
    storage[2] = inputVar[0] * inputVar[0];
    storage[3] = inputVar[1];
    storage[4] = log(inputVar[1]);
    
    break;
    
  default:
    cout << "Error: Features choice " << algParams.featuresChoice
	 << " not available." << endl;
    exit(1);
  }
  
}

double LinearArch::evalArchitecture(InputParams const &algParamsExternal, 
				    vector<double> const &state)
{

  /* Initialization. */
  double value(0.0);
  
  /* Evaluate all features on the current input state, and then sum
   * them weighed according to the coefficients. */
  evalAllFeatures(state, featuresVals);
  for (unsigned int i = 0; i < featuresCoefs.size(); i++)
    value += featuresCoefs[i] * featuresVals[i];
  
  return value;
  
}

void LinearArch::exportCoefs(vector<double> &coefsExternal)
{
  
  /* Extract the coefs to external storage. */
  for (unsigned int i = 0; i < featuresCoefs.size(); i++)
    coefsExternal[i] = featuresCoefs[i];
  
}



/*********************************************************************
 *********************************************************************/
