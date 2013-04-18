#include "linearArch.h"

/*********************************************************************
 * LinearArch
 *********************************************************************/

LinearArch::LinearArch(Controls const &refControls)
{
  
  /* Initializations. */
  primary = refControls;

  featuresCoefs = new double[primary.nFeatures];
  stateSample = new double[primary.nStatesDim];

  /* Total order polynomial features. */
  if (primary.featuresChoice == 1)
  {
    refTable = new int*[primary.nFeatures];
    refTable[0] = new int[primary.nFeatures * primary.nStatesDim];
    for (int i = 1; i < primary.nFeatures; i++)
      refTable[i] = refTable[i - 1] + primary.nStatesDim;

    /* Make refTable. */
    int ctrCoords(0);
    int * tempCoords = new int[primary.nStatesDim];
    
    permPolyOrders(primary.nStatesDim, primary.pOrder, ctrCoords, tempCoords);

    delete [] tempCoords;
  }
  
  /* Construct coefficients via construction. */
  if (primary.coefsConstructionMethod == 1)
  {
    /* Memory allocation for regression method. Note ATrans is stored
     * in its transposed form from the algebraic form of the matrix to
     * accomodate LAPACK (which uses its transpose). Note ATrans is
     * overwritten by right singular vectors upon output of linear
     * least squares computation. LHS (in non-transposed form) will be
     * used to store the original matrix. */
    ATrans = new double*[primary.nFeatures];
    LHS = new double*[primary.nRegressionSamples];
    ATrans[0] = new double[primary.nFeatures * primary.nRegressionSamples];
    LHS[0] = new double[primary.nRegressionSamples * primary.nFeatures];
    for (int i = 1; i < primary.nFeatures; i++)
      ATrans[i] = ATrans[i - 1] + primary.nRegressionSamples;
    for (int i = 1; i < primary.nRegressionSamples; i++)
      LHS[i] = LHS[i - 1] + primary.nFeatures;
    B = new double[primary.nRegressionSamples];
    RHS = new double[primary.nRegressionSamples];
    /* Size of soln array required to be at least
     * primary.nRegressionSamples for LAPACK. */
    soln = new double[primary.nFeatures];

    singularValues = new double[primary.nFeatures];
    lwork = 3 * primary.nFeatures + primary.nRegressionSamples;
    work = new double[lwork];
  }
  
  /* Random number generator initialization. */
  rngType = gsl_rng_ranlxs0;
  generator = gsl_rng_alloc(rngType);
  gsl_rng_env_setup();

}

LinearArch::~LinearArch()
{

  /* Free memory. */
  delete [] featuresCoefs;
  delete [] stateSample;

  if (primary.featuresChoice == 1)
  {
    delete [] refTable[0];
    delete [] refTable;
  }

  if (primary.coefsConstructionMethod == 1)
  {
    delete [] ATrans[0];
    delete [] ATrans;
    delete [] LHS[0];
    delete [] LHS;
    delete [] B;
    delete [] RHS;
    delete [] soln;
    delete [] singularValues;
    delete [] work;
  }
  
  gsl_rng_free(generator);

}

LinearArch::LinearArch(LinearArch const &other)
{

  /* Initializations. */
  primary = other.primary;

  featuresCoefs = new double[primary.nFeatures];
  stateSample = new double[primary.nStatesDim];

  /* Total order polynomial features. */
  if (primary.featuresChoice == 1)
  {
    refTable = new int*[primary.nFeatures];
    refTable[0] = new int[primary.nFeatures * primary.nStatesDim];
    for (int i = 1; i < primary.nFeatures; i++)
      refTable[i] = refTable[i - 1] + primary.nStatesDim;

    /* Copy the same refTable from input, assuming it is the same. */
    for (int i = 0; i < primary.nFeatures * primary.nStatesDim; i++)
      refTable[0][i] = other.refTable[0][i];
  }
  
  /* Construct coefficients via construction. */
  if (primary.coefsConstructionMethod == 1)
  {
    /* Memory allocation for regression method. */
    /* Memory allocation for regression method. Note A is stored in
     * its transposed form to accomodate LAPACK (which uses its
     * transpose). Note A is overwritten by right singular vectors
     * upon output of linear least squares computation. LHS will be
     * used to store the original matrix. */
    ATrans = new double*[primary.nFeatures];
    LHS = new double*[primary.nRegressionSamples];
    ATrans[0] = new double[primary.nFeatures * primary.nRegressionSamples];
    LHS[0] = new double[primary.nRegressionSamples * primary.nFeatures];
    for (int i = 1; i < primary.nFeatures; i++)
      ATrans[i] = ATrans[i - 1] + primary.nRegressionSamples;
    for (int i = 1; i < primary.nRegressionSamples; i++)
      LHS[i] = LHS[i - 1] + primary.nFeatures;
    B = new double[primary.nRegressionSamples];
    RHS = new double[primary.nRegressionSamples];
    /* Size of soln array required to be at least
     * primary.nRegressionSamples for LAPACK. */
    soln = new double[primary.nFeatures];

    singularValues = new double[primary.nFeatures];
    lwork = 3 * primary.nFeatures + primary.nRegressionSamples;
    work = new double[lwork];
  }

  /* Random number generator initialization. */
  rngType = gsl_rng_ranlxs0;
  generator = gsl_rng_alloc(rngType);
  gsl_rng_env_setup();  
  
}

LinearArch& LinearArch::operator=(LinearArch const &rhs)
{

  /* Protect against invalid self-assignment. */
  if (this != &rhs)
  {

    /* Initializations. */
    primary = rhs.primary;

    featuresCoefs = new double[primary.nFeatures];
    stateSample = new double[primary.nStatesDim];

    /* Total order polynomial features. */
    if (primary.featuresChoice == 1)
    {
      refTable = new int*[primary.nFeatures];
      refTable[0] = new int[primary.nFeatures * primary.nStatesDim];
      for (int i = 1; i < primary.nFeatures; i++)
	refTable[i] = refTable[i - 1] + primary.nStatesDim;

      /* Copy the same refTable from input, assuming it is the same. */
      for (int i = 0; i < primary.nFeatures * primary.nStatesDim; i++)
	refTable[0][i] = rhs.refTable[0][i];
    }
  
    /* Construct coefficients via construction. */
    if (primary.coefsConstructionMethod == 1)
    {
      /* Memory allocation for regression method. */
      /* Memory allocation for regression method. Note A is stored in
       * its transposed form to accomodate LAPACK (which uses its
       * transpose). Note A is overwritten by right singular vectors
       * upon output of linear least squares computation. LHS will be
       * used to store the original matrix. */
      ATrans = new double*[primary.nFeatures];
      LHS = new double*[primary.nRegressionSamples];
      ATrans[0] = new double[primary.nFeatures * primary.nRegressionSamples];
      LHS[0] = new double[primary.nRegressionSamples * primary.nFeatures];
      for (int i = 1; i < primary.nFeatures; i++)
	ATrans[i] = ATrans[i - 1] + primary.nRegressionSamples;
      for (int i = 1; i < primary.nRegressionSamples; i++)
	LHS[i] = LHS[i - 1] + primary.nFeatures;
      B = new double[primary.nRegressionSamples];
      RHS = new double[primary.nRegressionSamples];
      /* Size of soln array required to be at least
       * primary.nRegressionSamples for LAPACK. */
      soln = new double[primary.nFeatures];

      singularValues = new double[primary.nFeatures];
      lwork = 3 * primary.nFeatures + primary.nRegressionSamples;
      work = new double[lwork];
    }
    
    /* Random number generator initialization. */
    rngType = gsl_rng_ranlxs0;
    generator = gsl_rng_alloc(rngType);
    gsl_rng_env_setup();  

  }
  
  /* By convention, always return *this. */
  return *this;
  
}

void LinearArch::permPolyOrders(int zetaRemain, int const upperLimit, 
				int &ctrCoords, int* const tempCoords)
{

  /* Base case. */
  if (zetaRemain == 1)
  {
    for (int z = 0; z < upperLimit + 1; z++)
    {
      tempCoords[primary.nStatesDim - 1] = z;
      for (int i = 0; i < primary.nStatesDim; i++)
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
      tempCoords[primary.nStatesDim - 1 - zetaRemain] = z;
      permPolyOrders(zetaRemain, upperLimit - z, ctrCoords, tempCoords);
    }
  }

}

void LinearArch::makeCoefs(double (*trueFcnRef)(Controls const &, GenericInputs &, 
						double const * const),
			   GenericInputs &trueFcnInputsRef)
{

  /* Set pointer to true function. */
  trueFcn = trueFcnRef;
  trueFcnInputs = trueFcnInputsRef;

  switch (primary.coefsConstructionMethod)
  {
  case 1:
    /* Regression. */

    /* Sample regression data points. */
    for (int p = 0; p < primary.nRegressionSamples; p++)
    {

      /* !!! Currently sample according to prior only, in future want
       * to sample according to the approximated state measure,
       * perhaps adaptively. */
      stateSample[0] = gsl_ran_gaussian(generator, 1.0) 
	* sqrt(primary.initialState[1]) 
	+ primary.initialState[0];
      stateSample[1] = gsl_rng_uniform(generator)
      	* (primary.initialState[1] - 1.0e-5)
      	+ 1.0e-5;

      /* Evaluate the features on this state, and construct the LHS
       * matrix. */
      evalAllFeatures(stateSample, LHS[p]); 
      for (int i = 0; i < primary.nFeatures; i++)
	ATrans[i][p] = LHS[p][i];
      
      /* Evaluate the true function values, and construct the RHS
       * vector. */
      B[p] = trueFcn(primary, trueFcnInputs, stateSample);
      RHS[p] = B[p];
    }

    for (int i = 0; i < primary.nRegressionSamples; i++)
    {
      for (int j = 0; j < primary.nFeatures; j++)
	cout << LHS[i][j] << "  ";
      cout << endl;
    }
    cout << endl;

    for (int i = 0; i < primary.nRegressionSamples; i++)
      cout << RHS[i] << endl;
    cout << endl;

    /* Solve linear least squares problem. Linear least square solve
     * verified with MATLAB. */
    linearLeastSquares(primary.nRegressionSamples, primary.nFeatures, 
    		       ATrans, singularValues, B, work, lwork, soln);
    
    /* Store features coefficients. */
    for (int i = 0; i < primary.nFeatures; i++)
      featuresCoefs[i] = soln[i];
    
    break;
    
  default:
    cout << "Error: Approximation function coefficient construction method " 
	 << primary.coefsConstructionMethod << " not available." << endl;
    exit(1);
  }
  
}

void LinearArch::evalAllFeatures(double const * const inputVar, 
				 double * const storage)
{
  
  switch (primary.featuresChoice)
  {
  case 1:
    /* Total order polynomial. */
    for (int i = 0; i < primary.nFeatures; i++)
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
    cout << "Error: Features choice " << primary.featuresChoice
	 << " not available." << endl;
    exit(1);
  }
  
}

double LinearArch::evalArchitecture(Controls const &primaryExternal, 
				    double const * const state)
{

  /* Initialization. */
  double value(0.0);

  /* Evaluate all features on the current input state, and then sum
   * them weighed according to the coefficients. */
  evalAllFeatures(state, featuresVals);
  for (int i = 0; i < primary.nFeatures; i++)
    value += featuresCoefs[i] * featuresVals[i];
  
  return value;
  
}

void LinearArch::exportCoefs(double* const coefsExternal)
{
  
  /* Extract the coefs to external storage. */
  for (int i = 0; i < primary.nFeatures; i++)
    coefsExternal[i] = featuresCoefs[i];
  
}



/*********************************************************************
 *********************************************************************/
