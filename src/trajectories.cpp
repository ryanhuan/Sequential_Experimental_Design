#include "trajectories.h"

/*********************************************************************
 * Trajectories
 *********************************************************************/

Trajectories::Trajectories()
{
  initialized = 0;
}

Trajectories::Trajectories(InputParams const &refAlgParams)
{
  initialize(refAlgParams);
}

Trajectories::Trajectories(Trajectories const &other)
{

  /* Initializations. */
  algParams = other.algParams;

  /* Allocate space for master list of trajectories for master node
   * only. */
  if (algParams.rank == 0)
  {
    thetaSamplesAll = new double*[algParams.nTrajectories];
    thetaSamplesAll[0] = new double[algParams.nTrajectories 
				    * algParams.nInferenceParamsDim];
    statesAll = new double*[algParams.nTrajectories];
    statesAll[0] = new double[algParams.nTrajectories 
			      * (algParams.nStages + 1) * algParams.nStatesDim];
    controlsAll = new double*[algParams.nTrajectories];
    controlsAll[0] = new double[algParams.nTrajectories 
				* algParams.nStages * algParams.nControlsDim];
    disturbanceAll = new double*[algParams.nTrajectories];
    disturbanceAll[0] = new double[algParams.nTrajectories
				   * algParams.nStages * algParams.nDisturbanceDim];
    finalRewardsAll = new double[algParams.nTrajectories];
    for (int i = 1; i < algParams.nTrajectories; i++)
    {
      thetaSamplesAll[i] = thetaSamplesAll[i - 1] + algParams.nInferenceParamsDim;
      statesAll[i] = statesAll[i - 1] + (algParams.nStages + 1) * algParams.nStatesDim;
      controlsAll[i] = controlsAll[i - 1] + algParams.nStages * algParams.nControlsDim;
      disturbanceAll[i] = disturbanceAll[i - 1] 
	+ algParams.nStages * algParams.nDisturbanceDim;
    }

    /* Copy over the data. */
    for (int i = 0; i < algParams.nTrajectories * algParams.nInferenceParamsDim; i++)
      thetaSamplesAll[0][i] = other.thetaSamplesAll[0][i];
    for (int i = 0; i < algParams.nTrajectories * (algParams.nStages + 1) 
	   * algParams.nStatesDim; i++)
      statesAll[0][i] = other.statesAll[0][i];
    for (int i = 0; i < algParams.nTrajectories * algParams.nStages 
	   * algParams.nControlsDim; i++)
      controlsAll[0][i] = other.controlsAll[0][i];
    for (int i = 0; i < algParams.nTrajectories * algParams.nStages 
	   * algParams.nDisturbanceDim; i++)
      disturbanceAll[0][i] = other.disturbanceAll[0][i];
    for (int i = 0; i < algParams.nTrajectories; i++)
      finalRewardsAll[i] = other.finalRewardsAll[i];
  }
  
  /* Compute local task allocation. */
  nLocalTrajectoriesAll = other.nLocalTrajectoriesAll;
  nLocalTrajectoriesAllSum = other.nLocalTrajectoriesAllSum;
  nLocalTrajectories = nLocalTrajectoriesAll[algParams.rank];

  /* Allocate space for local list of trajectories. */
  thetaSamples = new double*[nLocalTrajectories];
  thetaSamples[0] = new double[nLocalTrajectories * algParams.nInferenceParamsDim];
  states = new double*[nLocalTrajectories];
  states[0] = new double[nLocalTrajectories 
			 * (algParams.nStages + 1) * algParams.nStatesDim];
  controls = new double*[nLocalTrajectories];
  controls[0] = new double[nLocalTrajectories 
			   * algParams.nStages * algParams.nControlsDim];
  disturbance = new double*[nLocalTrajectories];
  disturbance[0] = new double[nLocalTrajectories 
			      * algParams.nStages * algParams.nDisturbanceDim];
  finalRewards = new double[nLocalTrajectories];
  for (int i = 1; i < nLocalTrajectories; i++)
  {
    thetaSamples[i] = thetaSamples[i - 1] + algParams.nInferenceParamsDim;
    states[i] = states[i - 1] + (algParams.nStages + 1) * algParams.nStatesDim;
    controls[i] = controls[i - 1] + algParams.nStages * algParams.nControlsDim;
    disturbance[i] = disturbance[i - 1] + algParams.nStages * algParams.nDisturbanceDim;
  }

  /* Copy over the data. */
  for (int i = 0; i < nLocalTrajectories * algParams.nInferenceParamsDim; i++)
    thetaSamples[0][i] = other.thetaSamples[0][i];
  for (int i = 0; i < nLocalTrajectories * (algParams.nStages + 1) 
	 * algParams.nStatesDim; i++)
    states[0][i] = other.states[0][i];
  for (int i = 0; i < nLocalTrajectories * algParams.nStages 
	 * algParams.nControlsDim; i++)
    controls[0][i] = other.controls[0][i];
  for (int i = 0; i < nLocalTrajectories * algParams.nStages 
	   * algParams.nDisturbanceDim; i++)
    disturbance[0][i] = other.disturbance[0][i];
  for (int i = 0; i < nLocalTrajectories; i++)
    finalRewards[i] = other.finalRewards[i];
  tempTheta = other.tempTheta;
  tempState = other.tempState;
  tempControl = other.tempControl;
  tempDisturbance = other.tempDisturbance;
  tempNewState = other.tempNewState;
  
  /* Random number generator initialization. */
  rngType = gsl_rng_ranlxs0;
  generator = gsl_rng_alloc(rngType);
  gsl_rng_env_setup();
  gsl_rng_set(generator, rand() + algParams.rank);
  
  initialized = other.initialized;

}

Trajectories::~Trajectories()
{

  if (initialized)
  {
    /* Free memory. */
    if (algParams.rank == 0)
    {
      delete [] thetaSamplesAll[0];
      delete [] thetaSamplesAll;
      delete [] statesAll[0];
      delete [] statesAll;
      delete [] controlsAll[0];
      delete [] controlsAll;
      delete [] disturbanceAll[0];
      delete [] disturbanceAll;
      delete [] finalRewardsAll;
    }
    delete [] thetaSamples[0];
    delete [] thetaSamples;
    delete [] states[0];
    delete [] states;
    delete [] controls[0];
    delete [] controls;
    delete [] disturbance[0];
    delete [] disturbance;
    delete [] finalRewards;
    gsl_rng_free(generator);
  }
  
}

Trajectories& Trajectories::operator=(Trajectories const &rhs)
{

  /* Protect against invalid self-assignment. */
  if (this != &rhs)
  {

    /* Initializations. */
    algParams = rhs.algParams;

    /* Master list of trajectories for master node only. */
    if (algParams.rank == 0)
    {
      /* Copy over the data. */
      for (int i = 0; i < algParams.nTrajectories * algParams.nInferenceParamsDim; i++)
	thetaSamplesAll[0][i] = rhs.thetaSamplesAll[0][i];
      for (int i = 0; i < algParams.nTrajectories * (algParams.nStages + 1) 
	     * algParams.nStatesDim; i++)
	statesAll[0][i] = rhs.statesAll[0][i];
      for (int i = 0; i < algParams.nTrajectories * algParams.nStages 
	     * algParams.nControlsDim; i++)
	controlsAll[0][i] = rhs.controlsAll[0][i];
      for (int i = 0; i < algParams.nTrajectories * algParams.nStages 
	     * algParams.nDisturbanceDim; i++)
	disturbanceAll[0][i] = rhs.disturbanceAll[0][i];
      for (int i = 0; i < algParams.nTrajectories; i++)
	finalRewardsAll[i] = rhs.finalRewardsAll[i];
    }
  
    /* Compute local task allocation. */
    nLocalTrajectoriesAll = rhs.nLocalTrajectoriesAll;
    nLocalTrajectoriesAllSum = rhs.nLocalTrajectoriesAllSum;
    nLocalTrajectories = nLocalTrajectoriesAll[algParams.rank];

    /* Copy over the data. */
    for (int i = 0; i < nLocalTrajectories * algParams.nInferenceParamsDim; i++)
      thetaSamples[0][i] = rhs.thetaSamples[0][i];
    for (int i = 0; i < nLocalTrajectories * (algParams.nStages + 1) 
	   * algParams.nStatesDim; i++)
      states[0][i] = rhs.states[0][i];
    for (int i = 0; i < nLocalTrajectories * algParams.nStages 
	   * algParams.nControlsDim; i++)
      controls[0][i] = rhs.controls[0][i];
    for (int i = 0; i < nLocalTrajectories * algParams.nStages 
	   * algParams.nDisturbanceDim; i++)
      disturbance[0][i] = rhs.disturbance[0][i];
    for (int i = 0; i < nLocalTrajectories; i++)
      finalRewards[i] = rhs.finalRewards[i];
    tempTheta = rhs.tempTheta;
    tempState = rhs.tempState;
    tempControl = rhs.tempControl;
    tempDisturbance = rhs.tempDisturbance;
    tempNewState = rhs.tempNewState;
  
    /* Random number generator initialization. */
    rngType = gsl_rng_ranlxs0;
    generator = gsl_rng_alloc(rngType);
    gsl_rng_env_setup();
    gsl_rng_set(generator, rand() + algParams.rank);

    initialized = rhs.initialized;

  }
  
  /* By convention, always return *this. */
  return *this;
  
}

void Trajectories::initialize(InputParams const &refAlgParams)
{
  
  /* Initializations. */
  algParams = refAlgParams;

  /* Allocate space for master list of trajectories for master node
   * only. */
  if (algParams.rank == 0)
  {
    thetaSamplesAll = new double*[algParams.nTrajectories];
    thetaSamplesAll[0] = new double[algParams.nTrajectories 
				    * algParams.nInferenceParamsDim];
    statesAll = new double*[algParams.nTrajectories];
    statesAll[0] = new double[algParams.nTrajectories 
			      * (algParams.nStages + 1) * algParams.nStatesDim];
    controlsAll = new double*[algParams.nTrajectories];
    controlsAll[0] = new double[algParams.nTrajectories 
				* algParams.nStages * algParams.nControlsDim];
    disturbanceAll = new double*[algParams.nTrajectories];
    disturbanceAll[0] = new double[algParams.nTrajectories
				   * algParams.nStages * algParams.nDisturbanceDim];
    finalRewardsAll = new double[algParams.nTrajectories];
    for (int i = 1; i < algParams.nTrajectories; i++)
    {
      thetaSamplesAll[i] = thetaSamplesAll[i - 1] + algParams.nInferenceParamsDim;
      statesAll[i] = statesAll[i - 1] + (algParams.nStages + 1) * algParams.nStatesDim;
      controlsAll[i] = controlsAll[i - 1] + algParams.nStages * algParams.nControlsDim;
      disturbanceAll[i] = disturbanceAll[i - 1] 
	+ algParams.nStages * algParams.nDisturbanceDim;
    }
  }

  /* Compute local task allocation. */
  nLocalTrajectoriesAll.assign(algParams.nTasks, 0);
  nLocalTrajectoriesAllSum.assign(algParams.nTasks, 0);
  int base = algParams.nTrajectories / algParams.nTasks;
  unsigned int nExtras = algParams.nTrajectories % algParams.nTasks;
  for (unsigned int i = 0; i < nLocalTrajectoriesAll.size(); i++)
  {
    /* Assign local trajectory numbers. */
    if (i < nExtras)
      nLocalTrajectoriesAll[i] = base + 1;
    else
      nLocalTrajectoriesAll[i] = base;

    /* Summation. */
    if (i == 0)
      nLocalTrajectoriesAllSum[i] = nLocalTrajectoriesAll[i];
    else
      nLocalTrajectoriesAllSum[i] = nLocalTrajectoriesAllSum[i - 1] + nLocalTrajectoriesAll[i];
  }
  if (nLocalTrajectoriesAllSum[nLocalTrajectoriesAllSum.size() - 1] != algParams.nTrajectories)
  {
    cout << "Error: check local task allocation for trajectories simulation." << endl;
    exit(1);
  }
  nLocalTrajectories = nLocalTrajectoriesAll[algParams.rank];

  /* Allocate space for local list of trajectories. */
  thetaSamples = new double*[nLocalTrajectories];
  thetaSamples[0] = new double[nLocalTrajectories * algParams.nInferenceParamsDim];
  states = new double*[nLocalTrajectories];
  states[0] = new double[nLocalTrajectories 
			 * (algParams.nStages + 1) * algParams.nStatesDim];
  controls = new double*[nLocalTrajectories];
  controls[0] = new double[nLocalTrajectories 
			   * algParams.nStages * algParams.nControlsDim];
  disturbance = new double*[nLocalTrajectories];
  disturbance[0] = new double[nLocalTrajectories 
			      * algParams.nStages * algParams.nDisturbanceDim];
  finalRewards = new double[nLocalTrajectories];
  for (int i = 1; i < nLocalTrajectories; i++)
  {
    thetaSamples[i] = thetaSamples[i - 1] + algParams.nInferenceParamsDim;
    states[i] = states[i - 1] + (algParams.nStages + 1) * algParams.nStatesDim;
    controls[i] = controls[i - 1] + algParams.nStages * algParams.nControlsDim;
    disturbance[i] = disturbance[i - 1] + algParams.nStages * algParams.nDisturbanceDim;
  }
  tempTheta.assign(algParams.nInferenceParamsDim, 0.0);
  tempState.assign(algParams.nStatesDim, 0.0);
  tempControl.assign(algParams.nControlsDim, 0.0);
  tempDisturbance.assign(algParams.nDisturbanceDim, 0.0);
  tempNewState.assign(algParams.nStatesDim, 0.0);
  
  /* Random number generator initialization. */
  rngType = gsl_rng_ranlxs0;
  generator = gsl_rng_alloc(rngType);
  gsl_rng_env_setup();
  gsl_rng_set(generator, rand() + algParams.rank);
  
  initialized = 1;
  
}

void Trajectories::simulateTrajectories(vector<LinearArch> &arch)
{

  if (!initialized)
  {
    cout << "Error: Attempting to use Trajectories without initializing, please check code." 
	 << endl;
    exit(1);
  }

  /* Initializations. */
  GenericInputs maxExpInputs;

  for (int t = 0; t < nLocalTrajectories; t++)
  {
    
    /* Sample local parameter instances from the prior. */
    //!!! For future need to generalize to non-Gaussians. 
    for (int i = 0; i < algParams.nInferenceParamsDim; i++)
    {
      thetaSamples[t][i] = algParams.initialState[0]
	+ sqrt(algParams.initialState[1]) * gsl_ran_gaussian(generator, 1.0);
      tempTheta[i] = thetaSamples[t][i];
    }

    /* Reset initial state. */
    finalRewards[t] = 0.0;
    for (int i = 0; i < algParams.nStatesDim; i++)
      states[t][i] = algParams.initialState[i];
    
    /* Evaluate the cost-to-go functions at the realized states. */
    for (int k = 0; k < algParams.nStages; k++)
    {
      
      if (k == (algParams.nStages - 1))
      {
	/* For when J_{k+1} is terminal reward. */
	maxExpInputs.systemEqnPtr = &systemEquation;
	maxExpInputs.stageFcnPtr = &stageCost;
	maxExpInputs.futureFcnPtr = &evalTerminalReward;
	maxExpInputs.futureLinearArchClassPtr = NULL;
	for (unsigned int i = 0; i < tempState.size(); i++)
	  tempState[i] = states[t][k * algParams.nStatesDim + i];
      }
      else
      {
	/* For when J_{k+1} is the previous approximation function. */
	maxExpInputs.systemEqnPtr = &systemEquation;
	maxExpInputs.stageFcnPtr = &stageCost;
	maxExpInputs.futureFcnPtr = NULL;
	maxExpInputs.futureLinearArchClassPtr = &(arch[k]);
	for (unsigned int i = 0; i < tempState.size(); i++)
	  tempState[i] = states[t][k * algParams.nStatesDim + i];
      }	

      /* Extract and record optimal control. */
      maxExpectation(algParams, maxExpInputs, tempState, tempControl);
      for (unsigned int i = 0; i < tempControl.size(); i++)
	controls[t][k * algParams.nControlsDim + i] = tempControl[i];
      
      /* Generate and record disturbance. */
      generateDisturbance(algParams, tempTheta, tempState, tempControl, generator, 
			  tempDisturbance);
      for (unsigned int i = 0; i < tempDisturbance.size(); i++)
	disturbance[t][k * algParams.nDisturbanceDim + i] = tempDisturbance[i];
      
      /* Evaluate system equation to compute state at stage k+1. */
      systemEquation(algParams, tempState, tempControl, tempDisturbance, 
		     tempNewState);
      for (unsigned int i = 0; i < tempNewState.size(); i++)
	states[t][(k + 1) * algParams.nStatesDim + i] = tempNewState[i];
      
      if (k == (algParams.nStages - 1))
	/* Evaluate and accumulate terminal reward. */
	finalRewards[t] += evalTerminalReward(algParams, tempState);
      else
	/* Evaluate and accumulate stage cost. */
	finalRewards[t] -= stageCost(algParams, tempState, tempControl, 
				     tempDisturbance);
      
    }
    
    /* Collect results on master node. */
    if (algParams.rank != 0)
    {
      MPI_Send(thetaSamples[0], sizeof(double) * nLocalTrajectories 
	       * algParams.nInferenceParamsDim, 
	       MPI_CHAR, 0, algParams.rank, MPI_COMM_WORLD);

      MPI_Send(states[0], sizeof(double) * nLocalTrajectories 
	       * (algParams.nStages + 1) * algParams.nStatesDim, 
	       MPI_CHAR, 0, algParams.rank + algParams.nTasks, MPI_COMM_WORLD);

      MPI_Send(controls[0], sizeof(double) * nLocalTrajectories 
	       * algParams.nStages * algParams.nControlsDim, 
	       MPI_CHAR, 0, algParams.rank + 2 * algParams.nTasks, MPI_COMM_WORLD);

      MPI_Send(disturbance[0], sizeof(double) * nLocalTrajectories 
	       * algParams.nStages * algParams.nDisturbanceDim, 
	       MPI_CHAR, 0, algParams.rank + 3 * algParams.nTasks, MPI_COMM_WORLD);

      MPI_Send(finalRewards, sizeof(double) * nLocalTrajectories, 
	       MPI_CHAR, 0, algParams.rank + 4 * algParams.nTasks, MPI_COMM_WORLD);
    }
    else
    {
      
      /* Store local trajectories computed on master node. */
      for (int i = 0; i < nLocalTrajectories * algParams.nInferenceParamsDim; i++)
	thetaSamplesAll[0][i] = thetaSamples[0][i];
      for (int i = 0; i < nLocalTrajectories * (algParams.nStages + 1) 
	     * algParams.nStatesDim; i++)
	statesAll[0][i] = states[0][i];
      for (int i = 0; i < nLocalTrajectories * algParams.nStages 
	     * algParams.nControlsDim; i++)
	controlsAll[0][i] = controls[0][i];
      for (int i = 0; i < nLocalTrajectories * algParams.nStages 
	     * algParams.nDisturbanceDim; i++)
	disturbanceAll[0][i] = disturbance[0][i];
      for (int i = 0; i < nLocalTrajectories; i++)
	finalRewardsAll[i] = finalRewards[i];

      /* Receive and store trajectories computed from other nodes. */
      for (int r = 1; r < algParams.nTasks; r++)
      {
	MPI_Recv(thetaSamplesAll[nLocalTrajectoriesAllSum[r - 1]], 
		 sizeof(double) * nLocalTrajectoriesAll[r] 
		 * algParams.nInferenceParamsDim, MPI_CHAR, 
		 r, r, MPI_COMM_WORLD, &algParams.status);

	MPI_Recv(statesAll[nLocalTrajectoriesAllSum[r - 1]], 
		 sizeof(double) * nLocalTrajectoriesAll[r] 
		 * (algParams.nStages + 1) * algParams.nStatesDim, MPI_CHAR, 
		 r, r + algParams.nTasks, MPI_COMM_WORLD, &algParams.status);

	MPI_Recv(controlsAll[nLocalTrajectoriesAllSum[r - 1]], 
		 sizeof(double) * nLocalTrajectoriesAll[r] 
		 * algParams.nStages * algParams.nControlsDim, MPI_CHAR, 
		 r, r + 2 * algParams.nTasks, MPI_COMM_WORLD, &algParams.status);

	MPI_Recv(disturbanceAll[nLocalTrajectoriesAllSum[r - 1]], 
		 sizeof(double) * nLocalTrajectoriesAll[r] 
		 * algParams.nStages * algParams.nDisturbanceDim, MPI_CHAR, 
		 r, r + 3 * algParams.nTasks, MPI_COMM_WORLD, &algParams.status);

	MPI_Recv(&(finalRewardsAll[nLocalTrajectoriesAllSum[r - 1]]), 
		 sizeof(double) * nLocalTrajectoriesAll[r], MPI_CHAR, 
		 r, r + 4 * algParams.nTasks, MPI_COMM_WORLD, &algParams.status);
      }
      
    }
  }
  
}

void Trajectories::writeTrajectoriesToFile(string const fName)
{

  /* Do this on master node only. */
  if (algParams.rank == 0)
  {
    
    /* Open file. */
    ofstream f;
    
    /* Write as text file. */
    f.open(fName.c_str());
  
    /* Set output format. */
    f.setf(ios::scientific,ios::floatfield);
    f.precision(16);

    /* Print labels. */
    f << "% [True Theta Sample] + [x_k, u_k, w_k] for k=1...N-1 + [x_N, final cost] " << endl;
    
    for (int t = 0; t < algParams.nTrajectories; t++)
    {
      /* Write theta sample for that run. */
      for (int i = 0; i < algParams.nInferenceParamsDim; i++)
	f << thetaSamplesAll[t][i] << "  ";

      /* Write all intermediate states, controls, and disturbance. */
      for (int k = 0; k < algParams.nStages; k++)
      {
	for (int i = 0; i < algParams.nStatesDim; i++)
	  f << statesAll[t][k * algParams.nStatesDim + i] << "  ";
	for (int i = 0; i < algParams.nControlsDim; i++)
	  f << controlsAll[t][k * algParams.nControlsDim + i] << "  ";
	for (int i = 0; i < algParams.nDisturbanceDim; i++)
	  f << disturbanceAll[t][k * algParams.nDisturbanceDim + i] << "  ";
      }
      
      /* Write final states and rewards. */
      for (int i = 0; i < algParams.nStatesDim; i++)
	f << statesAll[t][algParams.nStages * algParams.nStatesDim + i] << "  ";

      f << finalRewardsAll[t] << endl;
    }
  
    /* Close file. */
    f.close();

  }

}



/*********************************************************************
 *********************************************************************/
