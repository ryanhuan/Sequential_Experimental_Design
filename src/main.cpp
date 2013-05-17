#include "main.h"

int main(int argc, char **argv)
{

  /* Start time. */
  clock_t t0, t1;
  t0 = time(NULL);
  
  /* MPI initializations. */
  int nTasks, rank, rc; 
  rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS)
  {
    cout << "Error: MPI program failed to start. " << endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  /* Construct and set input parameter values from input file. */
  InputParams algParams("inputParams.inp");
  algParams.addMPIInfo(nTasks, rank);

  /* Linear architecture initialization. */
  vector<double> tempArchCoefs(algParams.nFeatures, 0.0);
  vector<LinearArch> arch(algParams.nStages - 1, LinearArch(algParams));
  GenericInputs maxExpInputs;
  string fName;
  bool needMakeCoefs(0);

  if (algParams.simulateTrajectories)
  {
    /* Read in coefficients from files. */
    for (unsigned int k = arch.size(); k > 0; k--)
    {
      fName = string("TJ_") + num2string<int>(k) + string("_coefs.dat");
      if (arch[k - 1].readCoefsFromFile(fName))
      {
	needMakeCoefs = 1;
	break;
      }
    }
  }
  else
    needMakeCoefs = 1;
  
  /* Make coefficients. */
  if (needMakeCoefs)
  {
    
    /* DP algorithm (backward induction), first construct all the
     * surrogate functions using linear architecture. */
    for (unsigned int k = arch.size(); k > 0; k--)
    {
      if (algParams.rank == 0)
	cout << "Making TJ_" << k << "... " << endl;
    
      //!!! for now, we are going to use the exact function for J2, ie
      //!!! no approximation. In future, this would be needed for
      //!!! general, non-Gaussian distributions. Note the first if
      //!!! condition is approximating the MAX EXPECTED J2.
      /* Pass in function pointer to the true function to be
       * approximate. */
      if (k == arch.size())
      {
	/* For when true function is terminal reward. */
	maxExpInputs.systemEqnPtr = &systemEquation;
	maxExpInputs.stageFcnPtr = &stageCost;
	maxExpInputs.futureFcnPtr = &evalTerminalReward;
	maxExpInputs.futureLinearArchClassPtr = NULL;
	arch[k - 1].makeCoefs(maxExpectation, maxExpInputs);
      }
      else
      {
	/* For when true function is the previous approximation function. */
	maxExpInputs.systemEqnPtr = &systemEquation;
	maxExpInputs.stageFcnPtr = &stageCost;
	maxExpInputs.futureFcnPtr = NULL;
	maxExpInputs.futureLinearArchClassPtr = &(arch[k]);
	arch[k - 1].makeCoefs(maxExpectation, maxExpInputs);
      }

      /* Retrieve and output coefficients. */
      arch[k - 1].exportCoefs(tempArchCoefs);
      // for (unsigned int i = 0; i < tempArchCoefs.size(); i++)
      //   cout << tempArchCoefs[i] << endl;

      /* Write coefficients to file. */
      fName = string("TJ_") + num2string<int>(k) + string("_coefs.dat");
      arch[k - 1].writeCoefsToFile(fName);    
    }

    if (algParams.rank == 0)
      cout << "Done making all TJ_k." << endl;

  }

  /* Simulate trajectories and write them to file. */
  if (algParams.simulateTrajectories)
  {
    Trajectories traj(algParams);
    traj.simulateTrajectories(arch);
    fName = string("trajectories.dat");
    traj.writeTrajectoriesToFile(fName);
  }

  /* Finalize MPI. */
  MPI_Finalize();

  /* Computing time. */
  t1 = time(NULL);
  if (algParams.rank == 0)
    cout << "% Estimated wall clock time = " << t1-t0 << " s." << endl;
  
  return 0;

}

