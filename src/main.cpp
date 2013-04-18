#include "main.h"

int main(int argc, char **argv)
{
  
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
  
  /* Construct and set control values. Make changes through input
   * file. */
  Controls primary;
  initializeControls(primary);
  addMPIInfoControls(primary, nTasks, rank);

  /* Linear architecture initialization. */
  vector<double> tempArchCoefs(primary.nFeatures, 0.0);
  vector<LinearArch> arch(primary.nStages - 1, primary);
  GenericInputs maxExpInputs;
  
  /* DP algorithm (backward induction), first construct all the
   * surrogate functions using linear architecture. */
  for (int k = arch.size() - 1; k > -1; k--)
  {
    if (primary.rank == 0)
      cout << "Making TJ_" << k + 1 << ". " << endl;
    
    //!!! for now, we are going to use the exact function for J2, ie
    //!!! no approximation. In future, this would be needed for
    //!!! general, non-Gaussian distributions. Note the first if
    //!!! condition is approximating the MAX EXPECTED J2.
    /* Pass in function pointer to the true function to be
     * approximate. */
    if (k == int(arch.size() - 1))
    {
      /* For when true function is terminal reward. */
      maxExpInputs.futureFcnPtr = evalTerminalRewardFromKm1;
      maxExpInputs.stageFcnPtr = quadraticControlsCost;
      arch[k].makeCoefs(maxExpectation, maxExpInputs);
    }
    else
    {
      /* For when true function is the previous approximation function. */
      maxExpInputs.futureFcnPtr = arch[k+1].evalArchitecture;
      maxExpInputs.stageFcnPtr = quadraticControlsCost;
      arch[k].makeCoefs(maxExpectation, maxExpInputs);
    }

    arch[k].exportCoefs(tempArchCoefs);

    for (int i = 0; i < tempArchCoefs.size(); i++)
      cout << tempArchCoefs[i] << endl;
  }
  
  /* Free memory. */
  finalizeControls(primary);

  /* Finalize MPI. */
  MPI_Finalize();

  return 0;

}

