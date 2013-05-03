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
  
  /* Construct and set input parameter values from input file. */
  InputParams algParams("inputParams.inp");
  algParams.addMPIInfo(nTasks, rank);

  /* Linear architecture initialization. */
  vector<double> tempArchCoefs(algParams.nFeatures, 0.0);
  vector<LinearArch> arch(algParams.nStages - 1, LinearArch(algParams));
  GenericInputs maxExpInputs;
  
  /* DP algorithm (backward induction), first construct all the
   * surrogate functions using linear architecture. */
  for (unsigned int k = arch.size(); k > 0; k--)
  {
    if (algParams.rank == 0)
      cout << "Making TJ_" << k << ". " << endl;
    
    //!!! for now, we are going to use the exact function for J2, ie
    //!!! no approximation. In future, this would be needed for
    //!!! general, non-Gaussian distributions. Note the first if
    //!!! condition is approximating the MAX EXPECTED J2.
    /* Pass in function pointer to the true function to be
     * approximate. */
    if (k == arch.size())
    {
      /* For when true function is terminal reward. */
      maxExpInputs.systemEqnPtr = systemEquation;
      maxExpInputs.stageFcnPtr = stageCost;
      maxExpInputs.futureFcnPtr = evalTerminalReward;
      arch[k - 1].makeCoefs(maxExpectation, maxExpInputs);
    }
    else
    {
      /* For when true function is the previous approximation function. */
      //!!! this does not work yet, since our pointer is for pointing
      //!!! to ordinary functions and not class member functions.

      // maxExpInputs.systemEqnPtr = systemEquation;
      // maxExpInputs.stageFcnPtr = stageCost;
      // maxExpInputs.futureFcnPtr = arch[k].evalArchitecture;
      // arch[k - 1].makeCoefs(maxExpectation, maxExpInputs);
    }

    arch[k - 1].exportCoefs(tempArchCoefs);

    for (unsigned int i = 0; i < tempArchCoefs.size(); i++)
      cout << tempArchCoefs[i] << endl;
  }

  /* Finalize MPI. */
  MPI_Finalize();

  return 0;

}

