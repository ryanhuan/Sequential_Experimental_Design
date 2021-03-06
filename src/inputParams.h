/*! \file inputParams.h 

  \brief Class and functions concerning the input parameters of the
  algorithm.

*/
#ifndef _INPUTPARAMS_H
#define _INPUTPARAMS_H

#include <iostream>
#include <string.h>
#include <vector>

#include "mpi.h"

#include "fileRead.h"
#include "tools.h"

using namespace std;

/*! \class InputParams

  \brief Class for algorithm input parameters.  

  This class stores the input and control parameters for this
  algorithm. All members are public. All copies will be by value since
  we will not use pointers as members in this class.

*/
class InputParams
{
  
public:
  
  /* General. */
  double doubleCompTol;         //!< Comparison tolerance for doubles.

  /* MPI. */
  int nTasks;                   //!< Total number of CPUs.
  int rank;                     //!< Current CPU number.
  MPI_Status status;            //!< MPI status variable.

  /* Dynamic programming. */
  int nStages;                  //!< Horizon length.
  int nStatesDim;               //!< State space dimension.
  int nControlsDim;             //!< Control space dimension.
  int nDisturbanceDim;          //!< Disturbance space dimension.
  int nInferenceParamsDim;      //!< Inference parameters space dimension. Ensure consistent with nStatesDim.
  int nSearchDim;               //!< Search space dimension for optimization (same as nControlsDim in this case). Computed. 
  vector<double> initialState;  //!< Initial state vector. 
  vector<double> controlsLowerBounds;//!< Lower bounds of the control space, assumed box constraint.
  vector<double> controlsUpperBounds;//!< Upper bounds of the control space, assumed box constraint.
  int noiseModel;               //!< Choice of the noise model. 
  int forwardModel;             //!< Choice of the forward model. 
  int nStochOptPerState;        //!< Number of repeated stochastic optimization runs at every evaluation of J_k(x_k).
  bool simulateTrajectories;    //!< Flag to simulate trajectories.
  int nTrajectories;            //!< Number of trajectory simulations.
  
  /* ADP value function approximation. */
  int featuresChoice;           //!< 1-total order polynomial, 2-Gaussian KL terms.
  int pOrder;                   //!< Total order polynomial order (featuresChoice=1).
  int coefsConstructionMethod;  //!< 1-regression.
  int nRegressionSamples;       //!< Number of regression sample points (coefsConstructionMethod=1). 
  int nFeatures;                //!< Total number of features. Computed.

  /* Stochastic Optimization. */
  int optMethod;                //!< 1-SPSA, 2-ELRS, 3-NMNS, 4-SAA_NCG, 5-SAA_BFGS, 6-SARM.
  int maxOptIters;              //!< Maximum number of optimization iterations.
  int relXNormTerminateNormChoice;//!< Norm choice for computing relative position change.
  double relXNormTerminateTol;  //!< Termination tolerance of position norm change.
  int nConsecRelXNormTerminateTol;//!< Number of consecutive iterations hitting relXNormTerminateTol before terminating.
  int normChoice;               //!< Norm choice.
  int nObjMC;                   //!< MC sample size of the objective function.
  vector<double> userXInitial;  //!< User input of initial position in search space.
  bool randomizeXInitial;       //!< Randomize initial position in search space (uniform).
  int nFinalObjHighQualityMC;   //!< Monte Carlo size of high quality evaluation at the final optimized position. 
  bool displayOptProgress;      //!< Flag to display progress for stochastic optimization.
  bool displayOptSummary;       //!< Flag to display final summary for stochastic optimization.
  
  // /* SARM. */
  // bool checkInitialGradient;//!< Flag for checking initial analytic gradient against FD.
  // int SARMGainSeq;          //!< Gain sequence choice for SARM (1-1/k; 2-SPSA).
  // double gainMultiplier;    //!< Multiplier to the 1/k gain sequence. 
  // bool detectSPSAParams;    //!< Flag to automatically detect "optimal" SPSA parameters.
  // /*! \brief Input scaling factor in gain sequence a_k.

  //   Recommended such that a/(A+1)^alpha times |\hat{g}_0(x_0)| is
  //   approx equal to the smallest of the desired change magnitudes in
  //   the search space. 
  // */
  // double SPSAInputa;
  // /*! \brief Input scaling factor in gain sequence c_k.

  //   Recommended to be approximately the magnitude of noise standard
  //   deviation.
  // */
  // double SPSAInputc;
  // /*! \brief Input denominator term (stability constant) in gain
  //   sequence a_k.

  //   Recommended to be proportional (e.g., 10% of) to the max number of
  //   iterations allowed or expected. 
  // */
  // double SPSAInputA;
  // /*! \brief Input denominator exponent in gain sequence a_k. 
    
  //   Recommended to be 0.602. 
  // */
  // double SPSAInputalpha;
  // /*! \brief Input denominator exponent in gain sequence c_k. 

  //   Recommended to be 0.101. 
  // */
  // double SPSAInputgamma;

  /*! \fn InputParams();
    
    \brief Default constructor of the InputParams class.
  */
  InputParams();

  /*! \fn InputParams(string);
    
    \brief Constructor of the InputParams class.
    
    \param fName Input file name.
  */
  InputParams(string const);

  /*! \fn InputParams(InputParams const &);
    
    \brief Copy constructor of the InputParams class.
    
    \param other Reference to the source class to be copied from.
  */
  InputParams(InputParams const &);

  /*! \fn ~InputParams();
    
    \brief Destructor of the InputParams class.

  */
  ~InputParams();

  /*! \fn InputParams& operator=(InputParams const &);
    
    \brief The = operator of the InputParams class.

    \param rhs Reference to input class.
    
    \return Reference to output class.
  */
  InputParams& operator=(InputParams const &);

  /*! \fn void setDefaultValues();
    
    \brief Sets the default parameter values.
  */
  void setDefaultValues();

  /*! \fn void readInputParamsFile(string const);
    
    \brief Reads the input parameters from input file.

    \param fName Input file name.
  */
  void readInputParamsFile(string const);

  /*! \fn void addMPIInfo(int const, int const);
    
    \brief Stores the MPI variables in the Controls structure.
    
    \param nTasks Total number of CPUs.
    \param rank Current CPU rank.
  */
  void addMPIInfo(int const, int const);
  
};

#endif
