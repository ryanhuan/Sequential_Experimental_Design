/*! \file trajectories.h 

  \brief Class and functions concerning the forward trajectories
  (simulations, forward paths) of the experiments.

*/
#ifndef _TRAJECTORIES_H
#define _TRAJECTORIES_H

#include <iostream>
#include <vector>

#include "structDef.h"
#include "userFunctions.h"

using namespace std;

/*! \class Trajectories

  \brief Class for forward trajectories. We have elected to use
  pointers rather than vectors in this class to facilitate
  MPI. IMPORTANT: this class itself cannot be passed via MPI. Copy
  constructor involves a deep copy.

*/
class Trajectories
{
  
  /* General. */
  InputParams algParams;                //!< Algorithm parameters.
  gsl_rng_type const * rngType;         //!< GSL random number generator type.
  gsl_rng* generator;                   //!< GSL random number generator 1 for noise.
  bool initialized;                     //!< Flag to indicator whether the class has been initialized. 
  
  /* Master list of trajectories. */
  vector<int> nLocalTrajectoriesAll;    //!< All local numbers of trajectories. 
  vector<int> nLocalTrajectoriesAllSum; //!< Cumulative sums of all local numbers of trajectories. 
  double** thetaSamplesAll;             //!< All samples of "true" parameter values.
  double** statesAll;                   //!< All states over all stages of a trajectories, including x_0.
  double** controlsAll;                 //!< All controls over all stages of a trajectories.
  double** disturbanceAll;              //!< All disturbance over all stages of a trajectories.
  double* finalRewardsAll;              //!< All final rewards of the trajectories.

  /* Local list of trajectories. */
  int nLocalTrajectories;               //!< Local number of trajectories.
  double** thetaSamples;                //!< Local samples of "true" parameter values.
  double** states;                      //!< Local states over all stages of a trajectories, including x_0.
  double** controls;                    //!< Local controls over all stages of a trajectories.
  double** disturbance;                 //!< Local disturbance over all stages of a trajectories.
  double* finalRewards;                 //!< Local final rewards of the trajectories.
  vector<double> tempState;             //!< Temporary state vector. 
  vector<double> tempControl;           //!< Temporary control vector. 
  vector<double> tempDisturbance;       //!< Temporary disturbance vector. 
  vector<double> tempNewState;          //!< Temporary new state vector. 

public:
  
  /*! \fn Trajectories();
    
    \brief Default constructor of the Trajectories class. This should
    only be used if it is followed up with initialize.
  */
  Trajectories();

  /*! \fn Trajectories(InputParams const &);
    
    \brief Constructor of the Trajectories class taking in argument.
    
    \param algParams Reference to algorithm parameters.
  */
  Trajectories(InputParams const &);

  /*! \fn Trajectories(Trajectories const &);
    
    \brief Copy constructor of the Trajectories class. Deep copy is
    performed, with new memory allocations.
    
    \param other Reference to the source class to be copied from.
  */
  Trajectories(Trajectories const &);

  /*! \fn ~Trajectories();
    
    \brief Destructor of the Trajectories class.
  */
  ~Trajectories();

  /*! \fn Trajectories& operator=(Trajectories const &);
    
    \brief The = operator of the Trajectories class. Deep copy is
    performed, memory assumed to be already allocated (no new
    allocation).

    \param rhs Reference to input class.
    
    \return Reference to output class.
  */
  Trajectories& operator=(Trajectories const &);

  /*! \fn void initialize(InputParams const &);
    
    \brief Initializations of the Trajectories class.
    
    \param algParams Reference to algorithm parameters.
  */
  void initialize(InputParams const &);

  /*! \fn void simulateTrajectories(vector<LinearArch> &arch);
    
    \brief Simulate all local trajectories and pass back to master
    storage via MPI.

    \param arch Reference to vector of linear architectures of 
    value function approximations.
  */
  void simulateTrajectories(vector<LinearArch> &arch);

  /*! \fn void writeTrajectoriesToFile(string const);
    
    \brief Write the trajectories to file.
    
    \param fName File name in string form.
  */
  void writeTrajectoriesToFile(string const);
  
};

#endif
