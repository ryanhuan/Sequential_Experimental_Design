/*! \file trajectory.h 

  \brief Class and functions concerning the forward trajectory
  (simulation) of the experiments.

*/
#ifndef _TRAJECTORY_H
#define _TRAJECTORY_H

#include <iostream>
#include <vector>

using namespace std;

/*! \class Trajectory

  \brief Class for a forward trajectory. All information about this
  trajectory is stored in a contiguous chunk of memory. The are
  ordered as: thetaSample[], x_0[], u_0[], w_0[], x_1[], u_1[], w_1[],
  x_2[], final cost.

*/
class Trajectory
{
  
public:
  
  /* General. */
  InputParams algParams;                //!< Algorithm parameters.

  /* Trajectory. */
  vector<double> thetaSamples;          //!< Sample of "true" parameter value.
  vector< vector<double> > statesAll;   //!< States over all stages of a trajectory, including x_0.
  vector< vector<double> > controlsAll; //!< Controls over all stages of a trajectory.
  vector< vector<double> > disturbanceAll;//!< Disturbance over all stages of a trajectory.
  double finalReward;                   //!< Final reward of the trajectory.
  
};

#endif
