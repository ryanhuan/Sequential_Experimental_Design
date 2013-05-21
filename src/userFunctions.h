/*! \file userFunctions.h 
  
  \brief A collection of various functions needed for the algorithm,
defined by the user.
*/
#ifndef _USERFUNCTIONS_H
#define _USERFUNCTIONS_H

#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <vector>

#include "inputParams.h"
#include "stochasticSearch.h"
#include "structDef.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

using namespace std;

/*! \fn void setNoiseStdDev(InputParams const &, vector<double> &);
  
  \brief Computes the Gaussian noise standard deviation.
  
  \param algParams Reference to algorithm parameters.
  \param noiseStdDev Reference to noise standard deviation.
*/
void setNoiseStdDev(InputParams const &, vector<double> &);

/*! \fn void forwardModel(InputParams const &, vector<double> const &,
  vector<double> const &, vector<double> const &, vector<double> &);
  
  \brief Evaluates the forward model.
  
  \param algParams Reference to algorithm parameters.
  \param theta Reference to model parameters.
  \param state Reference to current state.
  \param control Reference to current control.
  \param modelOutputs Reference to forward model outputs.
*/
void forwardModel(InputParams const &, vector<double> const &, 
		  vector<double> const &, vector<double> const &, 
		  vector<double> &);

/*! \fn void generateDisturbance(InputParams const &, vector<double>
  const &, vector<double> const &, gsl_rng*, vector<double> &);
			 
  \brief Generates a single instance of the disturbance. (overloaded)
  
  \param algParams Reference to algorithm parameters.
  \param state Reference to current state.
  \param control Reference to current control.
  \param generator Random number generator.
  \param disturbance Reference to disturbance.
*/
void generateDisturbance(InputParams const &, vector<double> const &, 
			 vector<double> const &, gsl_rng*, vector<double> &);

/*! \fn void generateDisturbance(InputParams const &, vector<double>
  const &, vector<double> const &, vector<double> const &, gsl_rng*,
  vector<double> &);
			 
  \brief Generates a single instance of the disturbance. (overloaded)
  
  \param algParams Reference to algorithm parameters.
  \param theta Reference to true parameter.
  \param state Reference to current state.
  \param control Reference to current control.
  \param generator Random number generator.
  \param disturbance Reference to disturbance.
*/
void generateDisturbance(InputParams const &, vector<double> const &, 
			 vector<double> const &, vector<double> const &, 
			 gsl_rng*, vector<double> &);

/*! \fn void systemEquation(InputParams const &, vector<double> const
  &, vector<double> const &, vector<double> const &, vector<double>
  &);
  
  \brief Evaluates the system equation to obtain the next state. 
  
  \param algParams Reference to algorithm parameters.
  \param state Reference to current state.
  \param control Reference to current control.
  \param disturbance Reference to current disturbance.
  \param newState Reference to new state.
*/
void systemEquation(InputParams const &, vector<double> const &, 
		    vector<double> const &, vector<double> const &,
		    vector<double> &);

/*! \fn double stageCost(InputParams const &, vector<double> const &,
  vector<double> const &, vector<double> const &);
  
  \brief Evaluates the current stage cost.
  
  \param algParams Reference to algorithm parameters.
  \param state Reference to current state.
  \param control Reference to current control.
  \param disturbance Reference to current disturbance.

  \return Current stage cost. 
*/
double stageCost(InputParams const &, vector<double> const &, 
		 vector<double> const &, vector<double> const &);

/*! \fn double maxExpectation(InputParams const &, GenericInputs const
  &, vector<double> const &);
  
  \brief Evaluates the current value function (overloaded).
  
  \param algParams Reference to algorithm parameters.
  \param allInputs Reference to generic input.
  \param state Reference to current state.

  \return Current value function value. 
*/
double maxExpectation(InputParams const &, GenericInputs const &, 
		      vector<double> const &);

/*! \fn double maxExpectation(InputParams const &, GenericInputs const
  &, vector<double> const &);
  
  \brief Evaluates the current value function (overloaded).
  
  \param algParams Reference to algorithm parameters.
  \param allInputs Reference to generic input.
  \param state Reference to current state.
  \param finalPositionExternal Reference to output final optimization 
  position. 

  \return Current value function value. 
*/
double maxExpectation(InputParams const &, GenericInputs const &, 
		      vector<double> const &, vector<double> &);

/*! \fn double evalTerminalReward(InputParams const &, vector<double>
  const &);
  
  \brief Evaluates the terminal reward function.
  
  \param algParams Reference to algorithm parameters.
  \param state Reference to state.
  
  \return Value of terminal reward.
*/
double evalTerminalReward(InputParams const &, vector<double> const &);

#endif
