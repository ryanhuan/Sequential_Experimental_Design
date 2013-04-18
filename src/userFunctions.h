/*! \file userFunctions.h 
  
  \brief A collection of various functions that can be evaluated for
  different purposes.  
*/
#ifndef _USERFUNCTIONS_H
#define _USERFUNCTIONS_H

#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>

//#include "stochasticSearch.h"
#include "structDef.h"
#include "tools.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

using namespace std;

/*! \fn double evalTerminalReward(Controls const &, double const *
  const);
  
  \brief Evaluates the terminal reward function.
  
  \param primary Reference to primary controls.
  \param state Input final state.
  
  \return Value of terminal reward.
*/
double evalTerminalReward(Controls const &, double const * const);

/*! \fn double evalJk(Controls const &, double const, int const,
  double const * const, int const * const * const)
  
  \brief Evaluates or approximates the J_k function using stochastic
  optimization at the specified abscissa (transformed to state).
  
  \param primary Reference to primary controls.
  \param abscissas Array of abscissas.
  \param currentStage The current DP stage.
  \param coefsJkp1Ref Coefs of the Jkp1 function.
  \param refTableJkp1Ref Reference table of the Jkp1 function.
  
  \return Value of J_k.
*/
double evalJk(Controls const &, double const * const, int const, 
	      double const * const, int const * const * const);

#endif
