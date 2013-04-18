/*! \file stochasticSearch.h 

  \brief Wrapper class and functions concerning the stochastic search
  algorithms.
*/
#ifndef _STOCHASTICSEARCH_H
#define _STOCHASTICSEARCH_H

#include <iomanip>
#include <iostream>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "SARM.h"
#include "structDef.h"
#include "tools.h"

using namespace std;

/*! \class StochasticSearch

  \brief Class for the stochastic search wrapper.

  This wrapper class includes some common elements that are needed for
  all stochastic search algorithms, then further includes classes that
  are specific to the different search algorithms.
*/
class StochasticSearch
{
  
  /* Controls. */
  Controls primary;             //!< Controls.
  
  /* Optimization. */
  /*! \brief Current iteration number (1 element array). 

    We want to use a pointer instead of a simple int so that the child
    classes within StochasticSearch can directly point to the same
    memory storage.
   */
  int* nIter;
  int* storageCtr;              //!< Usually equal to nIter, but if only selected iterations are stored, it can be different. 
  string dataLabels;            //!< Labels for data display.
  string terminationLabel;      //!< Label for termination condition.
  DPCostsInsideExpectation* DPCostsForHighQualityMC;//!< Class for DP costs inside the expectation (for high quality evaluation).
  double finalObjHighQualityMC; //!< Objective value on the final position of optimization from a high quality MC.
  gsl_rng_type const * rngType; //!< GSL random number generator type.
  gsl_rng* generatorNoiseForHighQualityMC; //!< GSL random number generator 1 for noise.

  /* History storage. */
  int* iterHistory;             //!< History of iterations (may skip some iterations in storage).
  double** XHistory;            //!< History of position.
  double* valHistory;           //!< History of objective value.
  double* gradNormHistory;      //!< History of gradient norm. 
  
  /* Algorithm-specific classes. */
  SARM* robbinsMonro;           //!< Pointer to SARM class (optMeth = 6).

  /* Termination conditions. */
  int terminationCondition;     //!< Termination condition indicator. Details see stochasticSearch.cpp. 
  int nConsecRelXNorm;          //!< Counter for consecutive relative position norm change satisfying tolerance.
  double previousXNorm;         //!< Previous position norm.
  double curXNorm;              //!< Current position norm.
  bool isLocalMin;              //!< Indicator for reaching a local minimum (for determinisitic algorithms such as SAA_NCG).
  
public:

  /*! \fn StochasticSearch(Controls const &);
    
    \brief Constructor for StochasticSearch class.
    
    \param refControls Reference to primary controls.
  */
  StochasticSearch(Controls const &);

  /*! \fn ~StochasticSearch();
    
    \brief Destructor for the StochasticSearch class.
  */
  ~StochasticSearch();

  /*! \fn void makeLabels();
    
    \brief Makes the data labels for progress and text file.
  */
  void makeLabels();
  
  /*! \fn void initialPosition(double (*)(Controls const &, double
    const * const), double const * const);
    
    \brief General initialization of the search.

    \param innerFcn Pointer to the inner function. 
    \param state State (fixed for each optimization process).
  */
  void initialPosition(double (*)(Controls const &, double const * const), 
		       double const * const);

  /*! \fn void search(double (*)(Controls const &, double const *
    const), double const * const);
    
    \brief Wrapper function to perform the desired search algorithm.

    \param innerFcn Pointer to the inner function. 
    \param state State (fixed for each optimization process).
  */
  void search(double (*)(Controls const &, double const * const), 
	      double const * const);

  /*! \fn int checkTermination();
    
    \brief Checks whether any of the applicable termination conditions
    are true.

    \return 0 - condition unsatisfied, 1 - iterations, 2 -
    algorithm-relevant function evaluations, 3 - total function
    evaluations, 4 - relative position norm, 5 - relative objective
    value, 6 - gradient norm, 7 - simplex diameter.
  */
  int checkTermination();
  
  /*! \fn void displaySummary();
    
    \brief Displays the summary stats.
  */
  void displaySummary();

  /*! \fn double exportFinalObjHighQualityMC();
    
    \brief Exports the final objective function estimate from the
    stochastic optimization process, based on a high-quality Monte
    Carlo estimate.

    \return Final high-quality objective function estimate.
  */
  double exportFinalObjHighQualityMC();

  /*! \fn void exportFinalPosition(double * const);
    
    \brief Exports the final position from the stochastic optimization
    process.
    
    \param finalXExternal External storage for the final position.
  */
  void exportFinalPosition(double * const);

  /*! \fn void makeTerminationLabel();
    
    \brief Makes the label for the termination condition.
  */
  void makeTerminationLabel();

};

#endif
