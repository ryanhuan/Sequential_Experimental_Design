/*! \file stochasticSearch.h 

  \brief Wrapper class and functions concerning the stochastic search
  algorithms.
  
*/
#ifndef _STOCHASTICSEARCH_H
#define _STOCHASTICSEARCH_H

#include <iomanip>
#include <iostream>
#include <vector>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "inputParams.h"
//#include "SARM.h"
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
  
  /* Parameters. */
  InputParams algParams;        //!< Algorithm parameters.
  bool initialized;             //!< Flag to indicator whether the class has been initialized. 
  bool positionInitialized;     //!< Flag to indicator whether the class has been position initialized. 
  
  /* Optimization. */
  int nIter;                    //!< Current iteration number.
  int storageCtr;               //!< Usually equal to nIter, but if only selected iterations are stored, it can be different. 
  vector<double> XInitial;      //!< Initial search position.
  string dataLabels;            //!< Labels for data display.
  string terminationLabel;      //!< Label for termination condition.
//  DPCostsInsideExpectation* DPCostsForHighQualityMC;//!< Class for DP costs inside the expectation (for high quality evaluation).
  double finalObjHighQualityMC; //!< Objective value on the final position of optimization from a high quality MC.
  gsl_rng_type const * rngType; //!< GSL random number generator type.
  gsl_rng* generatorNoiseForHighQualityMC; //!< GSL random number generator 1 for noise.

  /* History storage. */
  vector<int> iterHistory;      //!< History of iterations (may skip some iterations in storage).
  vector< vector<double> > XHistory;//!< History of position.
  vector<double> valHistory;    //!< History of objective value.
  vector<double> gradNormHistory;//!< History of gradient norm. 
  
  /* Algorithm-specific classes. */
//  SARM robbinsMonro;            //!< SARM class (optMeth = 6).

  /* Termination conditions. */
  int terminationCondition;     //!< Termination condition indicator. Details see stochasticSearch.cpp. 
  int nConsecRelXNorm;          //!< Counter for consecutive relative position norm change satisfying tolerance.
  double previousXNorm;         //!< Previous position norm.
  double curXNorm;              //!< Current position norm.
  bool isLocalMin;              //!< Indicator for reaching a local minimum (for determinisitic algorithms such as SAA_NCG).
  
public:

  /*! \fn StochasticSearch();
    
    \brief Default constructor for StochasticSearch class.
  */
  StochasticSearch();

  /*! \fn StochasticSearch(InputParams const &);
    
    \brief Constructor for StochasticSearch class.
    
    \param algParams Reference to algorithm parameters.
  */
  StochasticSearch(InputParams const &);

  /*! \fn ~StochasticSearch();
    
    \brief Destructor for the StochasticSearch class.
  */
  ~StochasticSearch();

  /*! \fn void initialize(InputParams const &);
    
    \brief Initializer for StochasticSearch class.
    
    \param algParams Reference to algorithm parameters.
  */
  void initialize(InputParams const &);

  /*! \fn void makeLabels();
    
    \brief Makes the data labels for progress and text file.
  */
  void makeLabels();
  
  /*! \fn void initializePosition(vector<double> const &,
    vector<double> const &);
    
    \brief General initialization of the search.

    \param XInitialRef Reference to initial search position.
  */
  void initializePosition(vector<double> const &);

  /*! \fn void search(vector<double> const &, GeneritcInputs);
    
    \brief Wrapper function to perform the desired search algorithm.

    \param state Reference to state (fixed for each optimization process).
    \param allInputs Generic inputs. 
  */
  void search(vector<double> const &, GenericInputs);

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

  /*! \fn void exportFinalPosition(vector<double> &);
    
    \brief Exports the final position from the stochastic optimization
    process.
    
    \param finalXExternal Reference to external storage for the final position.
  */
  void exportFinalPosition(vector<double> &);

  /*! \fn void makeTerminationLabel();
    
    \brief Makes the label for the termination condition.
  */
  void makeTerminationLabel();

};

#endif
