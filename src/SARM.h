/*! \file SARM.h 

  \brief Class and functions concerning the SAA with BFGS algorithm.  
*/
#ifndef _SARM_H
#define _SARM_H

#include <iomanip>
#include <iostream>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "DPCostsInsideExpectation.h"
#include "structDef.h"

using namespace std;

/*! \class SARM

  \brief Class for the Stochastic Approximation Robins-Monro 
  algorithm.
*/
class SARM
{
  
  /* Controls. */
  Controls primary;         //!< Controls.
  
  /* History storage. */
  /*! \brief Pointer to iteration counter in StochasticSearch class
    (or some external source).
    
    The intention is for the external source to be more easily
    accessible to the user, and thus not stored in this class. Hence,
    no new memory is allocated or deleted for this pointer in this
    class.
  */
  int* nIter;
  int* storageCtr;          //!< Pointer to storage counter.
  int* iterHistory;         //!< Pointer to iteration history.
  double** XHistory;        //!< Pointer to position history.
  double* valHistory;       //!< Pointer to objective value history.
  double* gradNormHistory;  //!< Pointer to gradient norm history.
  bool flagCheckGradientFD; //!< Flag for internal checking gradient against FD at every computation.
  DPCostsInsideExpectation* DPCosts;//!< Class for DP costs inside the expectation (single sample).
  
  /* Optimization. */
  double* XCur;             //!< Current (unnormalized) position.
  double valTemp;           //!< Temporary objective value.
  double valCur;            //!< Current objective value.
  double* gradTemp;         //!< Temporary (unnormalized) gradient storage.
  double* gradTemp2;        //!< Temporary (unnormalized) gradient storage 2.
  double* gradCur;          //!< Current (unnormalized) gradient.
  gsl_rng_type const * rngType; //!< GSL random number generator type.
  gsl_rng* generatorNoise1; //!< GSL random number generator 1 for noise.
  gsl_rng* generatorNoise2; //!< GSL random number generator 2 for noise.
  gsl_rng* generatorNoiseFDCheck;//!< GSL random number generator for noise, for gradient checking using FD.
  int CRNSeed;              //!< CRN seed storage.
  
  /* Parameters. */
  double a_k;               //!< Gain sequence.
  /*! \brief Scaling factor in gain sequence a_k.
    
    Recommended such that a/(A+1)^alpha times |\hat{g}_0(x_0)| is
    approx equal to the smallest of the desired change magnitudes in
    the search space. 
  */
  double SPSAa;
  /*! \brief Denominator term (stability constant) in gain sequence a_k.
    
    Recommended to be proportional (e.g., 10% of) to the max number of
    iterations allowed or expected. 
  */
  double SPSAA;
  /*! \brief Denominator exponent in gain sequence a_k. 
    
    Recommended to be 0.602. 
  */
  double SPSAalpha;
  
public:
  
  /*! \fn SARM(Controls const &, int* const, int* const, int* const,
    double** const, double* const, double* const);
    
    \brief Constructor for SARM class.
    
    \param refControls Reference to primary controls.
    \param nIterRef Pointer to the nIter element in StochasticSearch
    class (or some external source).
    \param storageCtrRef Pointer to the storageCtr element.
    \param iterHistoryRef Pointer to the iterHistory element.
    \param XHistoryRef Pointer to the XHistory element.
    \param valHistoryRef Pointer to the valHistory element.
    \param gradNormHistoryRef Pointer to the gradNormHistory element.
  */
  SARM(Controls const &, int* const, int* const, int* const, double** const, 
       double* const, double* const);
  
  /*! \fn ~SARM();
    
    \brief Destructor for the SARM class. Note that pointers to
    elements in StochasticSearch class (or some external source) must
    not be freed here.
  */
  ~SARM();
  
  /*! \fn void initialPosition(double const * const, int const * const
    * const, double const * const);
    
    \brief Performs initializations that are specific to this search
    algorithm.

    \param coefsJkp1Ref Coefs of the Jkp1 function.
    \param refTableJkp1Ref Reference table of the Jkp1 function.
    \param stateKRef State of stage k (fixed for each optimization process).
  */
  void initialPosition(double const * const, int const * const * const, 
		       double const * const);

  /*! \fn void iterateOnce(double const * const, int const * const *
    const, double const * const);
    
    \brief Takes one step/iteration for the SARM algorithm. This is
    the bulk of the SARM algorithm.

    \param coefsJkp1Ref Coefs of the Jkp1 function.
    \param refTableJkp1Ref Reference table of the Jkp1 function.
    \param stateKRef State of stage k (fixed for each optimization process).
  */
  void iterateOnce(double const * const, int const * const * const, 
		   double const * const);
  
  /*! \fn void update_a_k();
    
    \brief Updates the gain sequence at the current iteration.
  */
  void update_a_k();

};

#endif
