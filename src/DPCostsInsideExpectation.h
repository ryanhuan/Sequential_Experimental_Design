/*! \file DPCostsInsideExpectation.h 

  \brief Class and functions concerning the evaluation of the costs
  inside the expectation in the dynamic programming formulation.  
*/
#ifndef _DPCOSTSINSIDEEXPECTATION_H
#define _DPCOSTSINSIDEEXPECTATION_H

#include <iostream>
#include <math.h>
#include <vector>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "inputParams.h"

using namespace std;

/*! \class DPCostsInsideExpectation

  \brief Class for evaluating the costs inside the expectation in the
  dynamic programming formulation.
*/
class DPCostsInsideExpectation
{
  
  /* Parameters. */
  InputParams algParams;        //!< Algorithm parameters.

  /* Dynamic programming. */
  double* stateKp1;         //!< Pointer to the state at stage k+1 (new memory needed).
  double* noiseK;           //!< Pointer to the noise at stage k (new memory needed).
  double* thetaSample;      //!< Underlying parameter sample.
  
  /* PC expansion. */
  double* forwardModelOutputs;//!< Storage for forward model outputs.
  double* noiseStdDev;      //!< Noise standard deviations.
  double* z;                //!< Storage for standard normal samples that are used to create the noise data points.
  double* zetas;            //!< Storage for zeta vector (argument to the Jkp1 PC expansion).
  double** psi;             //!< Reference basis function (psi) values for a particular realization of PCE inputs.
  double** dpsi;            //!< Reference basis function (psi) derivative values for a particular realization of PCE inputs.
  double psiProd;           //!< Product of psi functions for different experiments.
  double kC;                //!< Kahan summation variable for PCE evaluation.
  double kY;                //!< Kahan summation variable for PCE evaluation.
  double kT;                //!< Kahan summation variable for PCE evaluation.
  double* dValuedStateKp1Storage; //!< Gradient of J at k+1 with respect to x at k+1.
  double** dStateKp1dControl; //!< Jacobian of x at k+1 with respect to u at k.
  
public:

  /*! \fn DPCostsInsideExpectation(Controls const &);
    
    \brief Constructor for DPCostsInsideExpectation class.
    
    \param refControls Reference to primary controls.
  */
  DPCostsInsideExpectation(Controls const &);

  /*! \fn ~DPCostsInsideExpectation();
    
    \brief Destructor for the DPCostsInsideExpectation class.
  */
  ~DPCostsInsideExpectation();

  /*! \fn void setNoiseStdDev();
    
    \brief Sets the error standard deviations.
    
    The noise are assumed to have zero-mean Gaussian noise from the
    forward model.
  */
  void setNoiseStdDev();

  /*! \fn double sampleOnce(double const * const, int const * const *
    const, double const * const, double const * const, gsl_rng*);
    
    \brief Computes a one-sample MC estimate to the costs inside the
    expectation.
    
    \param coefsJkp1Ref Coefs of the Jkp1 function.
    \param refTableJkp1Ref Reference table of the Jkp1 function.
    \param stateKRef State of stage k (fixed for each optimization 
    process).
    \param controlKRef Control of stage k (varies throughout the 
    optimization process).
    \param generatorAll Pointer to all-purpose random number 
    generator.

    \return A one-sample MC estimate to the costs inside the
    expectation.
  */
  double sampleOnce(double const * const, int const * const * const, 
		    double const * const, double const * const, 
		    gsl_rng*);

  /*! \fn void evalSystemEquation(double const * const, double const *
    const, double const * const);
    
    \brief Evaluates the system equation to obtain the state at the
    next stage k+1 give the state, control, and noise realization of
    this stage k.
    
    \param stateKRef State of stage k (fixed for each optimization 
    process).
    \param controlKRef Control of stage k (varies throughout the 
    optimization process).
    \param noiseKRef Noise generated at stage k.
  */
  void evalSystemEquation(double const * const, double const * const, 
			  double const * const);
  
  /*! \fn double evalJkp1(double const * const, int const * const *
    const, double const * const);
    
    \brief Evaluates the PCE for Jkp1 function.
    
    \param coefsJkp1Ref Coefs of the Jkp1 function.
    \param refTableJkp1Ref Reference table of the Jkp1 function.
    \param stateKp1Ref State of stage k+1.
    
    \return The value of Jkp1. 
  */
  double evalJkp1(double const * const, int const * const * const, 
		  double const * const);

  /*! \fn void sampleOnceGradient(double const * const, int const *
    const * const, double const * const, double const * const, double
    * const);
    
    \brief Computes the gradient estimator of the expected immediate
    and future costs, using analytical formula.
    
    \param coefsJkp1Ref Coefs of the Jkp1 function.
    \param refTableJkp1Ref Reference table of the Jkp1 function.
    \param stateKRef State of stage k.
    \param controlKRef Control of stage k (to be perturbed).
    \param grad Storage for computed gradient.
  */
  void sampleOnceGradient(double const * const, int const * const * const, 
			  double const * const, double const * const, 
			  double * const);
  
  /*! \fn void dValuedStateKp1(double const * const, int const * const
    * const);
    
    \brief Computes the gradient of the value function J at k+1 with
    respect to the state variable x at k+1.
    
    \param coefsJkp1Ref Coefs of the Jkp1 function.
    \param refTableJkp1Ref Reference table of the Jkp1 function.
  */
  void dValuedStateKp1(double const * const, int const * const * const);

  /*! \fn void dSystemEquationdControl(double const * const, double
    const * const, double const * const);
    
    \brief Computes the Jacobian of the system equation (i.e., the
    state variable x at k+1) with respect to the control variable.
    
    \param stateKRef State of stage k (fixed for each optimization 
    process).
    \param controlKRef Control of stage k (varies throughout the 
    optimization process).
    \param noiseKRef Noise generated at stage k.
  */
  void dSystemEquationdControl(double const * const, double const * const, 
			       double const * const);

  /*! \fn void sampleOnceGradientFD(double const * const, int const *
    const * const, double const * const, double const * const,
    gsl_rng*, double * const);
    
    \brief Computes the gradient estimator of the expected immediate
    and future costs, using finite difference with the generator
    provided.
    
    \param coefsJkp1Ref Coefs of the Jkp1 function.
    \param refTableJkp1Ref Reference table of the Jkp1 function.
    \param stateKp1Ref State of stage k+1.
    \param controlKRef Control of stage k (to be perturbed).
    \param generator Pointer to random number generator with the 
    appropriate seed.
    \param grad Storage for computed gradient.
  */
  void sampleOnceGradientFD(double const * const, int const * const * const, 
			    double const * const, double const * const, 
			    gsl_rng*, double * const);
  
};

#endif
