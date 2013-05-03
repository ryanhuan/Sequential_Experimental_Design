/*! \file linearArch.h 
  
  \brief Class and functions concerning the linear architecture of
  approximation.

*/
#ifndef _LINEARARCH_H
#define _LINEARARCH_H

#include <math.h>
#include <vector>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "inputParams.h"
#include "tools.h"
#include "structDef.h"

/*! \class LinearArch

  \brief Class for the linear architecture of approximation.

  This class contains elements to construct the linear architecture as
  an approximation function (i.e., a linear expansion of basis
  function or features).

*/
class LinearArch
{
  
  /* Parameters. */
  InputParams algParams;            //!< Algorithm parameters.
  gsl_rng_type const * rngType;     //!< GSL random number generator type.
  gsl_rng* generator;               //!< GSL random number generator 1 for noise.
  bool initialized;                 //!< Flag to indicator whether the class has been initialized. 
  
  /* Linear architecture. */
  /*! \brief Pointer to the original true function to be approximated.
  */
  double (*trueFcn)(InputParams const &, GenericInputs const &, vector<double> const &);
  GenericInputs trueFcnInputs;      //!< Inputs for the true function. 
  vector<double> featuresCoefs;     //!< Architecture coefficients.
  vector<double> featuresVals;      //!< Evaluated feature values storages.
  vector<double> stateSample;       //!< Input variable array.
  vector< vector<int> > refTable;   //!< Refernce table for total order polynomial.
  vector< vector<double> > ATrans;  //!< Working transpose of the LHS matrix for regression, altered upon exit from LAPACK.
  vector< vector<double> > LHS;     //!< Reference LHS matrix for regression.
  vector<double> B;                 //!< Working RHS vector for regression, altered upon exit from LAPACK.
  vector<double> RHS;               //!< Reference RHS vector for regression.
  vector<double> soln;              //!< Solution vector for regression.
  vector<double> singularValues;    //!< Storage for singular values for regression.
  int lwork;                        //!< Size for work array for regression.
  vector<double> work;              //!< Work array for regression.

public:

  /*! \fn LinearArch();
    
    \brief Default constructor of the LinearArch class. This should
    only be used if it is followed up with initialize.
  */
  LinearArch();

  /*! \fn LinearArch(InputParams const &);
    
    \brief Constructor of the LinearArch class taking in argument.
    
    \param algParams Reference to algorithm parameters.
  */
  LinearArch(InputParams const &);

  /*! \fn LinearArch(LinearArch const &);
    
    \brief Copy constructor of the LinearArch class.
    
    \param other Reference to the source class to be copied from.
  */
  LinearArch(LinearArch const &);

  /*! \fn ~LinearArch();
    
    \brief Destructor of the LinearArch class.
  */
  ~LinearArch();

  /*! \fn LinearArch& operator=(LinearArch const &);
    
    \brief The = operator of the LinearArch class.

    \param rhs Reference to input class.
    
    \return Reference to output class.
  */
  LinearArch& operator=(LinearArch const &);

  /*! \fn void initialize(InputParams const &);
    
    \brief Initializations of the LinearArch class.
    
    \param algParams Reference to algorithm parameters.
  */
  void initialize(InputParams const &);

  /*! \fn void permPolyOrders(int, int const, int &, vector<int> &);
    
    \brief Computes all the permutation of indexes for the i vector
    according to the upper limit maximum order. It is used for
    computing all the terms (eg, x_1^2, x_1x_2, etc) in a
    polynomial. This is a recursive function.
    
    \param zetaRemain Zeta dimension tracker in the recursion.
    \param upperLimit The upper limit of the L1-norm of a coordinate.
    This changes at each recursion stage.
    \param ctrCoords A reference to how many coordinates have been
    produced so far in the recursion.
    \param tempCoords A vector to temporarily store an instance of
    the coordinate.
  */
  void permPolyOrders(int, int const, int &, vector<int> &);

  /*! \fn void makeCoefs(double (*)(InputParams const &, GenericInputs
    const &, vector<double> const &), GenericInputs const &);

    \brief Computes the linear architecture coefficients.

    \param trueFcnRef Pointer to the true function. 
    \param trueFcnInputsRef Reference to the true function input struct.
  */
  void makeCoefs(double (*)(InputParams const &, GenericInputs const &, 
			    vector<double> const &), GenericInputs const &);

  /*! \fn void evalAllFeatures(vector<double> const &, vector<double>
    &);
    
    \brief Evaluates the feature functions at the given input variable
    vector.

    \param inputVar Input variable vector.
    \param storage Storage vector for evaluated feature values. 
  */
  void evalAllFeatures(vector<double> const &, vector<double> &);

  /*! \fn double evalArchitecture(InputParams const &, vector<double>
    const &);
    
    \brief Evaluates the feature functions at the given input variable
    value.

    \param algParamsExternal Reference to external algorithm parameters. 
    This should not be used since class already has input parameters 
    information.
    \param inputVar Input variable value. 

    \return Value of architecture.
  */
  double evalArchitecture(InputParams const &, vector<double> const &);

  /*! \fn void exportCoefs(vector<double> &);
    
    \brief Export the coefs to external storage.
    
    \param coefsExternal External coefs storage vector.
  */
  void exportCoefs(vector<double> &);

};

#endif
