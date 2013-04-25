/*! \file linearArch.h 
  
  \brief Class and functions concerning the linear architecture of
  approximation.

*/
#ifndef _LINEARARCH_H
#define _LINEARARCH_H

#include "inputParams.h"
//#include "tools.h"

/*! \class LinearArch

  \brief Class for the linear architecture of approximation.

  This class contains elements to construct the linear architecture as
  an approximation function (i.e., a linear expansion of basis
  function or features).
*/
class LinearArch
{
  
  /* Controls. */
  InputParams algParams;            //!< Algorithm parameters.
  gsl_rng_type const * rngType;     //!< GSL random number generator type.
  gsl_rng* generator;               //!< GSL random number generator 1 for noise.
  
  /* Linear architecture. */
  /*! \brief Pointer to the original true function to be approximated.
    
    The original true function must take the argument of a double
    pointer of the location of evaluation, and returns the function
    value.
  */
  double (*trueFcn)(InputParams const &, double const * const);
  vector<double> featuresCoefs;     //!< Architecture coefficients.
  vector<double> featuresVals;      //!< Evaluated feature values storages.
  vector<double> stateSample;       //!< Input variable array.
  vector< vector<int> > refTable;   //!< Refernce table for total order polynomial.
  // double** ATrans;                  //!< Working transpose of the LHS matrix for regression, altered upon exit from LAPACK.
  // double** LHS;                     //!< Reference LHS matrix for regression.
  // double* B;                        //!< Working RHS vector for regression, altered upon exit from LAPACK.
  // double* RHS;                      //!< Reference RHS vector for regression.
  // double* soln;                     //!< Solution vector for regression.
  // double* singularValues;           //!< Storage for singular values for regression.
  // int lwork;                        //!< Size for work array for regression.
  // double* work;                     //!< Work array for regression.

public:

  /*! \fn LinearArch(InputParams const &);
    
    \brief Constructor of the LinearArch class.
    
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

  /*! \fn void permPolyOrders(int, int const, int &, int* const);
    
    \brief Computes all the permutation of indexes for the i vector
    according to the upper limit maximum order. It is used for
    computing all the terms (eg, x_1^2, x_1x_2, etc) in a
    polynomial. This is a recursive function.
    
    \param zetaRemain Zeta dimension tracker in the recursion.
    \param upperLimit The upper limit of the L1-norm of a coordinate.
    This changes at each recursion stage.
    \param ctrCoords A reference to how many coordinates have been
    produced so far in the recursion.
    \param tempCoords A 1D array to temporarily store an instance of
    the coordinate.
  */
  void permPolyOrders(int, int const, int &, int* const);

  /*! \fn void makeCoefs(double (*)(Controls const &, GenericInputs &,
    double const * const), GenericInputs &);

    \brief Computes the linear architecture coefficients.

    \param trueFcnRef Pointer to the true function. 
    \param trueFcnInputsRef Reference to the true function input struct.
  */
  void makeCoefs(double (*)(Controls const &, GenericInputs &, 
			    double const * const), GenericInputs &);

  /*! \fn void evalAllFeatures(double const * const, double * const);
    
    \brief Evaluates the feature functions at the given input variable
    value.

    \param inputVar Input variable value.
    \param storage Storage array for evaluated feature values. 
  */
  void evalAllFeatures(double const * const, double * const);

  /*! \fn double evalAllFeatures(double const * const, double * const);
    
    \brief Evaluates the feature functions at the given input variable
    value.

    \param primaryExternal Reference to primary controls. This should 
    not be needed since class already has controls information.
    \param inputVar Input variable value. 

    \return Value of architecture.
  */
  double evalArchitecture(Controls const &, double const * const);

  /*! \fn void exportCoefs(double* const);
    
    \brief Export the coefs to external storage.
    
    \param coefsExternal External coefs storage.
  */
  void exportCoefs(double* const);

};

#endif
