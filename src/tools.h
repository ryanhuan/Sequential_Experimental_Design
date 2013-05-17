/*! \file tools.h 

  \brief Various tools.

*/
#ifndef _TOOLS_H
#define _TOOLS_H

#include <iostream>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <vector>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
//#include "gtest/gtest.h"

using namespace std;

/*! \fn void randMultivariateNormal(vector< vector<double> > const &,
  gsl_rng const *, vector<double> &);
  
  \brief Generates a sample of a zero-mean multivariate normal random
  vector from a square-root of the covariance matrix.
  
  \param sqrtCov Reference to the square-root of the covariance matrix.
  \param generator GSL random number generator.
  \param sample Reference to the resulting sample vector.
*/
void randMultivariateNormal(vector< vector<double> > const &, gsl_rng const *, 
			    vector<double> &);

/*! \fn double vectorNorm(vector<double> const &, int const);
  
  \brief Computes a norm of a 1D double array.
  
  \param array Reference to the double vector.
  \param normChoice Choice of norm (0-infinity norm, 1-1 norm, 2-2 norm).

  \return Norm of the vector.
*/
double vectorNorm(vector<double> const &, int const);

/*! \fn int nCr(int const, int const);

  \brief Computes the combination term, n choose r.
  
  \param n Total elements.
  \param r Number of elements for combination.
  
  \return The computed combination. If this <0, an overflow has
  occurred.
*/
int nCr(int const, int const);

/*! \fn int factorial(int const);

  \brief Computes the factorial of an integer.
  
  \param n Factorial argument.
  
  \return The computed factorial. If this <0, an overflow has
  occurred.
*/
int factorial(int const);

/*! \fn int binarySearch(double const, double const * const, int
  const, double const);
  
  \biref Performs a binary search of a query value in a sorted 1D
  array with distinct entries. The corresponding index value is
  returned.
  
  \param query Search query.
  \param refList Sorted 1D array with distinct entries.
  \param nLength Length of 1D array.
  \param dCompTol Comparison tolerance for double type.
  
  \return Returns the index of the query. If not found, returns -1.
*/
int binarySearch(double const, double const * const, int const, double const);

/*! \fn void threeTermLegendre(double const, int const, double*
  const);
  
  \brief Computes the Legendre polynomial values evaluated for some x
  via the 3-term recursive formula. This is for the weight function of
  unity between -1 to 1.

  \param x Point of evaluation.
  \param maxDegree Maximum degree of polynomial to evaluate.
  \param P Array storing function values of all degrees up to the 
  maximum.
*/
void threeTermLegendre(double const, int const, double* const);

/*! \fn void threeTermHermite(double const, int const, double*
  const);
  
  \brief Computes the probabilists' Hermite polynomial values
  evaluated for some x via the 3-term recursive formula. This is for
  the weight function of exp(-x^2/2).

  \param x Point of evaluation.
  \param maxDegree Maximum degree of polynomial to evaluate.
  \param P Array storing function values of all degrees up to the 
  maximum.
*/
void threeTermHermite(double const, int const, double* const);

/*! \fn void threeTermLegendreDerivatives(double const, int const,
  double const * const, double* const);
  
  \brief Computes the Legendre polynomial derivative values evaluated
  for some x by directly differentiating the 3-term recursive formula
  (i.e., not the 3-term recurrence differential relations, since that
  becomes singular at |x|=1). This is for the weight function of unity
  between -1 to 1.

  \param x Point of evaluation.
  \param maxDegree Maximum degree of polynomial to evaluate.
  \param P Input array of polynomial evaluations of all degrees up to
  the maximum.
  \param dP Array storing derivative results of all degrees up to the 
  maximum.
*/
void threeTermLegendreDerivatives(double const, int const, 
				  double const * const, double* const);

/*! \fn void threeTermHermiteDerivatives(double const, int const,
  double const * const, double* const);
  
  \brief Computes the probabilists' Hermite polynomial derivative
  values evaluated for some x by the 3-term differential recursive
  formula. This is for the weight function of exp(-x^2/2).

  \param x Point of evaluation.
  \param maxDegree Maximum degree of polynomial to evaluate.
  \param P Input array of polynomial evaluations of all degrees up to
  the maximum.
  \param dP Array storing derivative results of all degrees up to the 
  maximum.
*/
void threeTermHermiteDerivatives(double const, int const, 
				 double const * const, double* const);

/*! \fn void eigSymTridiag(int const &, double* const, double* const,
  double const &, double* const, double* const, int* const, double*
  const, int* const);

  \brief Computes the eigenvalues and eigenvectors of a real symmetric
  tridiagonal matrix.
 
  \param nLength Length of the square matrix.
  \param diagonal 1D array of diagonal entries.
  \param subdiagonal 1D array of subdiagonal entries.
  \param abstol Absolute eigenvalue error tolerance, obtained from 
  DLAMCH LAPACK function.
  \param eigenVals Array of resulting eigenvalues.
  \param eigenVecs 2D array of resulting eigenvectors.
  \param ifailOut Index numbers for which eigenvectors failed to converge.
  \param work Work double array.
  \param iwork Work int array.
*/
void eigSymTridiag(int const &, double* const, double* const, double const &, 
		   double* const, double* const, int* const, double* const, 
		   int* const);

/*! \fn void mvProduct(int const, int const, double* const * const,
  double* const, double* const);
  
  \brief A wrapper function that computes matrix-vector
  multiplication. DGEMV routine from BLAS is used. 

  \param nLength1 Number of rows of matrix.
  \param nLength2 Number of columns of matrix.
  \param A Matrix of interest.
  \param x Vector right-multiplying A.
  \param y Results vector.
*/
void mvProduct(int const, int const, double* const * const, 
	       double* const, double* const);

/*! \fn void choleskyDecomp(int const, double const * const * const,
  double* const * const);
  
  \brief A wrapper function that computes the Cholesky decomposition
  on a symmetric positive-definite matrix. The factored lower
  triangular matrix is produced. DPOTRF routine from LAPACK is used.

  \param nLength Length of the square matrix.
  \param A The symmetric positive-definite matrix.
  \param cholesky Storage for the Cholesky (upper triangular) matrix.
*/
void choleskyDecomp(int const, double const * const * const, 
		    double* const * const);

/*! \fn double invCondNumber(int const, double* const * const, int*
  const, double* const, int* const);
  
  \brief A wrapper function that estimates the inverse condition
  number of a matrix. DLANGE, DGETRF, and DGECON FORTRAN routines
  from LAPACK are used. 
  
  \param length Length of the matrix.
  \param A Matrix of interest.
  \param ipivOut Array for storing the temporary pivoting vector.
  \param work Work double array.
  \param iwork Work int array.
*/
double invCondNumber(int const, double* const * const, int* const, 
		     double* const, int* const);

/*! \fn void solveBandedLinearSystem(int const, int const, int const,
  double* const * const, int* const, double* const);
  
  \brief A wrapper function that solves a banded linear system. DGBSV
  routine from LAPACK is used.

  \param nLength Size of the square matrix.
  \param nSubdiagonals Number of nonzero subdiagonals.
  \param nSuperdiagonals Number of nonzero superdiagonals.
  \param A LHS matrix in standard (not LAPACK-required) format.
  \param ABIn Storage for LHS matrix in LAPACK-required format. This
  function will automatically create ABIn.
  \param ipivOut Storage for permutation vector.
  \param uIn RHS vector, which also contains the solution vector upon 
  exit.
*/
void solveBandedLinearSystem(int const, int const, int const, 
			     double const * const * const, double* const * const, 
			     int* const, double* const);

/*! \fn void linearLeastSquares(vector< vector<double> > &,
  vector<double> &, vector<double> &, vector<double> &, int const,
  vector<double> &);
  
  \brief A wrapper function that solves a banded linear system. DGBSV
  routine from LAPACK is used.

  \param ATrans Reference to the transpose of the vandermonde matrix 
  (LAPACK-required).
  \param sOut Reference to the array for storing computed singular values. 
  \param B Reference to the RHS vector. 
  \param work Reference to the work double array.
  \param workLength Work double array length.
  \param soln Reference to the solution vector. 
*/
void linearLeastSquares(vector< vector<double> > &, vector<double> &, 
			vector<double> &, vector<double> &, int const, 
			vector<double> &);

/*! \fn void quickSort(T* const, U* const, V** const, int const, int
  const);

  \brief Sorts a 1D array of elements in ascending order using a
  recursive implementation of quick sort. It can also sort a secondary
  1D array and a 2D array according to the sorting of the 1st 1D
  array.
  
  \param arrayA A 1D array to be sorted.
  \param arrayB Another 1D array to be sorted according to the sorting
  of arrayA.
  \param arrayC A 2D array to be sorted according to the sorting of
  arrayA.
  \param left Starting position to the section in arrayA to be sorted.
  \param right Ending position to the section in arrayA to be sorted.
*/
template <class T, class U, class V>
void quickSort(T* const, U* const, V** const, int const, int const);

/*! \fn int partition(T* const, U* const, V** const, int const, int
  const, int const);
  
  \brief The bulk of the quick sort. It sorts the pivot, and then
  partitions the array to be sorted into 2 subarrays from the pivot
  (which would be at the correct global position after the previous
  iteration).
  
  \param arrayA A 1D array to be sorted.
  \param arrayB Another 1D array to be sorted according to the sorting
  of arrayA.
  \param arrayC A 2D array to be sorted according to the sorting of
  arrayA.
  \param left Starting position to the section in arrayA to be sorted.
  \param right Ending position to the section in arrayA to be sorted.
  \param pivotIndex The element index selected as the pivot.
  
  \return The final (correct global) element index of the selected
  pivot.
*/
template <class T, class U, class V>
int partition(T* const, U* const, V** const, int const, int const, int const);

/*! \fn void maxOf1DArray(T const * const, int const, int&, T&);
  
  \brief Finds the maximum index and value in a 1D array.
  
  \param V The 1D array of interest.
  \param length Length of the array.
  \param maxIndex Reference to store the index of the maximum element.
  \param maxVal Reference to store the maximum value.
*/
template <class T>
void maxOf1DArray(T const * const, int const, int&, T&);

/*! \fn void minOf1DArray(T const * const, int const, T &, int &);
  
  \brief Finds the minimum value and its correponding index in an 1D
  array.
  
  \param V Pointer to the 1D array of interest. 
  \param length Length of the array.
  \param minVal Reference to the minimum value variable.
  \param minIndex Reference to the minimum index variable.
*/
template <class T>
void minOf1DArray(T const * const, int const, T &, int &);

/*! \fn void maxOfAbs2DArray(T const * const * const, int const, int
  const, int&, int&, T&);

  \brief Finds the index and the maximum absolute value in a 2D array.
  
  \param V The 2D array of interest.
  \param nLen1 Length of 1st dimension.
  \param nLen2 Length of 2nd dimension.
  \param maxIndex1 Reference to store the 1st dimension index of the
  maximum element.
  \param maxIndex2 Reference to store the 2nd dimension index of the
  maximum element.
  \param maxVal Reference to store the maximum value.
*/
template <class T>
void maxOfAbs2DArray(T const * const * const, int const, int const,
		     int&, int&, T&);

/*! \fn void swapPtr(T*&, T*&);
  
  \brief Swaps two pointers (used in sorting multiD arrays). Note that
  the original array must be allocated using a for-loop (ie, not a
  continuous chunk), or else there would be problems when freeing the
  memory after sorting.
  
  \param a Reference to the first pointer.
  \param b Reference to the second pointer.
*/
template <class T>
void swapPtr(T*&, T*&);

/*! \fn string num2string(T const);

  \brief Converts a numerical number to a string (like the inverse of
  atoi or atof).
  
  \param input Numerical number.

  \return String form of the input number.
*/
template <class T>
string num2string(T const);

/*! \fn void minOfVector(const T&, int&, U&);

  \brief Finds the minimum value and its correponding index in a
  vector.
  
  \param V Reference to the vector of interest.
  \param minIndex Reference to the minimum index variable.
  \param minVal Reference to the minimum value variable.
*/
template <class T, class U>
void minOfVector(const T&, int&, U&);

/*! \fn int sgn(T const d);
  
  \brief Sign function.
  
  \param d Input value.

  \return 1 if positive, -1 if negative, 0 if 0.
*/
template <class T>
int sgn(T const);

/*! \def EPS
  \brief EPS = 1.0e-15, some small value above machine precision.
*/
#define EPS 1.0e-15




/* Template definitions. */

template <class T, class U, class V>
void quickSort(T* const arrayA, U* const arrayB, V** const arrayC, 
	       int const left, int const right)
{
  
  /* Partitioning and sorting. */
  if (left < right){
    int pivotIndex = left;
    int pivotNewIndex = partition<T,U,V>(arrayA,arrayB,arrayC,left,right,pivotIndex);
    quickSort<T,U,V>(arrayA,arrayB,arrayC,left,pivotNewIndex-1);
    quickSort<T,U,V>(arrayA,arrayB,arrayC,pivotNewIndex+1,right);
  }

}

template <class T, class U, class V>
int partition(T* const arrayA, U* const arrayB, V** const arrayC, 
	      int const left, int const right, int const pivotIndex)
{

  /* Initializations. */
  T pivotValue(arrayA[pivotIndex]);

  swap<T>(arrayA[right],arrayA[pivotIndex]);
  if (arrayB!=NULL)
    swap<U>(arrayB[right],arrayB[pivotIndex]);
  if (arrayC!=NULL)
    swapPtr<V>(arrayC[right],arrayC[pivotIndex]);

  int storeIndex = left;
  for (int i=left; i<right; ++i) {
    if (arrayA[i]<pivotValue) {
      swap<T>(arrayA[i],arrayA[storeIndex]);
      if (arrayB!=NULL)
	swap<U>(arrayB[i],arrayB[storeIndex]);
      if (arrayC!=NULL)
	swapPtr<V>(arrayC[i],arrayC[storeIndex]);
      storeIndex++;

    }
  }
  swap<T>(arrayA[storeIndex],arrayA[right]);
  if (arrayB!=NULL)
    swap<U>(arrayB[storeIndex],arrayB[right]);
  if (arrayC!=NULL)
    swapPtr<V>(arrayC[storeIndex],arrayC[right]);
    
  return storeIndex;
}

template <class T>
void maxOf1DArray(T const * const V, int const length, int &maxIndex, 
		  T &maxVal)
{
  
  /* Initializations. */
  maxIndex = 0;
  maxVal = V[0];
  
  /* Find the maximum value and its index. */
  for (int i = 1; i < length; i++)
  {
    if (V[i] > maxVal)
    {
      maxVal = V[i];
      maxIndex = i;
    }
  }
  
}

template <class T>
void minOf1DArray(T const * const V, int const length, T &minVal, 
		  int &minIndex)
{

  /* Initializations. */
  minIndex = 0;
  minVal = V[0];

  /* Find the minimum value and its index. */
  for (int i = 1; i < length; i++)
  {
    if (V[i] < minVal)
    {
      minVal = V[i];
      minIndex = i;
    }
  }
  
}

template <class T>
void maxOfAbs2DArray(T const * const * const V, int const nLen1, int const nLen2,
		     int &maxIndex1, int &maxIndex2, T &maxVal)
{
  
  /* Initializations. */
  maxIndex1 = 0;
  maxIndex2 = 0;
  maxVal = T(fabs(V[0][0]));
  
  /* Find the maximum value and its index. */
  for (int i=0; i<nLen1; i++){
    for (int j=0; j<nLen2; j++){
      if (T(fabs(V[i][j]))>maxVal){
	maxVal = T(fabs(V[i][j]));
	maxIndex1 = i;
	maxIndex2 = j;
      }
    }
  }

}

template <class T>
void swapPtr(T* &a, T* &b)
{

  T* c(a);
  a = b; 
  b = c;
}

template <class T>
string num2string(T const input)
{

  ostringstream oss;
  oss << input;
  string output = oss.str();
  oss.str("");
  
  return output;

}

template <class T, class U>
void minOfVector(const T &V, int &minIndex, U &minVal)
{
  
  /* Initializations. */
  minIndex = 0;
  minVal = V[0];
  
  /* Find the minimum value and its index. */
  for (unsigned int i=1; i<V.size(); i++){
    if (V[i]<minVal){
      minVal = V[i];
      minIndex = i;
    }
  }
}

template <class T>
int sgn(T const d)
{
  return d < -T(EPS) ? -1 : d > T(EPS);
}



/* LAPACK function protoypes. */
extern "C"
{
  
  /*! \fn void dstevx_(char &, char &, int &, double*, double*, double
    &, double &, int &, int &, double &, int &, double*, double*, int
    &, double*, int*, int*, int &);

    \brief LAPACK function. Please see
    http://netlib.org/lapack/double/dstevx.f for inputs and outputs.
  */
  void dstevx_(char &, char &, int &, double*, double*, double &, double &, 
	       int &, int &, double &, int &, double*, double*, int &, 
	       double*, int*, int*, int &);
  
  /*! \fn double dlamch_(char &);
    
    \brief LAPACK function. Please see
    http://www.netlib.org/lapack/util/dlamch.f for inputs and outputs.
  */
  double dlamch_(char &);

  /*! \fn void dgemv_(char &, int &, int &, double &, double*, int &,
    double*, int &, double &, double*, int &);
    
    \brief BLAS function. Please see
    http://www.netlib.org/blas/dgemv.f for inputs and outputs.
  */
  void dgemv_(char &, int &, int &, double &, double*, int &, double*, 
	      int &, double &, double*, int &);

  /*! \fn void dpotrf_(char &, int &, double*, int &, int &);
    
    \brief LAPACK function. Please see
    http://netlib.org/lapack/double/dpotrf.f for inputs and outputs.
  */
  void dpotrf_(char &, int &, double*, int &, int &);

  /*! \fn double dlange_(char &, int &, int &, double*, int &,
    double*);
    
    \brief LAPACK function. Please see 
    http://www.netlib.org/lapack/explore-html/a00725_source.html
    for inputs and outputs.

  */
  double dlange_(char &, int &, int &, double*, int &, double*);
  
  /*! \fn void dgetrf_(int &, int &, double*, int &, int*, int &);
    
    \brief LAPACK function. Please see
    http://www.netlib.org/lapack/double/dgetrf.f for inputs and
    outputs.

  */
  void dgetrf_(int &, int &, double*, int &, int*, int &);

  /*! \fn void dgecon_(char &, int &, double*, int &, double &, double
    &, double*, int*, int &);
    
    \brief LAPACK function. Please see
    http://www.netlib.org/lapack/double/dgecon.f for inputs and
    outputs.
  */
  void dgecon_(char &, int &, double*, int &, double &, double &, double*,
	       int*, int &);
  
  /*! \fn void dgbsv_(int &, int &, int &, int &, double*, int &, int*,
      double*, int &, int &);
    
    \brief LAPACK function. Please see
    http://netlib.org/lapack/double/dgbsv.f for inputs and outputs.
  */
  void dgbsv_(int &, int &, int &, int &, double*, int &, int*, double*, 
	      int &, int &);

  /*! \fn void dgbtrs_(char &, int &, int &, int &, int &, double*,
    int &, int*, double*, int &, int &);
    
    \brief LAPACK function. Please see
    http://netlib.org/lapack/double/dgbtrs.f for inputs and outputs.
  */
  void dgbtrs_(char &, int &, int &, int &, int &, double*, int &, int*, 
	       double*, int &, int &);

  /*! \fn void dgelss_(int &, int &, int &, double*, int &, double*,
    int &, double*, double &, int &, double*, int &, int &);
    
    \brief LAPACK function. Please see
    http://netlib.org/lapack/double/dgelss.f for inputs and outputs.
  */
  void dgelss_(int &, int &, int &, double*, int &, double*, 
	       int &, double*, double &, int &, double*, int &, int &);
  
}

#endif
