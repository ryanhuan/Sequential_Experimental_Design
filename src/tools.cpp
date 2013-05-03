#include "tools.h"

/* Functions with templates are defined in header file. */

// void randMultivariateNormal(vector< vector<double> > const &sqrtCov, 
// 			    gsl_rng const * generator, vector<double> &sample)
// {
  
//   /* Initialization. */
//   unsigned int length(sqrtCov.size());
//   vector<double> indep1DGaussians;
//   indep1DGaussians.assign(length, 0.0);
  
//   /* Construct W, a vector of independent N(0,1) rv's. */
//   for (unsigned int i = 0; i < length; i++)
//     indep1DGaussians[i] = gsl_ran_gaussian(generator, 1.0);
  
//   /* Matrix-vector multiplication. */
//   mvProduct(sqrtCov, indep1DGaussians, sample);

// }

double vectorNorm(vector<double> const &vectorRef, int const normChoice)
{
  
  /* Initializations. */
  double normValue(0.0);
  unsigned int length(vectorRef.size());

  switch (normChoice)
  {
    
  case 0:
    /* Inifinity norm (does not return position of the maximum
     * element). */
    normValue = fabs(vectorRef[0]);
    for (unsigned int i = 1; i < length; i++)
      normValue = max<double>(normValue, fabs(vectorRef[i]));
    break;
    
  case 1:
    /* 1-norm. */
    for (unsigned int i = 0; i < length; i++)
      normValue += fabs(vectorRef[i]);
    break;
    
  case 2:
    /* 2-norm. */
    for (unsigned int i = 0; i < length; i++)
      normValue += pow(vectorRef[i], 2);
    normValue = sqrt(normValue);
    break;

  default:
    cout << "Error: vector norm choice " << normChoice
	 << " not available." << endl;
    exit(1);
  }
  
  return normValue;
  
}

int nCr(int const n, int const r)
{
  
  /* Initializations. */
  int largerDenom = max<int>(r, n - r);
  int smallerDenom = n - largerDenom;
  int result(1);

  /* Cancel the larger of the denominator terms. */
  for (int i = largerDenom; i < n; i++)
  {
    result *= i + 1;
    if (result < 0)
      return -1;
  }

  /* Final result. */
  result /= int(tgamma(smallerDenom + 1));
  
  return result;

}

// int binarySearch(double const query, double const * const refList, 
// 		 int const nLength, double const dCompTol)
// {
  
//   /* Initializations. */
//   int first(0), last(nLength-1);

//   while (first <= last)
//   {
//     int mid = (first + last) / 2;  /* Compute mid-point. */
//     if (query - refList[mid] > dCompTol)
//       first = mid + 1;    /* Repeat search in top half. */
//     else if (query - refList[mid] < -dCompTol) 
//       last = mid - 1;     /* Repeat search in bottom half. */
//     else
//       return mid;       /* Found it, return position. */
//   }
  
//   return -1;            /* Failed to find key. */
  
// }

// void threeTermLegendre(double const x, int const maxDegree, 
// 		       double* const P)
// {
  
//   /* Compute term via 3-term recursive formula. */
//   for (int i = 0; i < maxDegree + 1; i++)
//   {
//     if (i == 0)
//       P[0] = 1.0;
//     else if (i == 1)
//       P[1] = x;
//     else
//       P[i] = ((2.0 * double(i - 1) + 1.0) * x * P[i - 1] 
// 	      - double(i - 1) * P[i - 2]) / double(i);
//   }
  
// }

// void threeTermHermite(double const x, int const maxDegree, 
// 		      double* const P)
// {
  
//   /* Compute term via 3-term recursive formula. */
//   for (int i = 0; i < maxDegree + 1; i++)
//   {
//     if (i == 0)
//       P[0] = 1.0;
//     else if (i == 1)
//       P[1] = x;
//     else
//       P[i] = x * P[i - 1] - double(i - 1) * P[i - 2];
//   }

// }

// void threeTermLegendreDerivatives(double const x, int const maxDegree, 
// 				  double const * const P, double* const dP)
// {
  
//   /* Compute term via direct differentiation of the 3-term recursive
//    * formula. */
//   for (int i = 0; i < maxDegree + 1; i++)
//   {
//     if (i == 0)
//       dP[0] = 0.0;
//     else if (i == 1)
//       dP[1] = 1.0;
//     else
//       dP[i] = ((2.0 * double(i - 1) + 1.0) * P[i - 1] 
// 	       + (2.0 * double(i - 1) + 1.0) * x * dP[i - 1] 
// 	       - double(i - 1) * dP[i - 2]) / double(i);
//   }
  
// }

// void threeTermHermiteDerivatives(double const x, int const maxDegree, 
// 				 double const * const P, double* const dP)
// {
  
//   /* Compute term via direct differentiation of the 3-term recursive
//    * formula. */
//   for (int i = 0; i < maxDegree + 1; i++)
//   {
//     if (i == 0)
//       dP[0] = 0.0;
//     else if (i == 1)
//       dP[1] = 1.0;
//     else
//       dP[i] = double(i) * P[i - 1];
//   }
  
// }

// void eigSymTridiag(int const &nLength, double* const diagonal, 
// 		   double* const subdiagonal, double const &abstol, 
// 		   double* const eigenVals, double* const eigenVecs, 
// 		   int* const ifailOut, double* const work, int* const iwork)
// {
  
//   /* Memory allocation and initializations. */
//   char jobzIn('V'), rangeIn('A');
//   double vlIn(0.0), vuIn(0.0), abstolIn(abstol);
//   int nIn(nLength), ilIn(0), iuIn(0), ldzIn(nLength), mOut(0), infoOut(0);
  
//   /* Compute eigenvalue and eigenvectors. */
//   dstevx_(jobzIn, rangeIn, nIn, diagonal, subdiagonal, vlIn, vuIn, ilIn, iuIn,
// 	  abstolIn, mOut, eigenVals, eigenVecs, ldzIn, work, iwork, ifailOut, 
// 	  infoOut);
//   if (infoOut != 0)
//   {
//     cout << "Error: DSTEVX LAPACK function returned error. " << endl;
//     exit(1);
//   }
  
// }

// void mvProduct(int const nLength1, int const nLength2, double* const * const A, 
// 	       double* const xIn, double* const yOut)
// {
  
//   /* Initialization. */
//   char transIn('T');
//   int mIn(nLength1), nIn(nLength2), ldaIn(nLength1), incxIn(1), incyIn(1);
//   double alphaIn(1.0), betaIn(0.0);
  
//   /* Matrix-vector multiplication. */
//   dgemv_(transIn, mIn, nIn, alphaIn, A[0], ldaIn, xIn, incxIn, betaIn, 
// 	 yOut, incyIn);
  
// }

// void choleskyDecomp(int const nLength, double const * const * const A, 
// 		    double* const * const cholesky)
// {
  
//   /* Initializations. */
//   char uploIn('U');
//   int infoOut(0), nIn(nLength), ldaIn(nLength);
//   for (int i = 0; i < nLength; i++)
//   {
//     for (int j = 0; j < nLength; j++)
//     {
//       if (j <= i)
// 	cholesky[i][j] = A[i][j];
//       else
// 	cholesky[i][j] = 0.0;
//     }
//   }

//   /* Cholesky decomposition. */
//   dpotrf_(uploIn, nIn, cholesky[0], ldaIn, infoOut);
//   if (infoOut != 0)
//   {
//     cout << "Error: DPOTRF LAPACK function returned error. " << endl;
//     exit(1);
//   }
  
// }

// double invCondNumber(int const length, double* const * const A, 
// 		     int* const ipivOut, double* const work, int* const iwork)
// {
  
//   /* Initializations. */
//   char normIn('1');
//   double rcondOut(0.0);
//   int nIn(length), ldaIn(length), infoOut(0);

//   /* Compute norm of matrix. */
//   double anormIn = dlange_(normIn, nIn, nIn, A[0], ldaIn, work);
  
//   /* LU-factorize of matrix. */
//   dgetrf_(nIn, nIn, A[0], ldaIn, ipivOut, infoOut);
//   if (infoOut != 0)
//   {
//     cout << "Error: DGETRF LAPACK function returned error. " << endl;
//     exit(1);
//   }
  
//   /* Estimate inverse condition number. */
//   dgecon_(normIn, nIn, A[0], ldaIn, anormIn, rcondOut, work, iwork, 
// 	  infoOut);
//   if (infoOut != 0)
//   {
//     cout << "Error: DGECON LAPACK function returned error. " << endl;
//     exit(1);
//   }
  
//   return rcondOut;
  
// }

// void solveBandedLinearSystem(int const nLength, int const nSubdiagonals, 
// 			     int const nSuperdiagonals, double const * const * const A,
// 			     double* const * const ABIn, 
// 			     int* const ipivOut, double* const uIn)
// {
  
//   /* Initialization. */
//   int nIn(nLength), klIn(nSubdiagonals), kuIn(nSuperdiagonals), nrhsIn(1);
//   int ldabIn(2 * klIn + kuIn + 1), ldbIn(nLength), infoOut(0);
  
//   /* Create LAPACK-requested LHS matrix format. Note the transpose. */
//   for (int i = 0; i < nLength * (2 * klIn + kuIn + 1); i++)
//     ABIn[0][i] = 0.0;
//   for (int j = 0; j < nLength; j++)
//   {
//     for (int i = max(0, j - kuIn); i < min(nLength, j + 1 + klIn); i++)
//       ABIn[j][klIn + kuIn + i - j] = A[i][j];
//   }
  
//   /* Solve a banded linear system. */
//   dgbsv_(nIn, klIn, kuIn, nrhsIn, ABIn[0], ldabIn, ipivOut, uIn, ldbIn, 
// 	 infoOut);
//   if (infoOut != 0)
//   {
//     cout << "Error: DGBSV LAPACK function returned error code "
// 	 << infoOut << ". " << endl;
//     exit(1);
//   }
  
// }

// void solveBandedLinearSystemFromLU(int const nLength, int const nSubdiagonals, 
// 				   int const nSuperdiagonals, 
// 				   double * const * const AB,
// 				   int* const ipivIn, double* const bIn)
// {
  
//   /* Initialization. */
//   char trans('N');
//   int nIn(nLength), klIn(nSubdiagonals), kuIn(nSuperdiagonals), nrhsIn(1);
//   int ldabIn(2 * klIn + kuIn + 1), ldbIn(nLength), infoOut(0);
  
//   /* Solve a banded linear system from available LU. */
//   dgbtrs_(trans, nIn, klIn, kuIn, nrhsIn, AB[0], ldabIn, ipivIn, bIn, ldbIn, 
// 	  infoOut);
//   if (infoOut != 0)
//   {
//     cout << "Error: DGBTRS LAPACK function returned error code "
// 	 << infoOut << ". " << endl;
//     exit(1);
//   }
  
// }

/* Note: ATrans, B will be over-written. ATrans is also stored in the
 * transposed form of the actual algebraic form of the A matrix, due
 * to LAPACK needs. */
void linearLeastSquares(vector< vector<double> > &ATrans, 
			vector<double> &sOut, vector<double> &B, 
			vector<double> &work, int const workLength, 
			vector<double> &soln)
{

  /* Note that ATrans is the transposed matrix. */
  unsigned int nCols = ATrans.size();
  unsigned int nRows = ATrans[0].size();

  int mIn(nRows), nIn(nCols), nrhsIn(1), ldaIn(nRows), ldbIn(nRows);
  double rcondIn(-1.0);
  int rankOut(0), lwork(workLength), infoOut(0);

  static double** ATransDblPtr = new double*[nCols];
  static double* ATransPtr = new double[nCols * nRows];
  ATransDblPtr[0] = ATransPtr;
  for (unsigned int i = 1; i < nCols; i++)
    ATransDblPtr[i] = ATransDblPtr[i - 1] + nRows;
  for (unsigned int i = 0; i < nCols; i++)
  {
    for (unsigned int j = 0; j < nRows; j++)
      ATransDblPtr[i][j] = ATrans[i][j];
  }

  /* Solve linear least squares problem. */
  dgelss_(mIn, nIn, nrhsIn, ATransPtr, ldaIn, &(B[0]), ldbIn, &(sOut[0]), 
  	  rcondIn, rankOut, &(work[0]), lwork, infoOut);
  
  /* Extract the solution part from the B vector. */
  for (unsigned int i = 0; i < nCols; i++)
    soln[i] = B[i];

  if (infoOut != 0)
  {
    cout << "Error: DGELSS LAPACK function returned error code "
	 << infoOut << ". " << endl;
    exit(1);
  }
  
}
