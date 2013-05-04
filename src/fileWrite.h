/*! \file fileWrite.h 

  \brief File writing capabilities.

*/
#ifndef _FILEWRITE_H
#define _FILEWRITE_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "structDef.h"

using namespace std;

/*! \fn void write2DArray(const Controls&, const string, T const * const * const,
  const int, const int);

  \brief Writes a 2D array of data to file.
  
  \param primary Reference to primary controls.
  \param fName File name.
  \param input 2D array to be written to file.
  \param nLen1 Length of 1st dimension.
  \param nLen2 Length of 2nd dimension.
*/
template <class T>
void write2DArray(const Controls&, const string, T const * const * const,
		  const int, const int);



/* Template definitions. */

template <class T>
void write2DArray(const Controls &primary, const string fName, 
		  T const * const * const input, const int nLen1, const int nLen2)
{

  /* Open file. */
  ofstream f;

  /* Write as text file. */
  f.open(fName.c_str());

  /* Set output format. */
  f.setf(ios::scientific,ios::floatfield);
  f.precision(16);

  for (int i = 0; i < nLen1; i++)
  {
    for (int j = 0; j < nLen2; j++)
      f << input[i][j] << "  ";
    f << endl;
  }
  f << endl;

  /* Close file. */
  f.close();

}

#endif
