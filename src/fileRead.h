/*! \file fileRead.h 

  \brief File reading capabilities.
*/
#ifndef _FILEREAD_H
#define _FILEREAD_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>

#include "structDef.h"
#include "tools.h"

using namespace std;

/*! \fn void readFileControls(Controls &);

  \brief Reads the controls corrections from input file.
  
  \param primary Reference to controls structure of interest.
*/
void readFileControls(Controls &);

/*! \fn void readScalar(char* const, char const * const, int &);
  
  \brief Reads a single value (string) from file. The string can then
  be converted into integer, double, etc.
  
  \param val Pointer to store the value as a string.
  \param data Pointer to the full row data string (including variable
  name).
  \param pos Reference to current position in the data string.
*/
void readScalar(char* const, char const * const, int &);

/*! \fn void readVector(T* const, int const, int const, char const *
  const data, int &);

  \brief Reads a row vector (1D array) from file.

  \param var Pointer to store the vector.
  \param dataType Type of data T, options are INT_TYPE and DOUBLE_TYPE.
  \param length Length of the vector.
  \param data Pointer to the full row of data string (including variable name).
  \param pos Reference to current position in the data string.
*/
template <class T>
void readVector(T* const, int const, int const, char const * const data, int &);

/*! \def INT_TYPE

  \brief Indicator for int type.
*/
#define INT_TYPE 0

/*! \def DOUBLE_TYPE

  \brief Indicator for double type.
*/
#define DOUBLE_TYPE 1



/* Template definitions. */

template <class T>
void readVector(T* const array, int const dataType, int const length, 
		char const * const data, int &pos)
{

  /* Initializations. */
  char val[256];

  /* Extract each entry. */
  for (int i = 0; i < length; i++)
  {
    readScalar(val, data, pos);
    
    /* Assign values to arrayiable. */
    if (dataType == INT_TYPE)
      /* Integer type. */
      array[i] = atoi(val);
    else if (dataType == DOUBLE_TYPE)
      /* Doublereal type. */
      array[i] = atof(val);

    pos++;
  }
}

#endif
