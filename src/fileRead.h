/*! \file fileRead.h 

  \brief Fundamental file reading capabilities.
*/
#ifndef _FILEREAD_H
#define _FILEREAD_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <vector>

using namespace std;

/*! \fn void readScalarAsString(char* const, char const * const, int
    &);
  
  \brief Reads a single value (string) from file. The string can then
  be converted into integer, double, etc.
  
  \param val Pointer to store the value as a string.
  \param data Pointer to the full row data string (including variable
  name).
  \param pos Reference to current position in the data string.
*/
void readScalarAsString(char* const, char const * const, int &);

/*! \fn void readVector(vector &vals, int const, char const * const,
 *  int &);

  \brief Reads a (row) vector from file.

  \param vals Reference to the vector.
  \param dataType Type of data T, options are INT_TYPE and DOUBLE_TYPE.
  \param data Pointer to the full row of data string (including variable name).
  \param pos Reference to current position in the data string.
*/
template<class T>
void readVector(vector<T> &vals, int const, char const * const, int &);

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
void readVector(vector<T> &vals, int const dataType, char const * const data, 
		int &pos)
{

  /* Initializations. */
  char val[256];

  /* Extract each entry. */
  for (unsigned int i = 0; i < vals.size(); i++)
  {
    readScalarAsString(val, data, pos);
    
    /* Assign values to vector. */
    if (dataType == INT_TYPE)
      /* Integer type. */
      vals[i] = atoi(val);
    else if (dataType == DOUBLE_TYPE)
      /* Double type. */
      vals[i] = atof(val);

    pos++;
  }
  
}

#endif
