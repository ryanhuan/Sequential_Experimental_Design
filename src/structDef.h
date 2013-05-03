/*! \file structDef.h 

  \brief Definition of miscellaneous structures.

*/
#ifndef _STRUCTDEF_H
#define _STRUCTDEF_H

#include <vector>

#include "inputParams.h"

using namespace std;

/*! \struct GenericInputs

  \brief Struct for general input components. 

  Note that only shallow copies are used. This struct should only
  contain pointers and not actual arrays themselves.
*/
struct GenericInputs
{
  
  /* Function pointers. */
  void (*systemEqnPtr) (InputParams const &, vector<double> const &, 
			vector<double> const &, vector<double> const &,
			vector<double> &);
  double (*stageFcnPtr) (InputParams const &, vector<double> const &, 
  			 vector<double> const &, vector<double> const &);
  double (*futureFcnPtr) (InputParams const &, vector<double> const &);
  
};

#endif
