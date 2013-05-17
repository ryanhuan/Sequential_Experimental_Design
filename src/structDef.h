/*! \file structDef.h 

  \brief Definition of miscellaneous structures.

*/
#ifndef _STRUCTDEF_H
#define _STRUCTDEF_H

#include <vector>

#include "inputParams.h"

/* Note we cannot #include "linearArch.h" since that would lead to a
 * recursive definition of header files. Instead, we will do a forward
 * declaration. */
class LinearArch;

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
  LinearArch* futureLinearArchClassPtr;
  
};

#endif
