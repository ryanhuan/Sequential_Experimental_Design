#include "fileRead.h"

/* Functions with templates are defined in header file. */

void readScalarAsString(char* const val, char const * const data, int &pos)
{
  
  /* Initializations. */
  int ctrVal(0);
  
  /* Find starting point for numeric value. */
  while (data[pos] == ' ' || data[pos] == '=')
    pos++;
  
  /* Record numeric value. */
  while (data[pos] != ' ' && data[pos] != '\0')
  {
    val[ctrVal] = data[pos];
    ctrVal++;
    pos++;
  }
  val[ctrVal] = '\0';
  
}
