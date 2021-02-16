/*======================================================================
                    F -< MONOTONE(FF)

Turns a quasi-affine CAD into a monotone one.

\Input
  \parm{DD} is a CAD (which has quasi-affine cells)

  Output
  \parm{D} is the CAD modified so all cells are monotone

======================================================================*/
#include "qepcad.h"

Word QepcadCls::MONOTONE(Word DD)
{

Step1: /* Initialise */

    SWRITE("In monotone.");

Return:
    return DD;
}
