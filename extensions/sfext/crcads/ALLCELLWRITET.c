/*======================================================================
                 ALLECELLWRITET(D)

Write all cells as tarski formulas.

Inputs
  D : A CAD.

Side Effects
  Information for each cell in the CAD is written to the output
  stream.
======================================================================*/
#include "qepcad.h"

void QepcadCls::ALLCELLWRITET(Word D)
{
      Word T,F,c;

Step1: /* Get list of true and false cells, and write info for each element. */
      LISTOFCWTV(D,&T,&F);
      SWRITE("TRUE CELLS\n\n");
      while (T != NIL) {
        ADV(T,&c,&T);
	    CELLWRT(c);
      }
      SWRITE("FALSE  CELLS\n\n");
      while (F != NIL) {
        ADV(F,&c,&F);
	    CELLWRT(c);
      }

Return: /* Prepare to return. */
      return;
}
