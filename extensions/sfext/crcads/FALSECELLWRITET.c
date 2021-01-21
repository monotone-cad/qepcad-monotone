/*======================================================================
                 FALSECELLWRITET(D)

Write all false cells as tarski formulas.

Inputs
  D : A CAD.

Side Effects
  Information for each false cell in the CAD is written to the output 
  stream.
======================================================================*/
#include "qepcad.h"

void QepcadCls::FALSECELLWRITET(Word D)
{
      Word T,F,c;

Step1: /* Get list of true cells, and write info for each element. */
      LISTOFCWTV(D,&T,&F);
      while (F != NIL) {
	ADV(F,&c,&F);
	CELLWRT(c); }

Return: /* Prepare to return. */
      return;
}
