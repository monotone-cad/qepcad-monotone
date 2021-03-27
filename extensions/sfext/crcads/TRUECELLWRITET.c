/*======================================================================
                 TRUECELLWRITET(D)

Write all true cells as tarski formulas.

Inputs
  D : A CAD.

Side Effects
  Information for each true cell in the CAD is written to the output 
  stream.
======================================================================*/
#include "qepcad.h"

void QepcadCls::TRUECELLWRITET(Word D)
{
      Word T,F,c;

Step1: /* Get list of true cells, and write info for each element. */
      LISTOFCWTV(D,&T,&F);
      while (T != NIL) {
	ADV(T,&c,&T);
	CELLWRT(c); }

Return: /* Prepare to return. */
      return;
}
