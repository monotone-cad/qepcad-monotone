/*======================================================================
                      CELLIPLLDWR(V,A,x)

Cell signs; Integral Polynomial, List of Lists, Distributed Write.

\Input
  \parm{V} is a non-null list of $r$ variables.
  \parm{A} is a list $A_1,\ldots,A_r)$ where each $A_i$ is a list of
           $i$--variate integral polynomials.
  \parm{S} is a list of signs corresponding to the polynomials
======================================================================*/
#include "qepcad.h"

void CELLIPLLDWR(Word V, Word A, Word S)
{
       Word A1,A11,i,P,xp,s;
       /* hide i,j,n,r; */

Step1: /* Write. */
       i = 0;
       while (A != NIL)
         {
         ADV(A,&A1,&A);
	 i = i + 1;
         xp = LELTI(S,i);

         while (A1 != NIL)
           {
           ADV(A1,&A11,&A1);
           if (LELTI(A11,PO_STATUS) == PO_REMOVE || LELTI(A11,PO_TYPE) == PO_POINT)
	     continue;

	   P = LELTI(A11, PO_POLY);
           ADV(xp,&s,&xp);
	   IPDWRITE(i,P,V);
	   switch (s)
             {
             case POSITIVE: SWRITE(" > "); break;
             case ZERO: SWRITE(" = "); break;
             case NEGATIVE: SWRITE(" < "); break;
             case UNDET: SWRITE(" ? "); break;
             }
	   GWRITE(0);
           SWRITE("\n");
           }

         if (A != NIL)
           SWRITE("\n");
         }

Return: /* Prepare for return. */
       return;
}
