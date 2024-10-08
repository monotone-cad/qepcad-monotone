/*======================================================================
                      CELLIPLLDWR(V,A,S,k)

Cell signs; Integral Polynomial, List of Lists, Distributed Write.

\Input
  \parm{V} is a non-null list of $r$ variables.
  \parm{A} is a list $A_1,\ldots,A_r)$ where each $A_i$ is a list of
           $i$--variate integral polynomials.
  \parm{S} is a list of signs corresponding to the polynomials
======================================================================*/
#include "qepcad.h"

void CELLIPLLDWR(Word V, Word A, Word S, Word k)
{
    Word A1, A11, i, P, xp, s;
    /* hide i,j,n,r; */

Step1: /* Write. */
    i = 0;
    while (A != NIL)
    {
        if (k-i <= 0) break;

        ADV(A,&A1,&A);
        // get signs for level i (they come in reverse level order)
        xp = LELTI(S,k-i);
        i++;

        SWRITE("Level "); IWRITE(i); SWRITE("\n");
        while (A1 != NIL && xp != NIL)
        {
            ADV(A1,&A11,&A1);
            ADV(xp,&s,&xp);

            if (LELTI(A11,PO_STATUS) == PO_REMOVE || LELTI(A11,PO_TYPE) == PO_POINT)
                continue;

            // write polynomial
            P = LELTI(A11, PO_POLY);

            SWRITE("  ");
            IPDWRITE(i,P,V);

            // write sign (< 0, = 0, > 0)
            switch (s)
            {
                case POSITIVE: SWRITE(" > "); break;
                case ZERO: SWRITE(" = "); break;
                case NEGATIVE: SWRITE(" < "); break;
                default: SWRITE(" ? "); break;
            }

            GWRITE(0);
            SWRITE("\n");
        }
    }

Return: /* Prepare for return. */
    return;
}
