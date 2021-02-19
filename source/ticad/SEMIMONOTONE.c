/*======================================================================
                    F -< SEMIMONOTONE(FF)

Adds extra polynomials to projection factor set to ensure semi-monotone cells will be produced

\Input
  \parm{A} is a list a_1, ..., a_r where a_ is a list of i-level projection factors
  \parm{r} is the space in which D lives

  Output
  \parm{P} is a list p_1, ..., p_r where p_i is the list a_i with partial derivatives with respect to x_1 of a_i added

======================================================================*/
#include "qepcad.h"

Word QepcadCls::SEMIMONOTONE(Word A, Word r)
{
    Word P, L, L2, A1, p, d, W;

    if (r != 3) return A;
    // TODO handle r > 3
    printf("in semimonotone!!");fflush(0);

Step1: /* add derivatives of projection factors */
    P = A;
    L2 = NIL;
    int k = 0;
    while (A != NIL)  {
        ADV(A, &L, &A);
        k++;
        // skip level r and level 1
        if (k == 1 || k == r) continue;

        while (L != NIL) {
            ADV(L, &A1, &L);

            p = LELTI(A1, PO_POLY);

            d = IPDER(k, p, 1);

            if (IPCONST(k, d)) continue;

            W = MPOLY(d, NIL, NIL, PO_OTHER, PO_KEEP);

            IPDWRITE(k, d, GVVL); SWRITE("\n");
            L2 = COMP(W, L2);
        }

        printf("  %d\n", LENGTH(L2));
    }

    // using append from projection to handle the labels
    L2 = IPLFAC(r-1, L2);
    P = APPEND(P, r-1, L2);


Return: /* return */
    return P;
}
