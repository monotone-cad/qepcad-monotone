/*======================================================================
                    F -< SEMIMONOTONE(FF)

Adds extra polynomials to projection factor set to ensure semi-monotone cells will be produced
Does not recompute the cad

\Input
  \parm{A} is a list a_1, ..., a_r where a_ is a list of i-level projection factors
  \parm{AA} is a list a_1, ..., a_r where a_ is a list of i-level input polynomials
  \parm{D} is a CAD
  \parm{r} is the space in which D lives

======================================================================*/
#include "qepcad.h"

void QepcadCls::SEMIMONOTONE(Word A, Word AA, Word D, Word r)
{
    Word Ct, Cf, C, I, i ,j, A1, P, p, P1, p1, L;

Step2: /* calculate: looping through true cells with CELLDIM == 2 */
    LISTOFCWTV(D, &Ct, &Cf);
    while (Ct != NIL) {
        ADV(Ct, &C, &Ct);
        I = LEVELIDX(C);

        // choose an action based on dimension of cell
        if (LENGTH(I) > 2) {
            SWRITE("*** ERROR doesn't handle cells of dimension > 2");

            return;
        } else if (LENGTH(I) < 2) {
            printf("# cell dim < 2 - continue\n");

            continue;
        }

        /* We have a 2-dimensional cell */
        FIRST2(I, &i, &j);

        printf("# two dimensional index: %d %d\n", i, j);
        /* adding partials of level j polynomials with respect to variable i */
        L = NIL;
        A1 = LELTI(A, j);
        while (A1 != NIL) {
            ADV(A1, &P, &A1);
            p = LELTI(P, PO_POLY);
            p1 = IPDER(j, p, i);

            if (IPCONST(j, p1)) continue;

            P1 = MPOLY(p1, NIL, NIL, PO_OTHER, PO_KEEP);
            L = COMP(P1, L);
        }

        printf("# %d ", LENGTH(LELTI(A, i)));
        ADDPOLS(L, j, LFS("D"), &AA);
        A = APPEND(A, j, IPLFAC(j, L));
    }
}

