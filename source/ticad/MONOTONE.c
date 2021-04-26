/*======================================================================
                    F -< MONOTONE(FF)

Adds extra polynomials to projection factor set to ensure monotone cells will be produced
Does not recompute the cad

\Input
  \parm{A} is a list a_1, ..., a_r where a_ is a list of i-level projection factors
  \parm{AA} is a list a_1, ..., a_r where a_ is a list of i-level input polynomials
  \parm{D} is a CAD
  \parm{r} is the space in which D lives

======================================================================*/
#include "qepcad.h"

void QepcadCls::MONOTONE(Word A, Word AA, Word D, Word r)
{
    Word Ct, Cf, C, I, i ,j, A1, P, p, P1, p1, L;

Step2: /* calculate: looping through true cells */
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
        // TODO lagrange
    }
}


