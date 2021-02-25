/*======================================================================
                    F -< SEMIMONOTONE(FF)

Adds extra polynomials to projection factor set to ensure semi-monotone cells will be produced

\Input
  \parm{A} is a list a_1, ..., a_r where a_ is a list of i-level projection factors
  \parm{D} is a CAD
  \parm{r} is the space in which D lives

  Output
  \parm{D} new cad with semi-monotone cells

======================================================================*/
#include "qepcad.h"

// list of levels at which the cell is one-dimensional
Word LEVELIDX(Word C);

// factorise and null parents of a lest of level k polynomials
Word FAC(Word L, Word k);

Word QepcadCls::SEMIMONOTONE(Word A, Word DD, Word r)
{
    Word D, Ct, Cf, C, I, i ,j, A1, P, p, P1, p1, L, L1;

Step1: /* Initialise */
    D = DD;

Step2: /* calculate: looping through true cells */
    LISTOFCWTV(D, &Ct, &Cf);
    while (Ct != NIL) {
        ADV(Ct, &C, &Ct);
        I = LEVELIDX(C);
        if (LENGTH(I) > 2) {
            SWRITE("*** ERROR doesn't handle cells of dimension > 2");
            return D;
        } else if (LENGTH(I) < 2) {
            printf("# cell dim < 2 - continue\n");

            continue;
        }

        /* We have a 2-dimensional cell */
        i = LELTI(I, 1);
        j = LELTI(I, 2);

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
        A = APPEND(A, j, FAC(L, j));
        printf("%d\n", LENGTH(LELTI(A, i)));

    }

Return: /* returning */
    return D;
}

Word LEVELIDX(Word C)
{
    Word I, k, j, L;

    L = NIL;
    k = 0;
    I = LELTI(C, INDX);
    while (I != NIL) {
        ADV(I, &j, &I);
        k++;

        if (ODD(j)) {
            L = COMP(k, L);
        }
    }

    return INV(L);
}

Word FAC(Word L, Word k)
{
    Word LL, P;

    L = IPLFAC(k, L);

    LL = L;
    while (L != NIL) {
        ADV(L, &P, &L);
        SLELTI(P, PO_PARENT, NIL);
    }

    return LL;
}

