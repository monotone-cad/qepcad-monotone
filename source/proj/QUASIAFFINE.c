/*======================================================================
                    A -< QUASIAFFINE(A, J, r; P, J)

Adding first derivatives of projections onto coordinate axes, to result in 2d cells which are either increasing,
decreasing or constant along coordinate axes.

\Input
  \parm{AA} set of input polynomials (I_1,...,i_r), each A_i is a list of polynomials in Z[x_1,...,x_i]
  \parm{r} number of variables

Output
  \parm{A}: modified projection factors

SideEffect
  \parm{AA} is modified.

======================================================================*/
#include "qepcad.h"

// return the index immediately after I in {0,1,...,M}^k with respect to lex order.
// if no greater index exists, returns NIL
// helper for LexNext, to ensure it stops at maximum index
Word LexNext(Word I, Word M)
{
    Word m, I1, J1;
    m = FIRST(I);
    I1 = RED(I);

    // increment m, if we can
    if (m < M) {
        return COMP(m + 1, I1);
    }

    // maximal index
    if (I1 == NIL) {
        return NIL;
    }

    // try to roll over
    J1 = LexNext(I1, M);

    if (J1 == NIL) return NIL; // fail, maximal index

    return COMP(0, J1);
}

// return list 0^k
Word ZEROS(Word k)
{
    Word L = NIL;
    while (k > 0) {
        L = COMP(0, L);
        --k;
    }

    return L;
}

// list of all partials, sufficient for smooth stratification
Word PARTIALS(Word r, Word d, Word L)
{
    Word LL, k;
    k = LENGTH(L);

    // indices
    Word I = ZEROS(r);
    while (I != NIL) {
        LWRITE(I); SWRITE("\n");

        I = LexNext(I, d);
    }   SWRITE("\n");

    return NIL;
}

void QepcadCls::QUASIAFFINE(Word A, Word r, Word *A_)
{
    Word AA, A1, A11, P, k, d, r1, L;

Step1: /* decide based on dimension */
    if (r <= 1) {
        // nothing to do

        goto Return;
    }

Step3: /* dim >= 2: all partials of input polynomials */
    /* levels. */
    AA = A, k = NIL, k = 0, r1 = r, d = 0;
    while (AA != NIL) {
        ADV(AA, &A1, &AA);
        ++k; --r1;

        /* polynomials. */
        while (A1 != NIL) {
            ADV(A1, &A11, &A1);
            P = LELTI(A11, PO_POLY);

            d = MAX(d, PDEG(P));
            L = COMP(PADDVS(LELTI(A11, PO_POLY), r1), L);
        } /* END polynomials. */
    } /* END level. */

    PARTIALS(r, d, L);
    // factorise and append -- same function as used on input formula
    //ADDPOLS(IPLFAC(k, L),k,LFS("D"), &A);

Return: /* prepare for return */
    // put proj fac in correct order
    *A_ = A;
}

