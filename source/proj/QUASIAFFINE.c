/*======================================================================
  QUASIAFFINE(A, r, *A_);

Adds first partial derivatives of projections onto coordinate axes, to result in 2d cells which are either increasing,
decreasing or constant with respect to coordinate axes.

\Input
  \parm{A} set of projection factors (A_1,...,A_r), each A_i is a list of polynomials in Z[x_1,...,x_i]
  \parm{r} positive integer, number of variables

Output
  \parm{*A}: set of projection factors, modified to include first partial derivatives of each element with respect to
  each vaciable.

SideEffect
  \parm{A} is also modified.

======================================================================*/
#include "qepcad.h"

// all nonconstant partial derivatives of polynomial P, in Z[x_1,...,x_k]
Word PARTIALS(Word k, Word P)
{
    Word L = IPALLPARTIALS(k, LELTI(P, PO_POLY), 1, 1);

    Word LL = NIL, D;
    while (L != NIL) {
        ADV(L, &D, &L);

        // skip constant derivatives
        if (IPCONST(k, D)) continue;

	   LL = COMP(MPOLY(D,NIL,LIST1(LIST3(PO_DER,0,P)),PO_POLY,PO_KEEP), LL);
    }

    return LL;
}

void QepcadCls::QUASIAFFINE(Word A, Word r, Word *A_)
{
    Word AA, A1, A11, k, L;

    // nothing to do
    if (r <= 1) {
        *A_ = A;

        return;
    }

    /* consider projection factors at each level */
    AA = A; k = 0;
    while (AA != NIL) {
        ADV(AA, &A1, &AA);
        ++k;

        // consider each polynomial
        L = NIL;
        while (A1 != NIL) {
            ADV(A1, &A11, &A1);

            // adding each partial derivative
            L = CONC(L, PARTIALS(k, A11));
        }

        // and factorise (same function as used with the input formula)
	    ADDPOLS(IPLFAC(k, L),k,LFS("D"), &A);
    }

    /* prepare for return */
    *A_ = A;
}

