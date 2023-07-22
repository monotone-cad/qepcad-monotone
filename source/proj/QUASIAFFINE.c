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

Word PARTIALS(Word k, Word P) {
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

Step1: /* decide based on dimension */
    if (r <= 1) {
        // nothing to do

        goto Return;
    }

Step3: /* dim >= 2: all partials of input polynomials */
    /* levels. */
    AA = A; k = 0;
    L = NIL;
    while (AA != NIL) {
        ADV(AA, &A1, &AA);
        k++;

        /* polynomials. */
        while (A1 != NIL) {
            ADV(A1, &A11, &A1);

            L = CONC(L, PARTIALS(k, A11));
        } /* END polynomials. */

        // factorise and append -- same function as used on input formula
	    ADDPOLS(IPLFAC(k, L),k,LFS("D"), &A);
    } /* END level. */

Return: /* prepare for return */
    // put proj fac in correct order
    *A_ = A;
}

