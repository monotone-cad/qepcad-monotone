/*======================================================================
                    A -< QUASIAFFINE(A, J, r; P, J)

Adding first derivatives of projections onto coordinate axes, to result in 2d cells which are either increasing,
decreasing or constant along coordinate axes.

\Input
  \parm{AA} projection factor structure. ([A_1, ..., A_r], with A_i being i-level proj factors
  \parm{r} number of variables

Output
  \parm{A}: modified projection factors

SideEffect
  \parm{AA} is modified.

======================================================================*/
#include "qepcad.h"

// all (non-constant) partials of P (level k), factorised and formatted ready to append
Word PARTIALS(Word P, Word k);

void QepcadCls::QUASIAFFINE(Word A, Word r, Word* A_, Word *J_)
{
    Word AA, A1, A11, i, L;

Step1: /* decide based on dimension */
    if (r <= 1) {
        // nothing to do

        goto Return;
    }

    if (r > 3) { // TODO dim > 3
        SWRITE("Dimension > 3 not supported yet.\n");

        goto Return;
    }

Step3: /* dim >= 2: all partials of input polynomials */
    SWRITE("# adding derivs of projections.\n");

    /* levels. */
    AA = A; i = 0;
    while (AA != NIL) {
        ADV(AA, &A1, &AA);
        i++;

        /* polynomials. */
        L = NIL;
        while (A1 != NIL) {
            ADV(A1, &A11, &A1);
            // only work on input polynomials - skip projection factors
            if (LELTI(A11, PO_PARENT) != NIL) continue;

            L = CONC(L, PARTIALS(A11, i));
        } /* END polynomials. */

        // appending new level i polynomials
        ADDPOLS(L, i, LFS("D"), J_);
        A = APPEND(A, i, IPLFAC(i, L));
    } /* END level. */

Return: /* prepare for return */
    printf("> %d, %d\n", LENGTH(FIRST(*J_)), LENGTH(SECOND(*J_)));
    *A_ = A;

    GVPF = A;
    GVPJ = *J_;
}

Word PARTIALS(Word P, Word k)
{
    Word L = NIL;
    Word Pp = IPALLPARTIALS(k, LELTI(P, PO_POLY), 1, 1);
    IPWRITE(k, LELTI(P, PO_POLY), LIST3(LFS("x"), LFS("y"), LFS("z")));SWRITE("\n");

    Word D = NIL;
    /* derivatives. */
    while (Pp != NIL) {
        ADV(Pp, &D, &Pp);

        if (IPCONST(k, D)) continue;

        IPWRITE(k, D, LIST3(LFS("x"), LFS("y"), LFS("z")));SWRITE("\n");
        L = COMP(MPOLY(D, NIL, NIL, PO_OTHER, PO_KEEP), L);
    } /* END derivatives. */

    return L;
}

